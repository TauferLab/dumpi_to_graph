// Checks that the JSON Lines partitions + spark/jsonl_to_graphml.py reproduce
// igraph's GraphML output byte for byte on awkward values: XML special
// characters, control characters, non-UTF-8 bytes, embedded NUL, NaN/inf,
// tiny/huge doubles, large logical times, and the callstack/PAPI attributes.
//
// The program builds the same igraph object EventGraph::merge() would build
// and writes it with igraph_write_graph_graphml (expected.graphml), and writes
// the same data as two JSONL partitions with the helpers dumpi_to_graph uses.
// test/check_jsonl_format.sh builds it, runs the Spark job, and compares.
//
// Usage: jsonl_format_test <out_dir>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <sys/stat.h>

#include "igraph/igraph.h"

#include "JsonlFormat.hpp"

using namespace d2g_jsonl;

struct Vertex {
  std::string event_type;
  size_t lts;
  double wtime;
  int pid;
  std::string fn;
  std::string callstack;
  std::string papi;
};

int main( int argc, char** argv )
{
  if ( argc != 2 ) {
    std::fprintf( stderr, "usage: %s <out_dir>\n", argv[0] );
    return 2;
  }
  const std::string out = argv[1];
  mkdir( out.c_str(), 0775 );
  const std::string jdir = out + "/event_graph_jsonl";
  mkdir( jdir.c_str(), 0775 );
  mkdir( ( jdir + "/vertices" ).c_str(), 0775 );
  mkdir( ( jdir + "/edges" ).c_str(), 0775 );

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<Vertex> v = {
    { "init", 0, 0.0, 0, "MPI_Init", "main <a&b> \"q\" 'apos'", "1,2,3" },
    { "send", 1, 1.23456789e-05, 0, "MPI_Send", "tab\there\nnewline\r\\back", "" },
    { "recv", 9007199254740993ULL, 123456789.0, 0, "MPI_Recv", "caf\xc3\xa9 \xff\xfe raw", "x" },
    { "finalize", 3, nan, 0, "MPI_Finalize", std::string( "cut\0here", 8 ), "y" },
    { "init", 0, inf, 1, "MPI_Init", "\x01\x1f ctl", "" },
    { "barrier", 5, -inf, 1, "MPI_Barrier", "", "PAPI_TOT_INS=42" },
    { "", 6, 1e-300, 1, "", "&&&&", "<>" },
    { "finalize", 18446744073709551615ULL, 0.1, 1, "MPI_Finalize", "end", "end" },
  };
  const int n = v.size();
  // partition 0 owns vertices 0..3 (pid 0), partition 1 owns 4..7 (pid 1)
  std::vector<std::vector<std::pair<size_t,size_t>>> program = {
    { {0,1}, {1,2}, {2,3} }, { {4,5}, {5,6}, {6,7} } };
  std::vector<std::vector<std::pair<size_t,size_t>>> message = {
    { {1,6}, {2,5} }, { {5,2}, {4,0}, {7,3} } };

  // ---- igraph, exactly like EventGraph::merge() ----
  igraph_i_set_attribute_table( &igraph_cattribute_table );
  igraph_t graph;
  igraph_empty( &graph, n, true );
  for ( int vid = 0; vid < n; ++vid ) {
    igraph_cattribute_VAS_set( &graph, "event_type", vid, v[vid].event_type.c_str() );
    igraph_cattribute_VAN_set( &graph, "logical_time", vid, v[vid].lts );
    igraph_cattribute_VAN_set( &graph, "wall_time", vid, v[vid].wtime );
    igraph_cattribute_VAN_set( &graph, "process_id", vid, v[vid].pid );
    igraph_cattribute_VAS_set( &graph, "mpi_function", vid, v[vid].fn.c_str() );
    igraph_cattribute_VAS_set( &graph, "callstack", vid, v[vid].callstack.c_str() );
    igraph_cattribute_VAS_set( &graph, "PAPI_ctrs", vid, v[vid].papi.c_str() );
  }
  std::vector<std::pair<size_t,size_t>> all;
  for ( auto& p : program ) all.insert( all.end(), p.begin(), p.end() );
  for ( auto& m : message ) all.insert( all.end(), m.begin(), m.end() );
  igraph_vector_t edges;
  igraph_vector_init( &edges, 2 * all.size() );
  for ( size_t i = 0; i < all.size(); ++i ) {
    VECTOR(edges)[2*i] = all[i].first;
    VECTOR(edges)[2*i+1] = all[i].second;
  }
  igraph_add_edges( &graph, &edges, 0 );
  FILE* f = std::fopen( ( out + "/expected.graphml" ).c_str(), "w" );
  igraph_write_graph_graphml( &graph, f, false );
  std::fclose( f );

  // ---- JSONL partitions, with the helpers EventGraph::write_jsonl() uses ----
  for ( int part = 0; part < 2; ++part ) {
    char name[32];
    std::snprintf( name, sizeof(name), "part-%05d.jsonl", part );
    std::ofstream vf( jdir + "/vertices/" + name, std::ios::binary );
    for ( int vid = part * 4; vid < part * 4 + 4; ++vid ) {
      std::string line = "{\"id\":";
      append_uint( line, vid );
      line += ",\"event_type\":"; append_json_string( line, v[vid].event_type );
      line += ",\"logical_time\":"; append_uint( line, v[vid].lts );
      line += ",\"wall_time\":"; append_json_double( line, v[vid].wtime );
      line += ",\"process_id\":"; append_int( line, v[vid].pid );
      line += ",\"mpi_function\":"; append_json_string( line, v[vid].fn );
      line += ",\"callstack\":"; append_json_string( line, v[vid].callstack );
      line += ",\"PAPI_ctrs\":"; append_json_string( line, v[vid].papi );
      line += "}\n";
      vf << line;
    }
    std::ofstream ef( jdir + "/edges/" + name, std::ios::binary );
    const std::vector<std::pair<size_t,size_t>>* lists[2] = { &program[part], &message[part] };
    for ( int order = 0; order < 2; ++order ) {
      size_t seq = 0;
      for ( auto& e : *lists[order] ) {
        std::string line = "{\"src\":";
        append_uint( line, e.first );
        line += ",\"dst\":"; append_uint( line, e.second );
        line += ",\"order\":"; append_int( line, order );
        line += ",\"part\":"; append_int( line, part );
        line += ",\"seq\":"; append_uint( line, seq++ );
        line += "}\n";
        ef << line;
      }
    }
  }
  std::ofstream mf( jdir + "/_metadata.json" );
  mf << "{\"format\": \"dumpi_to_graph.event_graph.jsonl\", \"version\": 1, "
        "\"directed\": true, \"n_partitions\": 2, \"n_vertices\": " << n
     << ", \"n_edges\": " << all.size() << ", \"vertex_attributes\": ["
        "{\"name\": \"event_type\", \"type\": \"string\"}, "
        "{\"name\": \"logical_time\", \"type\": \"double\"}, "
        "{\"name\": \"wall_time\", \"type\": \"double\"}, "
        "{\"name\": \"process_id\", \"type\": \"double\"}, "
        "{\"name\": \"mpi_function\", \"type\": \"string\"}, "
        "{\"name\": \"callstack\", \"type\": \"string\"}, "
        "{\"name\": \"PAPI_ctrs\", \"type\": \"string\"}], "
        "\"partitions\": ["
        "{\"partition\": 0, \"vertices\": \"vertices/part-00000.jsonl\", \"edges\": \"edges/part-00000.jsonl\"}, "
        "{\"partition\": 1, \"vertices\": \"vertices/part-00001.jsonl\", \"edges\": \"edges/part-00001.jsonl\"}]}\n";
  return 0;
}
