// Distributed output path for EventGraph: each dumpi_to_graph rank writes its
// own partial view of the event graph as JSON Lines partitions instead of
// sending everything to rank 0 (see EventGraph::merge). The partitions are
// assembled into event_graph.graphml by spark/jsonl_to_graphml.py.
//
// Layout written under <trace_dir>/event_graph_jsonl/:
//
//   _metadata.json            written by rank 0: format version, attribute
//                             schema and one entry per partition (vertex-id
//                             range, edge counts, process ids it holds)
//   vertices/part-NNNNN.jsonl one vertex per line, NNNNN = dumpi_to_graph rank
//   edges/part-NNNNN.jsonl    one edge per line
//
// Vertex line:
//   {"id":12,"event_type":"recv","logical_time":4,"wall_time":7.8807e-05,
//    "process_id":0,"mpi_function":"MPI_Recv"[,"callstack":"..."]
//    [,"PAPI_ctrs":"..."]}
//
// Edge line:
//   {"src":3,"dst":4,"order":0,"part":0,"seq":3}
//   order: 0 = program order, 1 = message order (includes collectives)
//   part/seq: the writing rank and the edge's position in that rank's list.
//   (order, part, seq) reproduces the exact edge order merge() hands to igraph,
//   which is what makes the assembled GraphML byte-identical.
//
// Values mirror exactly what merge() would put into the igraph object,
// including its defaults for missing map entries (0 / empty string).
// Non-finite doubles are written as the strings "nan", "inf" and "-inf".
// Strings are cut at the first NUL byte (igraph stores C strings) and bytes
// >= 0x80 are passed through untouched, so the reader must not re-encode them.

#include "EventGraph.hpp"

#include <mpi.h>

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cinttypes>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <malloc.h>     // malloc_trim
#include <sys/stat.h>   // mkdir
#include <unistd.h>     // unlink

#include "boost/mpi.hpp"
#include "boost/serialization/vector.hpp"

#include "JsonlFormat.hpp"

using namespace d2g_jsonl;

namespace {

const int JSONL_FORMAT_VERSION = 1;

// Same mapping merge() uses for the event_type vertex attribute
const std::unordered_map<uint8_t, std::string>& event_type_names()
{
  static const std::unordered_map<uint8_t, std::string> names =
  {
    {0, "send"}, {1, "recv"}, {2, "init"}, {3, "finalize"}, {4, "barrier"},
    {5, "reduce"}, {6, "allreduce"}, {7, "bcast"}, {8, "alltoall"},
    {9, "gather"}, {10, "scatter"}
  };
  return names;
}

void make_dir( const std::string& path )
{
  if ( mkdir( path.c_str(), 0775 ) != 0 && errno != EEXIST ) {
    std::stringstream ss;
    ss << "Could not create directory " << path << ": " << std::strerror( errno );
    throw std::runtime_error( ss.str() );
  }
}

std::string part_name( int rank )
{
  char buf[32];
  std::snprintf( buf, sizeof(buf), "part-%05d.jsonl", rank );
  return std::string( buf );
}

// Writes to <path>.tmp and renames on close so readers never see a partial file
class JsonlWriter
{
public:
  explicit JsonlWriter( const std::string& path )
    : path_( path ), tmp_path_( path + ".tmp" )
  {
    out_.open( tmp_path_.c_str(), std::ios::out | std::ios::trunc | std::ios::binary );
    if ( !out_.is_open() ) {
      throw std::runtime_error( "Could not open " + tmp_path_ + " for writing" );
    }
    buffer_.reserve( 1 << 20 );
  }

  void write_line( const std::string& line )
  {
    buffer_ += line;
    buffer_.push_back( '\n' );
    if ( buffer_.size() >= ( 1 << 20 ) ) {
      flush();
    }
  }

  void close()
  {
    flush();
    out_.close();
    if ( out_.fail() ) {
      throw std::runtime_error( "Failed writing " + tmp_path_ );
    }
    if ( std::rename( tmp_path_.c_str(), path_.c_str() ) != 0 ) {
      throw std::runtime_error( "Could not rename " + tmp_path_ + " to " + path_ );
    }
  }

private:
  void flush()
  {
    out_.write( buffer_.data(), buffer_.size() );
    buffer_.clear();
  }

  std::string path_;
  std::string tmp_path_;
  std::ofstream out_;
  std::string buffer_;
};

} // namespace

void EventGraph::write_jsonl()
{
  boost::mpi::communicator world;
  const int rank = world.rank();
  const int n_parts = world.size();
  const bool has_csmpi = this->config.has_csmpi();
  const bool has_papi = this->config.get_papi_flag();

  const std::string trace_dir = this->config.get_trace_dirs()[0];
  const std::string out_dir = trace_dir + "/event_graph_jsonl";
  const std::string vertex_dir = out_dir + "/vertices";
  const std::string edge_dir = out_dir + "/edges";
  const std::string metadata_path = out_dir + "/_metadata.json";

  // Rank 0 prepares the directory layout and drops the index from any earlier
  // run, so a reader can never pair stale partitions with a new index
  if ( rank == 0 ) {
    make_dir( out_dir );
    make_dir( vertex_dir );
    make_dir( edge_dir );
    unlink( metadata_path.c_str() );
  }
  world.barrier();

  // ---- Vertices ------------------------------------------------------------
  const auto& type_names = event_type_names();
  unsigned long long n_vertices = 0;
  unsigned long long min_vertex_id = 0;
  unsigned long long max_vertex_id = 0;
  std::set<int> process_ids;
  {
    JsonlWriter writer( vertex_dir + "/" + part_name( rank ) );
    std::string line;
    for ( size_t vid : this->vertex_ids ) {
      // Same lookups (and same defaults for missing entries) as merge()
      std::string event_type_str;
      auto et = this->vertex_id_to_event_type.find( vid );
      if ( et != this->vertex_id_to_event_type.end() ) {
        auto name = type_names.find( et->second );
        if ( name != type_names.end() ) {
          event_type_str = name->second;
        }
      }
      auto lts_it = this->logical_timestamps.find( vid );
      size_t lts = lts_it != this->logical_timestamps.end() ? lts_it->second : 0;
      auto wt_it = this->vertex_id_to_wall_time.find( vid );
      double wtime = wt_it != this->vertex_id_to_wall_time.end() ? wt_it->second : 0.0;
      auto pid_it = this->vertex_id_to_pid.find( vid );
      int pid = pid_it != this->vertex_id_to_pid.end() ? pid_it->second : 0;
      auto fn_it = this->vertex_id_to_fn_call.find( vid );
      std::string mpi_fn = fn_it != this->vertex_id_to_fn_call.end() ? fn_it->second.first : std::string();

      line.clear();
      line += "{\"id\":";
      append_uint( line, vid );
      line += ",\"event_type\":";
      append_json_string( line, event_type_str );
      line += ",\"logical_time\":";
      append_uint( line, lts );
      line += ",\"wall_time\":";
      append_json_double( line, wtime );
      line += ",\"process_id\":";
      append_int( line, pid );
      line += ",\"mpi_function\":";
      append_json_string( line, mpi_fn );
      if ( has_csmpi ) {
        auto cs_it = this->vertex_id_to_callstack.find( vid );
        line += ",\"callstack\":";
        append_json_string( line, cs_it != this->vertex_id_to_callstack.end() ? cs_it->second : std::string() );
      }
      if ( has_papi ) {
        auto papi_it = this->vertex_id_to_papi.find( vid );
        line += ",\"PAPI_ctrs\":";
        append_json_string( line, papi_it != this->vertex_id_to_papi.end() ? papi_it->second : std::string() );
      }
      line += "}";
      writer.write_line( line );

      if ( n_vertices == 0 || vid < min_vertex_id ) { min_vertex_id = vid; }
      if ( n_vertices == 0 || vid > max_vertex_id ) { max_vertex_id = vid; }
      ++n_vertices;
      process_ids.insert( pid );
    }
    writer.close();
  }
  // Vertex data is no longer needed on this rank
  std::vector<size_t>().swap( this->vertex_ids );
  std::unordered_map<size_t,size_t>().swap( this->logical_timestamps );
  std::unordered_map<size_t,uint8_t>().swap( this->vertex_id_to_event_type );
  std::unordered_map<size_t,double>().swap( this->vertex_id_to_wall_time );
  std::unordered_map<size_t,int>().swap( this->vertex_id_to_pid );
  std::unordered_map<size_t,std::pair<std::string,size_t>>().swap( this->vertex_id_to_fn_call );
  std::unordered_map<size_t,std::string>().swap( this->vertex_id_to_callstack );
  std::unordered_map<size_t,std::string>().swap( this->vertex_id_to_papi );

  // ---- Edges ---------------------------------------------------------------
  // merge() hands igraph every rank's program order edges (rank 0 first),
  // then every rank's message order edges, so (order, part, seq) is the key
  // that reproduces igraph's edge ids
  const unsigned long long n_program_edges = this->program_order_edges.size();
  const unsigned long long n_message_edges = this->message_order_edges.size();
  {
    JsonlWriter writer( edge_dir + "/" + part_name( rank ) );
    std::string line;
    const std::vector<std::pair<size_t,size_t>>* lists[2] =
      { &this->program_order_edges, &this->message_order_edges };
    for ( int order = 0; order < 2; ++order ) {
      unsigned long long seq = 0;
      for ( const auto& edge : *lists[order] ) {
        line.clear();
        line += "{\"src\":";
        append_uint( line, edge.first );
        line += ",\"dst\":";
        append_uint( line, edge.second );
        line += ",\"order\":";
        append_int( line, order );
        line += ",\"part\":";
        append_int( line, rank );
        line += ",\"seq\":";
        append_uint( line, seq );
        line += "}";
        writer.write_line( line );
        ++seq;
      }
    }
    writer.close();
  }
  std::vector<std::pair<size_t,size_t>>().swap( this->program_order_edges );
  std::vector<std::pair<size_t,size_t>>().swap( this->message_order_edges );
  malloc_trim( 0 );

  // ---- Index -----------------------------------------------------------------
  // Only a handful of numbers per rank travel to rank 0
  std::vector<unsigned long long> stats = {
    n_vertices, min_vertex_id, max_vertex_id, n_program_edges, n_message_edges
  };
  std::vector<int> pids( process_ids.begin(), process_ids.end() );
  std::vector<std::vector<unsigned long long>> all_stats;
  std::vector<std::vector<int>> all_pids;
  boost::mpi::gather( world, stats, all_stats, 0 );
  boost::mpi::gather( world, pids, all_pids, 0 );

  if ( rank == 0 ) {
    unsigned long long total_vertices = 0;
    unsigned long long total_edges = 0;
    for ( const auto& s : all_stats ) {
      total_vertices += s[0];
      total_edges += s[3] + s[4];
    }

    std::string meta;
    meta += "{\n";
    meta += "  \"format\": \"dumpi_to_graph.event_graph.jsonl\",\n";
    meta += "  \"version\": "; append_int( meta, JSONL_FORMAT_VERSION ); meta += ",\n";
    meta += "  \"directed\": true,\n";
    meta += "  \"n_partitions\": "; append_int( meta, n_parts ); meta += ",\n";
    meta += "  \"n_vertices\": "; append_uint( meta, total_vertices ); meta += ",\n";
    meta += "  \"n_edges\": "; append_uint( meta, total_edges ); meta += ",\n";
    // Attribute order matches the order merge() creates them in igraph, which
    // is the order of the <key> elements in the GraphML file
    meta += "  \"vertex_attributes\": [\n";
    meta += "    {\"name\": \"event_type\", \"type\": \"string\"},\n";
    meta += "    {\"name\": \"logical_time\", \"type\": \"double\"},\n";
    meta += "    {\"name\": \"wall_time\", \"type\": \"double\"},\n";
    meta += "    {\"name\": \"process_id\", \"type\": \"double\"},\n";
    meta += "    {\"name\": \"mpi_function\", \"type\": \"string\"}";
    if ( has_csmpi ) {
      meta += ",\n    {\"name\": \"callstack\", \"type\": \"string\"}";
    }
    if ( has_papi ) {
      meta += ",\n    {\"name\": \"PAPI_ctrs\", \"type\": \"string\"}";
    }
    meta += "\n  ],\n";
    meta += "  \"edge_orders\": [\"program\", \"message\"],\n";
    meta += "  \"partitions\": [\n";
    for ( int p = 0; p < n_parts; ++p ) {
      const auto& s = all_stats[p];
      meta += "    {\"partition\": "; append_int( meta, p );
      meta += ", \"vertices\": \"vertices/" + part_name( p ) + "\"";
      meta += ", \"edges\": \"edges/" + part_name( p ) + "\"";
      meta += ", \"n_vertices\": "; append_uint( meta, s[0] );
      if ( s[0] > 0 ) {
        meta += ", \"min_vertex_id\": "; append_uint( meta, s[1] );
        meta += ", \"max_vertex_id\": "; append_uint( meta, s[2] );
      }
      meta += ", \"n_program_order_edges\": "; append_uint( meta, s[3] );
      meta += ", \"n_message_order_edges\": "; append_uint( meta, s[4] );
      meta += ", \"process_ids\": [";
      for ( size_t i = 0; i < all_pids[p].size(); ++i ) {
        if ( i ) { meta += ", "; }
        append_int( meta, all_pids[p][i] );
      }
      meta += "]}";
      meta += ( p + 1 < n_parts ) ? ",\n" : "\n";
    }
    meta += "  ]\n";
    meta += "}\n";

    JsonlWriter writer( metadata_path );
    writer.write_line( meta.substr( 0, meta.size() - 1 ) );
    writer.close();
    std::cout << "Wrote event graph JSON Lines partitions to: " << out_dir
              << " (" << n_parts << " partitions, " << total_vertices
              << " vertices, " << total_edges << " edges)" << std::endl;
  }
  world.barrier();
}
