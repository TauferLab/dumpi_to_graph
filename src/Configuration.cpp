#include <iostream> 
#include <algorithm>
#include <sstream> 
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <set>
#include <cassert>

#include "Configuration.hpp"
#include "Glob.hpp"
#include "EventGraphModel.hpp"

// Used to parse the configuration file 
#include "external/nlohmann/json.hpp"

// Boost 
#include "boost/mpi.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/serialization/vector.hpp" 
#include "boost/serialization/set.hpp" 
#include "boost/serialization/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/operations.hpp"

using json = nlohmann::json; 

/* Top-level function to 
 */
d2g_Configuration parse_args(int argc, char** argv) 
{
  // First parse the configuration file itself
  d2g_Configuration base_config = parse_config_file( argv[1] );
  // Next, parse all of the provided input directories of trace files
  std::vector<std::string> trace_dirs;
  for (int i=2; i<argc; ++i) {
    // Canonicalize path
    boost::filesystem::path trace_dir_path( argv[i] );
    boost::filesystem::path trace_dir_canonical_path = boost::filesystem::canonical( trace_dir_path );
    trace_dirs.push_back( trace_dir_canonical_path.string() );
  }
  base_config.set_trace_dirs( trace_dirs );
  return base_config;
}

d2g_Configuration parse_config_file( std::string config_file_path )
{
  // Parse the JSON configuration file
  json config_json;
  std::ifstream config_stream(config_file_path);
  config_stream >> config_json; 
  // Set which MPI functions we're modeling in the graph
  std::set<std::string> mpi_functions;
  for (auto elem : config_json["mpi_functions"]) {
    std::string function_name(elem);
    auto search = allowed_mpi_functions.find( function_name );
    if (search != allowed_mpi_functions.end()) {
      mpi_functions.insert( function_name );
    }
  }
  // Set which happens-before orders we're keeping track of
  std::set<std::string> happens_before_orders;
  for (auto elem : config_json["happens_before_orders"]) {
    std::string order(elem);
    auto search = allowed_happens_before_orders.find( order );
    if (search != allowed_happens_before_orders.end()) {
      happens_before_orders.insert( order );
    }
  }
  // Set which vertex labels we will assign
  std::set<std::string> vertex_labels;
  for (auto elem : config_json["vertex_labels"]) {
    std::string vertex_label(elem);
    auto search = allowed_vertex_labels.find( vertex_label );
    if (search != allowed_vertex_labels.end()) {
      vertex_labels.insert( vertex_label );
    }
  }
  // Set which edge labels we will assign
  std::set<std::string> edge_labels;
  for (auto elem : config_json["edge_labels"]) {
    std::string edge_label(elem);
    auto search = allowed_edge_labels.find( edge_label );
    if (search != allowed_edge_labels.end()) {
      edge_labels.insert( edge_label );
    }
  }
  // Construct the configuration 
  d2g_Configuration config( mpi_functions,
                            happens_before_orders,
                            vertex_labels,
                            edge_labels,
                            (bool)config_json["represent_unmatched_tests"],
                            (bool)config_json["condense_unmatched_tests"],
                            (bool)config_json["condense_matched_tests"] );
  return config; 
}

std::set<std::string> d2g_Configuration::get_mpi_functions() const
{
  return this->mpi_functions;
}

std::set<std::string> d2g_Configuration::get_happens_before_orders() const
{
  return this->happens_before_orders;
}

std::set<std::string> d2g_Configuration::get_vertex_labels() const
{
  return this->vertex_labels;
}

std::set<std::string> d2g_Configuration::get_edge_labels() const
{
  return this->edge_labels;
}

std::vector<std::string> d2g_Configuration::get_trace_dirs() const
{
  return this->trace_dirs;
}

std::unordered_map<std::string,std::vector<std::string> > d2g_Configuration::get_trace_files() const
{
  return this->trace_files;
}

bool d2g_Configuration::get_represent_unmatched_tests_flag() const
{
  return this->represent_unmatched_tests;
}

bool d2g_Configuration::get_condense_unmatched_tests_flag() const
{
  return this->condense_unmatched_tests;
}

bool d2g_Configuration::get_condense_matched_tests_flag() const
{
  return this->condense_matched_tests;
}

d2g_Configuration& d2g_Configuration::operator=(const d2g_Configuration& rhs)
{
  if (&rhs == this) {
    return *this;
  } 
  this->mpi_functions = rhs.get_mpi_functions();
  this->happens_before_orders = rhs.get_happens_before_orders();
  this->vertex_labels = rhs.get_vertex_labels();
  this->edge_labels = rhs.get_edge_labels();
  this->represent_unmatched_tests = rhs.get_represent_unmatched_tests_flag();
  this->condense_unmatched_tests = rhs.get_condense_unmatched_tests_flag();
  this->condense_matched_tests = rhs.get_condense_matched_tests_flag();
  this->trace_dirs = rhs.get_trace_dirs();
  this->trace_files = rhs.get_trace_files();
  return *this; 
}

void d2g_Configuration::set_trace_dirs(std::vector<std::string> trace_dirs) 
{
  this->trace_dirs = trace_dirs;
}

void broadcast_config( d2g_Configuration& config )
{
  boost::mpi::communicator world;
  int rank = world.rank();
  std::stringstream config_ss;
  boost::archive::text_oarchive config_archive { config_ss };
  if ( rank == 0 ) {
    config_archive << config;
  }
  std::string config_payload = config_ss.str();
  boost::mpi::broadcast( world, config_payload, 0 );
  if ( rank != 0 ) {
    std::istringstream config_iss ( config_payload );
    boost::archive::text_iarchive config_iarchive { config_iss };
    config_iarchive >> config;
  }
}

void d2g_Configuration::compute_trace_file_assignment()
{
  assert(this->trace_dirs.size() > 0);
  boost::mpi::communicator world;
  int n_procs = world.size();
  int rank = world.rank();
  for (auto trace_dir : this->trace_dirs) {
    std::string trace_file_pattern = trace_dir + "/dumpi-*.bin";
    std::vector<std::string> curr_trace_files = glob(trace_file_pattern);
    int n_trace_files = curr_trace_files.size();
    std::vector<std::string> assigned_trace_files;
    for ( int i=0; i<n_trace_files; ++i ) {
      if ( i % n_procs == rank ) {
        // Canonicalize path
        boost::filesystem::path trace_file_path( curr_trace_files[i] );
        boost::filesystem::path trace_file_canonical_path = boost::filesystem::canonical( trace_file_path );
        assigned_trace_files.push_back( trace_file_canonical_path.string() );
      }
    }
    this->trace_files.insert( { trace_dir, assigned_trace_files } );
  }
}

