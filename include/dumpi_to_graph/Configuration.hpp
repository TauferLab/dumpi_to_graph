#ifndef D2G_CONFIGURATION_H
#define D2G_CONFIGURATION_H

#include <iostream>
#include <string>
#include <vector> 
#include <set> 
#include <unordered_map>

// Boost
#include "boost/serialization/access.hpp"
#include "boost/serialization/string.hpp" 
#include "boost/serialization/vector.hpp" 
#include "boost/serialization/set.hpp" 
#include "boost/serialization/unordered_map.hpp" 
//#include "boost/filesystem.hpp"

class Configuration
{
public:
  Configuration() {}
  Configuration( 
    const std::set<std::string>& mpi_functions,
    const std::set<std::string>& happens_before_orders,
    const std::set<std::string>& barrier_fns,
    const std::set<std::string>& vertex_labels,
    const std::set<std::string>& edge_labels,
    bool represent_unmatched_tests,
    bool condense_unmatched_tests,
    bool condense_matched_tests,
    bool perf_counter ) : 
      mpi_functions(mpi_functions),
      happens_before_orders(happens_before_orders),
      barrier_fns( barrier_fns ),
      vertex_labels(vertex_labels),
      edge_labels(edge_labels),
      represent_unmatched_tests(represent_unmatched_tests),
      condense_unmatched_tests(condense_unmatched_tests),
      condense_matched_tests(condense_matched_tests),
      perf_counter(perf_counter)
      { 
        for ( auto vlabel : vertex_labels ) {
          if ( vlabel == "callstack" ) {
            this->csmpi_flag = true;
          }
        }
      }

  void compute_trace_file_assignment();
  void set_trace_dirs(std::vector<std::string> trace_dirs);

  Configuration& operator=(const Configuration& rhs);
  
  // Accessors for graph-building parameters
  std::set<std::string> get_mpi_functions() const;
  std::set<std::string> get_happens_before_orders() const;
  std::set<std::string> get_barrier_fns() const;
  std::set<std::string> get_vertex_labels() const;
  std::set<std::string> get_edge_labels() const;
  bool get_represent_unmatched_tests_flag() const;
  bool get_condense_unmatched_tests_flag() const;
  bool get_condense_matched_tests_flag() const;
  bool has_csmpi() const;
  bool get_papi_flag() const; 
  // Accessors for querying trace file to dumpi_to_graph assignment 
  std::unordered_map<std::string, std::unordered_map<int,std::vector<int>>> get_dir_to_trace_rank_assignments() const;
  std::unordered_map<std::string, std::unordered_map<int,std::vector<std::string>>> get_dir_to_trace_file_assignments() const;
  std::vector<std::string> get_trace_dirs() const;
  std::unordered_map<int, std::vector<std::string>> get_trace_file_assignment(std::string trace_dir) const;
  std::unordered_map<int, std::vector<int>> get_trace_rank_assignment(std::string trace_dir) const;
  std::unordered_map<int,int> get_trace_rank_to_owning_rank() const;

  int lookup_owning_rank( int trace_rank ) const;

  friend std::ostream& operator<< (std::ostream& out, 
                                   const Configuration& config)
  {
    out << std::endl << "DUMPI-to-Graph Configuration:" << std::endl; 
    out << "Event Graph Options:" << std::endl;
    out << "\t- Events represented in graph" << std::endl;
    for (auto fn : config.mpi_functions) {
      out << "\t\t- " << fn << std::endl;
    }
    out << "\t- Orders represented in graph" << std::endl;
    for (auto order : config.happens_before_orders) {
      out << "\t\t- " << order << std::endl;
    }
    out << "\t- Vertex labels:" << std::endl;
    for (auto vertex_label : config.vertex_labels) {
      out << "\t\t- " << vertex_label << std::endl;
    }
    out << "\t- Edge labels:" << std::endl;
    for (auto edge_label : config.edge_labels) {
      out << "\t\t- " << edge_label << std::endl;
    }
    out << "Graph Construction Options:" << std::endl;
    out << std::boolalpha;
    out << "\t- Represent unmatched tests? " 
        << config.represent_unmatched_tests << std::endl;
    out << "\t- Condense unmatched tests? "   
        << config.condense_unmatched_tests << std::endl;
    out << "\t- Condense matched tests? " 
        << config.condense_matched_tests << std::endl;
    out << "Trace Directories:" << std::endl;
    for (auto trace_dir : config.trace_dirs) {
      out << "\t- " << trace_dir << std::endl;
    }
    return out;
  }

private:
  // Each dumpi_to_graph process will eventually need a global view of which
  // dumpi_to_graph processes are handling which trace files
  std::unordered_map<std::string, std::unordered_map<int,std::vector<int>>> dir_to_trace_rank_assignments;
  std::unordered_map<std::string, std::unordered_map<int,std::vector<std::string>>> dir_to_trace_file_assignments;
 
  std::unordered_map<int,int> trace_rank_to_owning_rank;

  // Define data sources for the event graph construction
  std::vector<std::string> trace_dirs;
  
  // Define what events to represent in the event graph and how to represent 
  // them
  std::set<std::string> mpi_functions;
  
  // Define what MPI functions will be modeled as a barrier
  std::set<std::string> barrier_fns;

  // Define which ordering relationships between events will be represented by
  // edges in the event graph
  std::set<std::string> happens_before_orders;

  // Define which data we will associate with each vertex in the event graph
  std::set<std::string> vertex_labels;

  // Define which data we will associate with each edge in the event graph
  std::set<std::string> edge_labels; 

  // Are we going to represent unmatched tests at all? 
  bool represent_unmatched_tests = false; 
  // If we do, will each unmatched test, of which there could be very many, 
  // be represented by a separate vertex, or will we condense consecutive 
  // unmatched tests into a single vertex? 
  bool condense_unmatched_tests = false; 
  // Will we do a similar condensing into a single vertex for matching functions
  // that can result in more than one match? (e.g., MPI_Testsome) 
  bool condense_matched_tests = false;

  // A flag to determine whether we will ingest CSMPI traces in addition to 
  // DUMPI traces. This will be set to true if "callstack" is a requested 
  // vertex label. 
  bool csmpi_flag = false;
  
  // A flag to check for use of PAPI counter in vertex labelling. 
  bool perf_counter = false; 
  // Boost serialization stuff so that the root process can read in the config
  // and broadcast it to everyone else. 
  friend class boost::serialization::access; 
  template <typename Archive>
  void serialize( Archive& ar, const unsigned int version ) 
  {
    ar & mpi_functions; 
    ar & barrier_fns;
    ar & happens_before_orders;
    ar & vertex_labels;
    ar & edge_labels;
    ar & represent_unmatched_tests;
    ar & condense_unmatched_tests;
    ar & condense_matched_tests;
    ar & trace_dirs;
    ar & dir_to_trace_rank_assignments;
    ar & dir_to_trace_file_assignments;
    ar & trace_rank_to_owning_rank;
    ar & csmpi_flag;
    ar & perf_counter;
  }

};

// Helper function to parse the JSON configuration file
Configuration parse_config_file( std::string config_file_path );

// Called in main to construct a configuration object from command-line args
Configuration parse_args( int argc, char** argv );

// Used to distribute the configuration object to all dumpi_to_graph processes
// after the root process has read the configuration file and constructed it.
void broadcast_config( Configuration& config );


#endif // D2G_CONFIGURATION_H
