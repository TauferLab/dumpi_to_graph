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
#include "boost/filesystem.hpp"

class d2g_Configuration
{
public:
  d2g_Configuration() {}
  d2g_Configuration( 
    const std::set<std::string>& mpi_functions,
    const std::set<std::string>& happens_before_orders,
    const std::set<std::string>& vertex_labels,
    const std::set<std::string>& edge_labels,
    bool represent_unmatched_tests,
    bool condense_unmatched_tests,
    bool condense_matched_tests ) : 
      mpi_functions(mpi_functions),
      happens_before_orders(happens_before_orders),
      vertex_labels(vertex_labels),
      edge_labels(edge_labels),
      represent_unmatched_tests(represent_unmatched_tests),
      condense_unmatched_tests(condense_unmatched_tests),
      condense_matched_tests(condense_matched_tests)
      {}

  void compute_trace_file_assignment();
  void set_trace_dirs(std::vector<std::string> trace_dirs);

  d2g_Configuration& operator=(const d2g_Configuration& rhs);
  // Accessors for assignment operator implementation
  std::set<std::string> get_mpi_functions() const;
  std::set<std::string> get_happens_before_orders() const;
  std::set<std::string> get_vertex_labels() const;
  std::set<std::string> get_edge_labels() const;
  std::vector<std::string> get_trace_dirs() const;
  std::unordered_map<std::string, std::vector<std::string>> get_trace_files() const;
  bool get_represent_unmatched_tests_flag() const;
  bool get_condense_unmatched_tests_flag() const;
  bool get_condense_matched_tests_flag() const;
  
  friend std::ostream& operator<< (std::ostream& out, 
                                   const d2g_Configuration& config)
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
    out << "Assigned Trace Files:" << std::endl;
    for (auto kvp : config.trace_files) {
      std::cout << "\t- Trace files for trace dir:" << kvp.first << std::endl;
      for (auto trace_file : kvp.second) {
        out << "\t\t- " << trace_file << std::endl;
      }
    }
    return out;
  }

private:
  // Define data sources for the event graph construction
  std::vector<std::string> trace_dirs;
  std::unordered_map<std::string, std::vector<std::string> > trace_files;
  // Define what events to represent in the event graph and how to represent 
  // them
  std::set<std::string> mpi_functions;
  std::set<std::string> happens_before_orders;
  std::set<std::string> vertex_labels;
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
  
  // Boost serialization stuff so that the root process can read in the config
  // and broadcast it to everyone else. 
  friend class boost::serialization::access; 
  template <typename Archive>
  void serialize( Archive& ar, const unsigned int version ) 
  {
    ar & mpi_functions; 
    ar & happens_before_orders;
    ar & vertex_labels;
    ar & edge_labels;
    ar & represent_unmatched_tests;
    ar & condense_unmatched_tests;
    ar & condense_matched_tests;
    ar & trace_dirs;
    ar & trace_files;
  }
  
};

d2g_Configuration parse_config_file( std::string config_file_path );

d2g_Configuration parse_args( int argc, char** argv );

void broadcast_config( d2g_Configuration& config );


#endif // D2G_CONFIGURATION_H
