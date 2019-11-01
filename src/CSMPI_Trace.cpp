#include "CSMPI_Trace.hpp"

#include "mpi.h"

#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <regex>

// Boost
#include "boost/tokenizer.hpp"

// Internal
#include "Logging.hpp"
#include "CSMPI_Callstack.hpp"

CSMPI_Trace::CSMPI_Trace( std::string csmpi_trace_file,
                          int trace_file_rank )
{ 
  this->trace_rank = trace_file_rank;
  this->trace_file = csmpi_trace_file;

  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#ifdef REPORT_PROGRESS
  std::cout << "\tRank: " << rank
            << " ingesting CSMPI trace file: " << csmpi_trace_file
            << std::endl;
#endif
 
  std::regex mpi_function_regex("Callstacks for MPI Function: (MPI_\\w+)");

  std::regex callstack_regex("(\\d+), ([\\w+ ]+)");

  std::ifstream input( csmpi_trace_file );
  std::string current_line;
  std::smatch mpi_function_match;
 
  std::smatch callstack_match;


  std::getline( input, current_line );
  while (std::regex_match( current_line, mpi_function_match, mpi_function_regex ) ) {
    std::ssub_match func_submatch = mpi_function_match[1];
    std::string function_name = func_submatch.str();
    // Make entry in map from function to sequence of callstacks
    std::unordered_map<int,CSMPI_Callstack> idx_to_callstack;
    this->fn_to_idx_to_callstack.insert( { function_name, idx_to_callstack } );

    // Now consume callstacks until we hit another MPI function name 
    while( std::getline(input, current_line ) &&
           !std::regex_match( current_line, mpi_function_match, mpi_function_regex)) {

      if ( std::regex_match( current_line, callstack_match, callstack_regex ) ) {

        
        // Extract the call index
        std::ssub_match call_idx_submatch = callstack_match[1];
        std::string call_idx_str = call_idx_submatch.str();
        int call_idx = std::atoi( call_idx_str.c_str() );
        
        // Extract the addresses
        std::ssub_match addresses_submatch = callstack_match[2];
        std::string addresses_str = addresses_submatch.str();
        boost::char_separator<char> delim(" ");
        boost::tokenizer<boost::char_separator<char>> tokens( addresses_str, delim );
        std::vector<std::string> addresses;
        for ( auto token : tokens ) {
          addresses.push_back( token );
        }

        // Make callstack
        CSMPI_Callstack callstack( call_idx, addresses );
  
        // Associate with function
        //auto pair = std::make_pair( call_idx, callstack );
        //this->fn_to_callstack_seq.at( function_name ).push_back( pair );
        this->fn_to_idx_to_callstack.at( function_name ).insert( { call_idx, callstack } );
      }
    }
  }
}
  
std::string CSMPI_Trace::lookup_callstack( std::string fn, int call_idx ) const
{
  // A default value to return in the even that the specified (fn, call_idx)
  // does not have a callstack associated with it in the CSMPI trace
  std::string callstack_str = "";
  auto key = std::make_pair( fn, call_idx );
  auto fn_search = this->fn_to_idx_to_callstack.find( fn );
  if ( fn_search != this->fn_to_idx_to_callstack.end() ) {
    auto idx_to_callstack = this->fn_to_idx_to_callstack.at( fn );
    auto idx_search = idx_to_callstack.find( call_idx );
    if ( idx_search != idx_to_callstack.end() ) {
      auto callstack = idx_to_callstack.at( call_idx );
      callstack_str = callstack.str();
    }
  }
  return callstack_str;
}
