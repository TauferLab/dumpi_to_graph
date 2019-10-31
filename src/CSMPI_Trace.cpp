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
    std::vector<std::pair<int,CSMPI_Callstack>> callstack_seq;
    this->fn_to_callstack_seq.insert( { function_name, callstack_seq } );

    //if ( rank == REPORTING_RANK ) {
    //  std::cout << "Ingesting callstacks for function: " << function_name << std::endl;
    //}

    // Now consume callstacks until we hit another MPI function name 
    while( std::getline(input, current_line ) &&
           !std::regex_match( current_line, mpi_function_match, mpi_function_regex)) {

      //if ( rank == REPORTING_RANK ) {
      //  std::cout << "Callstack line: " << current_line << std::endl;
      //}

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
        auto pair = std::make_pair( call_idx, callstack );
        this->fn_to_callstack_seq.at( function_name ).push_back( pair );
      }
    }
  }

  //for ( auto kvp : fn_to_callstack_seq ) {
  //  std::cout << "Rank: " << rank
  //            << "\t" << kvp.first << " --> " 
  //            << kvp.second.size() << " callstacks"
  //            << std::endl;
  //}
}
  
