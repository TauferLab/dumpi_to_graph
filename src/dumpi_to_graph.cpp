// All hail the one true parallel runtime
#include <mpi.h>

#include <stdexcept>

// DUMPI
#include "dumpi/common/types.h"
#include "dumpi/libundumpi/callbacks.h" 
#include "dumpi/libundumpi/libundumpi.h"

// Internal 
#include "Logging.hpp"
#include "Configuration.hpp" 
#include "undumpi_callbacks.hpp"
#include "Trace.hpp" 
#include "Utilities.hpp"

int main(int argc, char** argv) 
{
  // Set up MPI
  int mpi_rc; 
  mpi_rc = MPI_Init(&argc, &argv);
  int dumpi_to_graph_rank, comm_size;
  mpi_rc = MPI_Comm_rank(MPI_COMM_WORLD, &dumpi_to_graph_rank);
  mpi_rc = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // Root dumpi_to_graph process parses command line arguments, constructs
  // Configuration object and broadcasts it to all other dumpi_to_graph processes
  d2g_Configuration config;
  if ( dumpi_to_graph_rank == 0 ) {
    config = parse_args(argc, argv);
    broadcast_config( config );
  } else {
    broadcast_config( config );
  }
  config.compute_trace_file_assignment(); 

  // Ensure that all dumpi_to_graph processes have the Configuration 
  mpi_rc = MPI_Barrier(MPI_COMM_WORLD);

#ifdef REPORT_PROGRESS
  if ( dumpi_to_graph_rank == REPORTING_RANK ) {
    std::cout << config << std::endl;
  }
#endif

  // Set up callbacks for parsing MPI event stream from trace files
  libundumpi_callbacks callbacks; 
  libundumpi_clear_callbacks(&callbacks);
  set_callbacks(&callbacks, config); 

  // Loop over the trace directories, constructing one event graph per
  // trace directory
  for ( auto kvp : config.get_trace_files() ) {
    auto trace_dir = kvp.first;
    auto assigned_trace_files = kvp.second;
    
    // Maintain a mapping between ranks of the traced application and their 
    // representation as a Trace object
    std::unordered_map<int, Trace*> rank_to_trace;
    
    // For each trace file that this dumpi_to_graph process is responsible for,
    // construct a Trace object that 
    for ( auto trace_file : assigned_trace_files ) {

      // Determine rank from the traced run that this trace file represents
      int trace_file_rank = trace_file_to_rank( trace_file );

#ifdef REPORT_PROGRESS
      if ( dumpi_to_graph_rank == REPORTING_RANK ) {
        std::cout << std::endl;
        std::cout << "dumpi_to_graph rank: " << dumpi_to_graph_rank 
                  << " will ingest trace file: " << trace_file << std::endl;
        std::cout << std::endl;
      }
#endif

      // Construct object to represent the entire trace file
      Trace trace( config, trace_file_rank, dumpi_to_graph_rank );

      // Check that the trace file is readable
      dumpi_profile* profile;
      profile = undumpi_open( trace_file.c_str() );
      if ( !profile ) {
        std::stringstream ss;
        ss << "Rank: " << dumpi_to_graph_rank 
           << " undumpi_open failed for trace file: " << trace_file 
           << std::endl;
        throw std::runtime_error( ss.str() );
      }   

      // Read the stream of MPI events
      int undumpi_rc;
      undumpi_rc = undumpi_read_stream( profile, &callbacks, &trace );

      // Close tracefile
      undumpi_close( profile );

#ifdef REPORT_PROGRESS
      if ( dumpi_to_graph_rank == REPORTING_RANK ) {
        std::cout << std::endl;
        std::cout << "dumpi_to_graph rank: " << dumpi_to_graph_rank 
                  << " successfully ingested trace file: " << trace_file 
                  << std::endl;
        std::cout << std::endl;
      }
#endif

      // Associate the trace file representation to its trace file
      rank_to_trace.insert( { trace_file_rank, &trace } );
    }
  }
  mpi_rc = MPI_Finalize();
}

