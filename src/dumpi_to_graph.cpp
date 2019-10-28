// All hail the one true parallel runtime
#include <mpi.h>

// Standard
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
#include "Debug.hpp"
#include "EventGraph.hpp"
#include "CommunicatorManager.hpp"

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
  Configuration config;
  if ( dumpi_to_graph_rank == 0 ) {
    config = parse_args(argc, argv);
    config.compute_trace_file_assignment(); 
    broadcast_config( config );
  } else {
    broadcast_config( config );
  }

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
  for ( auto trace_dir : config.get_trace_dirs() ) {
#ifdef REPORT_TIMINGS
    double ingest_trace_dir_start_time = MPI_Wtime();
    if ( dumpi_to_graph_rank == REPORTING_RANK ) {
      std::cout << "Beginning to Ingest Trace Dir: " << trace_dir << std::endl;
    }
    mpi_rc = MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    auto trace_file_assignment = config.get_trace_file_assignment( trace_dir );
    std::vector<std::string> assigned_trace_files;
    auto search = trace_file_assignment.find( dumpi_to_graph_rank );
    if (search != trace_file_assignment.end() ) {
      assigned_trace_files = trace_file_assignment.at(dumpi_to_graph_rank);
    } else {
      std::cout << "Rank: " << dumpi_to_graph_rank << " not found" << std::endl;
    }
    
    // Maintain a mapping between ranks of the traced application and their 
    // representation as a Trace object
    std::unordered_map<int, Trace*> rank_to_trace;
    
    // For each trace file that this dumpi_to_graph process is responsible for,
    // construct a Trace object that 
    for ( auto trace_file : assigned_trace_files ) {
#ifdef REPORT_TIMINGS
      double ingest_trace_file_start_time = MPI_Wtime();
#endif
      // Determine rank from the traced run that this trace file represents
      int trace_file_rank = trace_file_to_rank( trace_file );

      // Construct object to represent the entire trace file
      Trace* trace_ptr = new Trace( config, 
                                    trace_dir, 
                                    trace_file_rank, 
                                    dumpi_to_graph_rank );

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
      undumpi_rc = undumpi_read_stream( profile, &callbacks, trace_ptr );
      
      // Close tracefile
      undumpi_close( profile );

      // Associate the trace file representation to its trace file
      rank_to_trace.insert( { trace_file_rank, trace_ptr } );
        
#ifdef REPORT_TIMINGS
      double ingest_trace_file_elapsed_time = MPI_Wtime() - ingest_trace_file_start_time;
      std::cout << "\tRank: " << dumpi_to_graph_rank << ", "
                << "Ingested Trace File: " << trace_file << ", "
                << "In: " << ingest_trace_file_elapsed_time << " seconds" 
                << std::endl;

#endif
    } // end loop over trace files for this trace dir

    mpi_rc = MPI_Barrier(MPI_COMM_WORLD);

#ifdef REPORT_TIMINGS
    if ( dumpi_to_graph_rank == REPORTING_RANK ) {
      double ingest_trace_dir_elapsed_time = MPI_Wtime() - ingest_trace_dir_start_time;
      std::cout << "Rank: " << dumpi_to_graph_rank << ", "
                << "Trace Dir: " << trace_dir << ", "
                << "Total Trace Dir Ingestion Time: " 
                << ingest_trace_dir_elapsed_time << " seconds"
                << std::endl;
    }
#endif
    mpi_rc = MPI_Barrier( MPI_COMM_WORLD );

#ifdef SANITY_CHECK
    // Do a big sanity check on all of the trace contents
    for ( auto kvp : rank_to_trace ) {
      assert( validate_trace( *kvp.second ) );
    }
    //std::cout << "Rank: " << dumpi_to_graph_rank 
    //          << " traces validated" << std::endl;
#endif
    mpi_rc = MPI_Barrier( MPI_COMM_WORLD );
    
    //std::cout << "Rank: " << dumpi_to_graph_rank << " exiting" << std::endl;
    //exit(0);
    
    // Construct this dumpi_to_graph process's partial view of the entire event
    // graph.
    //std::cout << "Rank: " << dumpi_to_graph_rank 
    //          << " starting event graph construction" << std::endl;
    EventGraph event_graph( config, rank_to_trace );
    std::cout << "Rank: " << dumpi_to_graph_rank 
              << " finished event graph construction" << std::endl;
    
  

    // Apply logical timestamps
    std::cout << "Rank: " << dumpi_to_graph_rank 
              << " started applying logical timestamps" << std::endl;
    event_graph.apply_scalar_logical_clock();
    std::cout << "Rank: " << dumpi_to_graph_rank 
              << " finished applying logical timestamps" << std::endl;

    mpi_rc = MPI_Barrier( MPI_COMM_WORLD );
    
    // Merge all partial views of the event graph into a single igraph graph,
    // set vertex and edge attributes, and write to disk
    event_graph.merge_and_write();

  } // End of loop over trace directories

  mpi_rc = MPI_Finalize();
}

