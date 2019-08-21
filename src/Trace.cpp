#include "Trace.hpp"

#include <string>
#include <vector>
#include <set>
#include <unordered_map>

// DUMPI-specific 
#include "dumpi/common/argtypes.h" 

// Internal
#include "Logging.hpp"

int Trace::get_trace_rank()
{
  return this->trace_rank;
}

int Trace::get_dumpi_to_graph_rank()
{
  return this->dumpi_to_graph_rank;
}

void Trace::register_event( const dumpi_init init_event, 
                            const dumpi_time cpu_time,
                            const dumpi_time wall_time )
{
#ifdef REPORT_PROGRESS
  int mpi_rc, dumpi_to_graph_rank;
  mpi_rc = MPI_Comm_rank(MPI_COMM_WORLD, &dumpi_to_graph_rank);
  if ( dumpi_to_graph_rank == REPORTING_RANK ) {
    std::cout << "Rank: " << dumpi_to_graph_rank 
              << " in register_event for MPI_Init" << std::endl;
  }
#endif
}

void Trace::register_event( const dumpi_init_thread init_thread_event, 
                            const dumpi_time cpu_time,
                            const dumpi_time wall_time )
{
#ifdef REPORT_PROGRESS
  int mpi_rc, dumpi_to_graph_rank;
  mpi_rc = MPI_Comm_rank(MPI_COMM_WORLD, &dumpi_to_graph_rank);
  if ( dumpi_to_graph_rank == REPORTING_RANK ) {
    std::cout << "Rank: " << dumpi_to_graph_rank 
              << " in register_event for MPI_Init_thread" << std::endl;
  }
#endif
}

void Trace::register_event( const dumpi_finalize finalize_event, 
                            const dumpi_time cpu_time,
                            const dumpi_time wall_time )
{
#ifdef REPORT_PROGRESS
  int mpi_rc, dumpi_to_graph_rank;
  mpi_rc = MPI_Comm_rank(MPI_COMM_WORLD, &dumpi_to_graph_rank);
  if ( dumpi_to_graph_rank == REPORTING_RANK ) {
    std::cout << "Rank: " << dumpi_to_graph_rank 
              << " in register_event for MPI_Finalize" << std::endl;
  }
#endif
}
