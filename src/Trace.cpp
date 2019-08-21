#include "Trace.hpp"

#include <string>
#include <vector>
#include <set>
#include <unordered_map>

// DUMPI
#include "dumpi/common/argtypes.h" 

// Internal
#include "Logging.hpp"
#include "Channel.hpp"

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
  
void Trace::register_event( const dumpi_send send_event, 
                            const dumpi_time cpu_time,
                            const dumpi_time wall_time )
{
#ifdef REPORT_PROGRESS
  if ( dumpi_to_graph_rank == REPORTING_RANK ) {
    std::cout << "dumpi_to_graph_rank: " << this->get_dumpi_to_graph_rank() 
              << " handling trace_rank: " << this->get_trace_rank()
              << " in: register_event for MPI_Send" << std::endl;
  }
#endif
  // A send always has enough data to construct its corresponding channel, so we
  // do this now
  Channel channel( this->get_trace_rank(), send_event ); 
  // Register this send in the channel-specific send sequence
  auto channel_search = this->channel_to_send_seq.find( channel );
  if (channel_search != this->channel_to_send_seq.end()) {
    channel_search->second.push_back(0);
  } else {
    std::vector<int> sends;
    sends.push_back(0);
    this->channel_to_send_seq.insert( { channel, sends } );
  } 
}

void Trace::register_event( const dumpi_recv recv_event, 
                            const dumpi_time cpu_time,
                            const dumpi_time wall_time )
{
#ifdef REPORT_PROGRESS
  if ( dumpi_to_graph_rank == REPORTING_RANK ) {
    std::cout << "dumpi_to_graph_rank: " << this->get_dumpi_to_graph_rank() 
              << " handling trace_rank: " << this->get_trace_rank()
              << " in: register_event for MPI_Recv" << std::endl;
  }
#endif
}
  
void Trace::register_event( const dumpi_isend isend_event, 
                            const dumpi_time cpu_time,
                            const dumpi_time wall_time )
{
#ifdef REPORT_PROGRESS
  if ( dumpi_to_graph_rank == REPORTING_RANK ) {
    std::cout << "dumpi_to_graph_rank: " << this->get_dumpi_to_graph_rank() 
              << " handling trace_rank: " << this->get_trace_rank()
              << " in: register_event for MPI_Isend" << std::endl;
  }
#endif
  // Unpack isend event
  int request_id = isend_event.request;
  // A send always has enough data to construct its corresponding channel, so we
  // do this now
  Channel channel( this->get_trace_rank(), isend_event ); 
  // Register this send in the channel-specific send sequence
  auto channel_search = this->channel_to_send_seq.find( channel );
  if (channel_search != this->channel_to_send_seq.end()) {
    channel_search->second.push_back(0);
  } else {
    std::vector<int> sends;
    sends.push_back(0);
    this->channel_to_send_seq.insert( { channel, sends } );
  } 
}

void Trace::register_event( const dumpi_irecv irecv_event, 
                            const dumpi_time cpu_time,
                            const dumpi_time wall_time )
{
#ifdef REPORT_PROGRESS
  if ( dumpi_to_graph_rank == REPORTING_RANK ) {
    std::cout << "dumpi_to_graph_rank: " << this->get_dumpi_to_graph_rank() 
              << " handling trace_rank: " << this->get_trace_rank()
              << " in: register_event for MPI_Irecv" << std::endl;
  }
#endif
  // Unpack irecv event
  int request_id = irecv_event.request;
}
