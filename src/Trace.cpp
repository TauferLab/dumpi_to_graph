#include "Trace.hpp"

#include <string>
#include <vector>
#include <set>
#include <unordered_map>

// DUMPI
#include "dumpi/common/argtypes.h" 
#include "dumpi/common/constants.h" // DUMPI_ANY_SOURCE, DUMPI_ANY_TAG

// Internal
#include "Logging.hpp"
#include "Channel.hpp"
#include "Request.hpp"
#include "Event.hpp"

size_t Trace::get_next_vertex_id() 
{
  size_t vertex_id = this->curr_vertex_id;
  this->curr_vertex_id++;
  return vertex_id;
} 

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
  
  InitEvent event( this->get_next_vertex_id(), cpu_time, wall_time );
  this->event_seq.push_back( event );

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

  // Try to construct channel
  int dst_rank = this->get_trace_rank();

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
  
  // Construct an isend request
  IsendRequest request( request_id, channel );
  auto request_search = this->id_to_request.find( request_id );
  // Case 1: Request is not currently tracked so we add it to the map
  if ( request_search != this->id_to_request.end() ) {
    this->id_to_request.insert( { request_id, request } );
  }
  // Case 2: Request is currently being tracked. In this case, the trace is 
  // malformed and we should abort
  else {
    auto prev_request = request_search->second; 
    std::stringstream ss;
    ss << "dumpi_to_graph rank: " << this->get_dumpi_to_graph_rank()
       << " trying to map request ID: " << request_id 
       << " to request: " << request
       << " but request ID is already mapped to request: " << prev_request 
       << std::endl;
    throw std::runtime_error( ss.str() );
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

void Trace::register_event( const dumpi_waitall waitall_event,
                            const dumpi_time cpu_time,
                            const dumpi_time wall_time ) 
{
#ifdef REPORT_PROGRESS
  if ( dumpi_to_graph_rank == REPORTING_RANK ) {
    std::cout << "dumpi_to_graph_rank: " << this->get_dumpi_to_graph_rank() 
              << " handling trace_rank: " << this->get_trace_rank()
              << " in: register_event for MPI_Waitall" << std::endl;
  }
#endif
  
}
