#include "blocking_p2p_callbacks.hpp"

#include "mpi.h" 

// DUMPI
#include "dumpi/common/argtypes.h"
#include "dumpi/common/constants.h"

// Internal 
#include "Trace.hpp"
#include "Utilities.hpp"
#include "Channel.hpp"
#include "Request.hpp"

int cb_MPI_Recv(const dumpi_recv *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg) 
{
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  validate_dumpi_event(prm, cpu, wall, perf);
  dumpi_recv event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Since MPI_Recv is blocking, we will represent this event with a recv vertex
  // in the event graph. We first need to get the channel that this recv 
  // occurred in. 
  // FIXME: If the recv is a wildcard (i.e., the traced application used 
  // MPI_ANY_SOURCE or MPI_ANY_TAG) *and* the status was ignored (i.e., the 
  // application used MPI_STATUS_IGNORE) then it is not possible to determine
  // the channel of the receive and thus the graph cannot be constructed 
  // unambiguously.
  Channel channel = trace->determine_channel_of_recv( event );

  // Create the event and add it to the event sequence for this trace
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Associate this receive event with its channel
  trace->register_recv( channel, event_vertex_id );

  // Associate this receive event with a timestamp
  trace->register_dumpi_timestamp( wall_time );

  // Associate this receive event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Recv" );
  trace->associate_event_with_call( "MPI_Recv", event_vertex_id );
  
  // Return OK
  return 0;
}

int cb_MPI_Send(const dumpi_send *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg) 
{
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  validate_dumpi_event(prm, cpu, wall, perf);
  dumpi_send event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // A send always has an unambiguous channel at the call so just construct it
  Channel channel( trace->get_trace_rank(), event );

  // Create the event and add it to the event sequence for this trace
  //size_t predecessor_vertex_id = trace->get_event_seq().back()->get_vertex_id();
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Associate this send event with its channel
  trace->register_send( channel, event_vertex_id );
  
  // Associate this send event with a timestamp
  trace->register_dumpi_timestamp( wall_time );
  
  // Associate this send event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Send" );
  trace->associate_event_with_call( "MPI_Send", event_vertex_id );
  
  // Return OK
  return 0;
}
