#include "nonblocking_p2p_callbacks.hpp"

#include <mpi.h>

// DUMPI
#include "dumpi/common/argtypes.h"
#include "dumpi/common/constants.h"

// Internal 
#include "Trace.hpp"
#include "Utilities.hpp"
#include "Channel.hpp"
#include "Request.hpp"
#include "Event.hpp"

int cb_MPI_Irecv(const dumpi_irecv *prm, 
                 uint16_t thread, 
                 const dumpi_time *cpu, 
                 const dumpi_time *wall, 
                 const dumpi_perfinfo *perf, 
                 void *uarg)
{
  // Check that event data is OK 
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_irecv event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // We might not have the whole channel at the time of calling MPI_Irecv, but
  // construct one with the parts we're guaranteed to have, destination rank
  // and communicator ID
  Channel partial_channel( trace->get_trace_rank(), event.comm );

  // Might as well check if we complete the channel now
  if ( event.source != DUMPI_ANY_SOURCE ) {
    partial_channel.set_src( event.source );
  }
  if ( event.tag != DUMPI_ANY_TAG ) {
    partial_channel.set_tag( event.tag );
  } 
  
  // Since this is a non-blocking receive, we don't model it with a vertex, but
  // we do create a request so we can track its completion via a matching 
  // function later on
  int request_id = event.request;
  Request request(1, request_id, partial_channel );
  trace->register_request( request_id, request );

  // Return OK
  return 0;
}

int cb_MPI_Isend(const dumpi_isend *prm, 
                 uint16_t thread, 
                 const dumpi_time *cpu, 
                 const dumpi_time *wall, 
                 const dumpi_perfinfo *perf, 
                 void *uarg) 
{
  // Check that event data is OK 
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_isend event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // A send always has an unambiguous channel at the call so just construct it
  Channel channel( trace->get_trace_rank(), event );
  
  // Create the event and add it to the event sequence for this trace
  // FIXME: So technically... this send might not actually happen. It *might* 
  // have had MPI_Cancel called on its request... and then had a matching 
  // function called on that same request... and that might have happened fast
  // enough that the send wasn't able to transfer any of its data beforehand.
  // But cancelling isend requests is pure evil and we aren't going to support
  // modeling executions that engage in that kind of behavior. 
  //size_t predecessor_vertex_id = trace->get_event_seq().back()->get_vertex_id();
  size_t event_vertex_id = trace->get_next_vertex_id();
  
  // Associate this send event with its channel
  trace->register_send( channel, event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );

  // Since this is a non-blocking send, create and track its associated request
  int request_id = event.request;
  Request request( 0, request_id, channel );
  trace->register_request( request_id, request );
  
  // Return OK
  return 0;
}

