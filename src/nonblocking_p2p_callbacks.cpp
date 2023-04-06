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

int cb_MPI_Irecv(const dumpi_irecv *prm, 
                 uint16_t thread, 
                 const dumpi_time *cpu, 
                 const dumpi_time *wall, 
                 const dumpi_perfinfo *perf, 
                 void *uarg)
{
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type;
  long req_addr, event_num;
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_perfinfo counters;
  if(papi){
    dumpi_perfinfo counters = *perf;
  }
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
  long request_id = event.request;

  trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
  if(msg_type != 1 || call_type != 1){
      std::cerr << "Misaligned Pluto output in Irecv, found " << call_type << " " << event_num << std::endl;
  }

// JACK_ Modify to add idFK instead of ID
  Request request(1, req_addr, partial_channel );
  
  trace->register_request( req_addr, request );
  
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
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type;
  long req_addr, event_num;

  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_perfinfo counters; 
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

  // Associate this send event with a timestamp
  trace->register_dumpi_timestamp( wall_time );
  if(papi){
    counters = *perf;
    trace->register_papi_struct(counters);
  } 
  // Associate this send event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Isend" );
  trace->associate_event_with_call( "MPI_Isend", event_vertex_id );

  // Since this is a non-blocking send, create and track its associated request
  long request_id = event.request;

  trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
  if(msg_type != 0 || call_type != 0){
      std::cerr << "Misaligned Pluto output in isend found " << call_type << " " << event_num << " ID: " << req_addr << " rank: " <<trace->get_trace_rank() << std::endl;
  }

  Request request( 0, req_addr, channel );
  trace->register_request( req_addr, request );
  
  // Return OK
  return 0;
}

