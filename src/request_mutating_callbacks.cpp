#include "request_mutating_callbacks.hpp"

#include "mpi.h" 

// DUMPI
#include "dumpi/common/argtypes.h"

// Internal headers
#include "Trace.hpp"
#include "Utilities.hpp"

int cb_MPI_Cancel(const dumpi_cancel *prm,
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
  dumpi_cancel event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
  if(call_type != 7){
    std::cerr << "Misaligned Pluto output in Cancel, found " << call_type << " " << event_num << " ID: " << req_addr << " rank: " << trace->get_trace_rank() << std::endl;
  }

  // Mark the request as cancelled using the Pluto-recorded address so the
  // lookup hits the same key under which Isend/Irecv registered the Request
  trace->cancel_request( req_addr );

  // Return OK
  return 0;
}

int cb_MPI_Request_free(const dumpi_request_free *prm, 
                        uint16_t thread, 
                        const dumpi_time *cpu, 
                        const dumpi_time *wall, 
                        const dumpi_perfinfo *perf, 
                        void *uarg) 
{
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_request_free event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
    
  // Remove the request from the id_to_request map
  int request_id = event.request;
  trace->free_request( request_id );
  
  // Return OK
  return 0; 
}

int cb_MPI_Start(const dumpi_start *prm, 
                 uint16_t thread, 
                 const dumpi_time *cpu, 
                 const dumpi_time *wall, 
                 const dumpi_perfinfo *perf, 
                 void *uarg) 
{
  return 0;
}

int cb_MPI_Startall(const dumpi_startall *prm, 
                    uint16_t thread, 
                    const dumpi_time *cpu, 
                    const dumpi_time *wall, 
                    const dumpi_perfinfo *perf, 
                    void *uarg) 
{
  return 0;
}
