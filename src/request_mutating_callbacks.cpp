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
  validate_dumpi_event(prm, cpu, wall, perf);
  dumpi_cancel event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
    
  // Mark the request as cancelled
  int request_id = event.request;
  trace->cancel_request( request_id );
  
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
  validate_dumpi_event(prm, cpu, wall, perf);
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
