#include "init_callbacks.hpp"

// DUMPI-specific 
#include "dumpi/common/argtypes.h" 

// Internal 
#include "Trace.hpp"
#include "Utilities.hpp"

int cb_MPI_Init(const dumpi_init *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg) 
{
  // Check that event data is OK
  Trace* trace = (Trace*) uarg;
  bool papi = trace->get_papi_flag();
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
  dumpi_init event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  trace->update_call_idx( "MPI_Init" );
  if(papi){
    dumpi_perfinfo counters = *perf;
    trace->register_papi_struct(counters);
  }
  // Add the event to the event sequence for this trace  
  trace->register_init("MPI_Init");
  trace->register_initial_dumpi_timestamp( wall_time );

  
  // Return OK
  return 0;
}

int cb_MPI_Init_thread(const dumpi_init_thread *prm, 
                       uint16_t thread, 
                       const dumpi_time *cpu, 
                       const dumpi_time *wall, 
                       const dumpi_perfinfo *perf, 
                       void *uarg) 
{
  // Check that event data is OK
  Trace* trace = (Trace*) uarg;
  bool papi = trace->get_papi_flag();
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
  dumpi_init_thread event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  trace->update_call_idx( "MPI_Init_thread" );
  if(papi){
    trace->register_papi_struct(counters);
  }
  // Add the event to the event sequence for this trace  
  trace->register_init("MPI_Init_thread");
  trace->register_initial_dumpi_timestamp( wall_time );
  
  // Return OK
  return 0;
}
