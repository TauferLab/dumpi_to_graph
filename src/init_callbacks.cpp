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
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_init event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Add the event to the event sequence for this trace  
  trace->register_init();
  
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
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_init_thread event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Add the event to the event sequence for this trace  
  trace->register_init();
  
  // Return OK
  return 0;
}
