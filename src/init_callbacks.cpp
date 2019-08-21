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
  if ( validate_dumpi_event(prm, cpu, wall, perf) ) {
    // Get pointer to trace object in which this event will be represented
    Trace* trace = (Trace*) uarg;
    // Get DUMPI event data
    dumpi_init event = *prm;
    dumpi_time cpu_time = *cpu;
    dumpi_time wall_time = *wall;
    // Construct event representation in Trace
    trace->register_event(event, cpu_time, wall_time);
    return 0;
  } else {
    return -1;
  }
}

int cb_MPI_Init_thread(const dumpi_init_thread *prm, 
                       uint16_t thread, 
                       const dumpi_time *cpu, 
                       const dumpi_time *wall, 
                       const dumpi_perfinfo *perf, 
                       void *uarg) 
{
  if ( validate_dumpi_event(prm, cpu, wall, perf) ) {
    // Get pointer to trace object in which this event will be represented
    Trace* trace = (Trace*) uarg;
    // Get DUMPI event data
    dumpi_init_thread event = *prm;
    dumpi_time cpu_time = *cpu;
    dumpi_time wall_time = *wall;
    // Construct event representation in Trace
    trace->register_event(event, cpu_time, wall_time);
    return 0;
  } else {
    return -1;
  }
}
