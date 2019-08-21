#include "nonblocking_p2p_callbacks.hpp"

#include <mpi.h>

// DUMPI
#include "dumpi/common/argtypes.h"

// Internal 
#include "Trace.hpp"
#include "Utilities.hpp"

int cb_MPI_Irecv(const dumpi_irecv *prm, 
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
    dumpi_irecv event = *prm;
    dumpi_time cpu_time = *cpu;
    dumpi_time wall_time = *wall;
    // Construct event representation in Trace
    trace->register_event(event, cpu_time, wall_time);
    return 0;
  } else {
    return -1;
  }
}

int cb_MPI_Isend(const dumpi_isend *prm, 
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
    dumpi_isend event = *prm;
    dumpi_time cpu_time = *cpu;
    dumpi_time wall_time = *wall;
    // Construct event representation in Trace
    trace->register_event(event, cpu_time, wall_time);
    return 0;
  } else {
    return -1;
  }
}

