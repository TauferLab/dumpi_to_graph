#include "blocking_p2p_callbacks.hpp"

#include "mpi.h" 

// DUMPI
#include "dumpi/common/argtypes.h"
#include "dumpi/common/constants.h"

// Internal 
#include "Trace.hpp"
#include "Utilities.hpp"

int cb_MPI_Recv(const dumpi_recv *prm, 
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
    dumpi_recv event = *prm;
    dumpi_time cpu_time = *cpu;
    dumpi_time wall_time = *wall;
    // Construct event representation in Trace
    trace->register_event(event, cpu_time, wall_time);
    return 0;
  } else {
    return -1;
  }
}

int cb_MPI_Send(const dumpi_send *prm, 
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
    dumpi_send event = *prm;
    dumpi_time cpu_time = *cpu;
    dumpi_time wall_time = *wall;
    // Construct event representation in Trace
    trace->register_event(event, cpu_time, wall_time);
    return 0;
  } else {
    return -1;
  }
}
