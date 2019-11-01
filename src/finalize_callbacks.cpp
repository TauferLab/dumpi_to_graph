#include "finalize_callbacks.hpp"

// DUMPI-specific 
#include "dumpi/common/argtypes.h" 

// Internal 
#include "Trace.hpp"
#include "Utilities.hpp"
#include "Logging.hpp"

int cb_MPI_Finalize(const dumpi_finalize *prm, 
                    uint16_t thread, 
                    const dumpi_time *cpu, 
                    const dumpi_time *wall, 
                    const dumpi_perfinfo *perf, 
                    void *uarg) 
{
  // Check that event data is OK
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_finalize event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  trace->update_call_idx( "MPI_Finalize" );

  // Add the event to the event sequence for this trace  
  trace->register_finalize();
  trace->register_dumpi_timestamp( wall_time );

  (void)uarg;
  
  // Return OK
  return 0;
}

