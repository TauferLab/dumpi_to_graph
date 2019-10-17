#ifndef D2G_BLOCKING_COLLECTIVE_CALLBACKS_H
#define D2G_BLOCKING_COLLECTIVE_CALLBACKS_H

// DUMPI
#include "dumpi/common/argtypes.h"

int cb_MPI_Barrier(const dumpi_barrier *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg);

#endif // D2G_BLOCKING_COLLECTIVE_CALLBACKS_H
