#ifndef D2G_COMMUNICATOR_MANAGEMENT_CALLBACKS_H
#define D2G_COMMUNICATOR_MANAGEMENT_CALLBACKS_H

// DUMPI-specific headers
#include "dumpi/common/argtypes.h"

int cb_MPI_Comm_split(const dumpi_comm_split *prm, 
                      uint16_t thread, 
                      const dumpi_time *cpu, 
                      const dumpi_time *wall, 
                      const dumpi_perfinfo *perf, 
                      void *uarg);

int cb_MPI_Comm_rank(const dumpi_comm_rank *prm, 
                     uint16_t thread, 
                     const dumpi_time *cpu, 
                     const dumpi_time *wall, 
                     const dumpi_perfinfo *perf, 
                     void *uarg);

int cb_MPI_Comm_size(const dumpi_comm_size *prm, 
                     uint16_t thread, 
                     const dumpi_time *cpu, 
                     const dumpi_time *wall, 
                     const dumpi_perfinfo *perf, 
                     void *uarg);

#endif // D2G_COMMUNICATOR_MANAGEMENT_CALLBACKS_H
