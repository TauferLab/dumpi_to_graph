#ifndef NONBLOCKING_P2P_CALLBACKS_H
#define NONBLOCKING_P2P_CALLBACKS_H

// DUMPI
#include "dumpi/common/argtypes.h"

int cb_MPI_Isend(const dumpi_isend *prm, 
                 uint16_t thread, 
                 const dumpi_time *cpu, 
                 const dumpi_time *wall, 
                 const dumpi_perfinfo *perf, 
                 void *uarg);

int cb_MPI_Irecv(const dumpi_irecv *prm, 
                 uint16_t thread, 
                 const dumpi_time *cpu, 
                 const dumpi_time *wall, 
                 const dumpi_perfinfo *perf, 
                 void *uarg);


#endif // NONBLOCKING_P2P_CALLBACKS_H
