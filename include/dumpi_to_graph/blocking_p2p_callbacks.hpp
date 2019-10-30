#ifndef D2G_BLOCKING_P2P_CALLBACKS_H
#define D2G_BLOCKING_P2P_CALLBACKS_H

// DUMPI-specific headers
#include "dumpi/common/argtypes.h"

int cb_MPI_Send(const dumpi_send *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg);

int cb_MPI_Recv(const dumpi_recv *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg); 

#endif // D2G_BLOCKING_P2P_CALLBACKS_H
