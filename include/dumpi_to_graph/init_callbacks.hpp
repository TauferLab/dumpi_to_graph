#ifndef INIT_CALLBACKS_H
#define INIT_CALLBACKS_H

#include "dumpi/common/argtypes.h" 

int cb_MPI_Init(const dumpi_init *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg);

int cb_MPI_Init_thread(const dumpi_init_thread *prm, 
                       uint16_t thread, 
                       const dumpi_time *cpu, 
                       const dumpi_time *wall, 
                       const dumpi_perfinfo *perf, 
                       void *uarg);


#endif // INIT_CALLBACKS_H
