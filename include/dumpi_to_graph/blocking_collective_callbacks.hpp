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


int cb_MPI_Reduce(const dumpi_reduce *prm, 
                  uint16_t thread, 
                  const dumpi_time *cpu, 
                  const dumpi_time *wall, 
                  const dumpi_perfinfo *perf, 
                  void *uarg);

int cb_MPI_Allreduce(const dumpi_allreduce *prm, 
                     uint16_t thread, 
                     const dumpi_time *cpu, 
                     const dumpi_time *wall, 
                     const dumpi_perfinfo *perf, 
                     void *uarg);

int cb_MPI_Alltoall(const dumpi_alltoall *prm, 
                    uint16_t thread, 
                    const dumpi_time *cpu, 
                    const dumpi_time *wall, 
                    const dumpi_perfinfo *perf, 
                    void *uarg);

int cb_MPI_Alltoallv(const dumpi_alltoallv *prm, 
                     uint16_t thread, 
                     const dumpi_time *cpu, 
                     const dumpi_time *wall, 
                     const dumpi_perfinfo *perf, 
                     void *uarg); 

int cb_MPI_Bcast(const dumpi_bcast *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg);


#endif // D2G_BLOCKING_COLLECTIVE_CALLBACKS_H
