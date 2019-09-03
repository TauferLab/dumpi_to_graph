#ifndef D2G_REQUEST_MUTATING_CALLBACKS_H
#define D2G_REQUEST_MUTATING_CALLBACKS_H

#include "dumpi/common/argtypes.h" 

int cb_MPI_Cancel(const dumpi_cancel *prm,
                  uint16_t thread, 
                  const dumpi_time *cpu, 
                  const dumpi_time *wall, 
                  const dumpi_perfinfo *perf, 
                  void *uarg);

int cb_MPI_Request_free(const dumpi_request_free *prm, 
                        uint16_t thread, 
                        const dumpi_time *cpu, 
                        const dumpi_time *wall, 
                        const dumpi_perfinfo *perf, 
                        void *uarg);  

int cb_MPI_Start(const dumpi_start *prm, 
                 uint16_t thread, 
                 const dumpi_time *cpu, 
                 const dumpi_time *wall, 
                 const dumpi_perfinfo *perf, 
                 void *uarg);

int cb_MPI_Startall(const dumpi_startall *prm, 
                    uint16_t thread, 
                    const dumpi_time *cpu, 
                    const dumpi_time *wall, 
                    const dumpi_perfinfo *perf, 
                    void *uarg);

#endif // D2G_REQUEST_MUTATING_CALLBACKS_H
