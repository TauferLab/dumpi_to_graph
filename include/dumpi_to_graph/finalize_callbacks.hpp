#ifndef FINALIZE_CALLBACKS_H
#define FINALIZE_CALLBACKS_H

#include "dumpi/common/argtypes.h" 

int cb_MPI_Finalize(const dumpi_finalize *prm, 
                    uint16_t thread, 
                    const dumpi_time *cpu, 
                    const dumpi_time *wall, 
                    const dumpi_perfinfo *perf, 
                    void *uarg);

#endif // FINALIZE_CALLBACKS_H
