#ifndef MATCHING_FUNCTION_CALLBACKS_H
#define MATCHING_FUNCTION_CALLBACKS_H

// DUMPI
#include "dumpi/common/argtypes.h"
#include "dumpi/common/constants.h"

//// Internal headers
//#include "trace_monitor.hpp"
//#include "request_properties.hpp" 

//// Helper functions
//void complete_request( trace_monitor* tm, 
//                       int request_id, 
//                       dumpi_status* status_ptr,
//                       dumpi_time wtime
//                     );

// Callbacks for blocking matching functions
int cb_MPI_Wait(const dumpi_wait *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg);
int cb_MPI_Waitany(const dumpi_waitany *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg);
int cb_MPI_Waitsome(const dumpi_waitsome *prm, 
                    uint16_t thread, 
                    const dumpi_time *cpu, 
                    const dumpi_time *wall, 
                    const dumpi_perfinfo *perf, 
                    void *uarg); 
int cb_MPI_Waitall(const dumpi_waitall *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg);

// Callbacks for non-blocking matching functions
int cb_MPI_Test(const dumpi_test *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg);
int cb_MPI_Testany(const dumpi_testany *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg);
int cb_MPI_Testsome(const dumpi_testsome *prm, 
                    uint16_t thread, 
                    const dumpi_time *cpu, 
                    const dumpi_time *wall, 
                    const dumpi_perfinfo *perf, 
                    void *uarg);
int cb_MPI_Testall(const dumpi_testall *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg); 

#endif // MATCHING_FUNCTION_CALLBACKS_H
