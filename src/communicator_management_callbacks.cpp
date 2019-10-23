#include "communicator_management_callbacks.hpp"

// DUMPI-specific headers
#include "dumpi/common/argtypes.h"

int cb_MPI_Comm_split(const dumpi_comm_split *prm, 
                      uint16_t thread, 
                      const dumpi_time *cpu, 
                      const dumpi_time *wall, 
                      const dumpi_perfinfo *perf, 
                      void *uarg) 
{
  return 0;
}

int cb_MPI_Comm_rank(const dumpi_comm_rank *prm, 
                     uint16_t thread, 
                     const dumpi_time *cpu, 
                     const dumpi_time *wall, 
                     const dumpi_perfinfo *perf, 
                     void *uarg) 
{
  return 0;
}

int cb_MPI_Comm_size(const dumpi_comm_size *prm, 
                     uint16_t thread, 
                     const dumpi_time *cpu, 
                     const dumpi_time *wall, 
                     const dumpi_perfinfo *perf, 
                     void *uarg) 
{
  return 0;
}
