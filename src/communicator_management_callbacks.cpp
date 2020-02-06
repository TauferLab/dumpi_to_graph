#include "communicator_management_callbacks.hpp"

#include "mpi.h"

#include <iostream>

// DUMPI
#include "dumpi/common/argtypes.h"

// Internal
#include "Trace.hpp"
#include "Utilities.hpp"

int cb_MPI_Comm_split(const dumpi_comm_split *prm, 
                      uint16_t thread, 
                      const dumpi_time *cpu, 
                      const dumpi_time *wall, 
                      const dumpi_perfinfo *perf, 
                      void *uarg) 
{
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag()
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall)
  }
  dumpi_comm_split event = *prm;
  auto parent_comm_id = event.oldcomm;
  auto new_comm_id = event.newcomm;
  auto color = event.color;
  auto key = event.key; 
  trace->register_comm_split( parent_comm_id, new_comm_id, color, key );
  
  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  //std::cout << "Rank: " << rank << " in MPI_Comm_split -- "
  //          << "parent comm: " << parent_comm_id << " "
  //          << "new comm: " << new_comm_id << " "
  //          << "color: " << color << " "
  //          << "key: " << key 
  //          << std::endl;

  return 0;
}

int cb_MPI_Comm_rank(const dumpi_comm_rank *prm, 
                     uint16_t thread, 
                     const dumpi_time *cpu, 
                     const dumpi_time *wall, 
                     const dumpi_perfinfo *perf, 
                     void *uarg) 
{
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag()
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall)
  }
  dumpi_comm_rank event = *prm;

  auto comm_id = event.comm;
  auto rank = event.rank;
  //trace->register_communicator_rank( comm_id, rank );
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
