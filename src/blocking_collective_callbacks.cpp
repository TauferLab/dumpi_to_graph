#include "blocking_collective_callbacks.hpp"

// DUMPI
#include "dumpi/common/argtypes.h"

// Internal
#include "Trace.hpp"
#include "Utilities.hpp"

int cb_MPI_Barrier(const dumpi_barrier *prm, 
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
  dumpi_perfinfo counters = *perf;
  dumpi_barrier event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );
  if(papi){
    trace->register_papi_struct(counters);
  }
  // Associate this barrier event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Barrier" );
  trace->associate_event_with_call( "MPI_Barrier", event_vertex_id );

  // Return OK
  return 0;
}

int cb_MPI_Reduce(const dumpi_reduce *prm, 
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
  dumpi_reduce event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );
  
  // Associate this barrier event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Reduce" );
  trace->associate_event_with_call( "MPI_Reduce", event_vertex_id );
  
  // Return OK
  return 0;
}

int cb_MPI_Allreduce(const dumpi_allreduce *prm, 
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
  dumpi_allreduce event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );
  
  // Associate this barrier event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Allreduce" );
  trace->associate_event_with_call( "MPI_Allreduce", event_vertex_id );

  // Return OK
  return 0;
}

int cb_MPI_Alltoall(const dumpi_alltoall *prm, 
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
  dumpi_alltoall event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );
  
  // Associate this barrier event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Alltoall" );
  trace->associate_event_with_call( "MPI_Alltoall", event_vertex_id );

  // Return OK
  return 0;
}

int cb_MPI_Alltoallv(const dumpi_alltoallv *prm, 
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
  dumpi_alltoallv event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );

  // Associate this barrier event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Alltoallv" );
  trace->associate_event_with_call( "MPI_Alltoallv", event_vertex_id );

  // Return OK
  return 0;
}
