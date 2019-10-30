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
  // Check that event data is OK 
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_barrier event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );

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
  // Check that event data is OK 
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_reduce event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );

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
  // Check that event data is OK 
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_allreduce event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );

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
  // Check that event data is OK 
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_alltoall event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );

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
  // Check that event data is OK 
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_alltoallv event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );

  // Return OK
  return 0;
}
