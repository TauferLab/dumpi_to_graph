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
  std::cerr << "DumpiToGraph::" << __func__ << " detected " << " Rank: " << trace->get_trace_rank() << std::endl;
  int msg_type, call_type;
  long req_addr, event_num;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_perfinfo counters;
  dumpi_barrier event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );
  if(papi){
    counters = *perf;
    trace->register_papi_struct(counters);
  }
  // Associate this barrier event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Barrier" );
  trace->associate_event_with_call( "MPI_Barrier", event_vertex_id );


  trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
  if(msg_type != 2 || call_type != 4){
      std::cerr << "Misaligned Pluto output in Barrier, found " << call_type << " " << event_num << " Rank: " << trace->get_trace_rank() << std::endl;
  }

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
  bool papi = trace->get_papi_flag();
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_perfinfo counters;
  dumpi_reduce event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );
  if(papi){
    counters = *perf;
    trace->register_papi_struct(counters);
  }
  
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
  std::cerr << "DumpiToGraph::" << __func__ << " detected " << " Rank: " << trace->get_trace_rank() << std::endl;

  // Check that event data is OK 
  int msg_type, call_type;
  long req_addr, event_num;
  bool papi = trace->get_papi_flag();
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_perfinfo counters;
  dumpi_allreduce event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );
  if(papi){
    counters = *perf;
    trace->register_papi_struct(counters);
  }
  
  // Associate this barrier event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Allreduce" );
  trace->associate_event_with_call( "MPI_Allreduce", event_vertex_id );
  trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
  if(msg_type != 2 || call_type != 5){
      std::cerr << "Misaligned Pluto output in AllReduce, found " << call_type << " " << event_num << " Rank: " << trace->get_trace_rank() << std::endl;
  }
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
  bool papi = trace->get_papi_flag();
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_perfinfo counters;
  dumpi_alltoall event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );
  if(papi){
    counters = *perf;
    trace->register_papi_struct(counters);
  }
  
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
  
  bool papi = trace->get_papi_flag();
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_perfinfo counters;
  dumpi_alltoallv event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_barrier( event_vertex_id );
  trace->register_dumpi_timestamp( wall_time );
  if(papi){
    counters = *perf;
    trace->register_papi_struct(counters);
  }

  // Associate this barrier event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Alltoallv" );
  trace->associate_event_with_call( "MPI_Alltoallv", event_vertex_id );

  // Return OK
  return 0;
}

int cb_MPI_Bcast(const dumpi_bcast *prm, 
                     uint16_t thread, 
                     const dumpi_time *cpu, 
                     const dumpi_time *wall, 
                     const dumpi_perfinfo *perf, 
                     void *uarg) 
{
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  
  bool papi = trace->get_papi_flag();
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }

  dumpi_perfinfo counters;
  dumpi_bcast event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // A send always has an unambiguous channel at the call so just construct it
  //Channel channel( trace->get_trace_rank(), event );

  // Get the vertex ID for the vertex that will represent this event
  size_t event_vertex_id = trace->get_next_vertex_id();

  // Add the vertex to the event sequence
  trace->register_bcast( event_vertex_id );

  // Associate this send event with a timestamp
  trace->register_dumpi_timestamp( wall_time );
  if(papi){
    counters = *perf;
    trace->register_papi_struct(counters);
  }

  // Associate this barrier event with the MPI function call that generated it
  trace->update_call_idx( "MPI_Bcast" );
  trace->associate_event_with_call( "MPI_Bcast", event_vertex_id );

  // Return OK
  return 0;
}
