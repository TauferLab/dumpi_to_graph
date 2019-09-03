#include "matching_function_callbacks.hpp" 

#include "mpi.h" 

// Standard
#include <sstream>
#include <stdexcept> 
#include <set> 

// DUMPI
#include "dumpi/common/argtypes.h"
#include "dumpi/common/constants.h"

// Internal 
#include "Trace.hpp"
#include "Channel.hpp" 
#include "Request.hpp" 
#include "Utilities.hpp"


int cb_MPI_Wait(const dumpi_wait *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg) 
{
  // Check that event data is OK
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_wait event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Complete the request
  int request_id = event.request;
  dumpi_status* status_ptr = event.status;
  trace->complete_request( request_id, status_ptr, cpu_time, wall_time );

  // Return OK
  return 0;
} 

int cb_MPI_Waitany(const dumpi_waitany *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg) 
{
  // Check that event data is OK
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_waitany event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Determine which request was completed
  int request_idx = event.index;
  int request_id = *(event.requests + request_idx);
  dumpi_status* status_ptr = event.status;

  // Complete the request
  trace->complete_request( request_id, status_ptr, cpu_time, wall_time );
  
  // Return OK
  return 0; 
}

int cb_MPI_Waitsome(const dumpi_waitsome *prm, 
                    uint16_t thread, 
                    const dumpi_time *cpu, 
                    const dumpi_time *wall, 
                    const dumpi_perfinfo *perf, 
                    void *uarg) 
{
  // Check that event data is OK
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_waitsome event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Determine which requests were completed
  int n_requests_completed = event.outcount;
  int* indices_ptr = event.indices;
  std::set<int> request_indices;
  for ( int i=0; i<n_requests_completed; ++i ) {
    request_indices.insert( *indices_ptr );
    indices_ptr++;
  }
  
  // Iterate over the requests, handling them if they were completed
  int n_requests_tested = event.count;
  dumpi_request* requests_ptr = event.requests;
  dumpi_status* statuses_ptr = event.statuses;
  for ( int i=0; i<n_requests_tested; ++i ) {
    int request_id = *requests_ptr;
    auto search = request_indices.find( i );
    if (search != request_indices.end()) {
      trace->complete_request( request_id, statuses_ptr, cpu_time, wall_time); 
      statuses_ptr++;
    }
    requests_ptr++;
  }

  // Return OK
  return 0;
}

int cb_MPI_Waitall(const dumpi_waitall *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg) 
{
  // Check that event data is OK
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_waitall event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // All requests must have been completed
  int n_requests = event.count;
  dumpi_request* requests_ptr = event.requests;
  dumpi_status* statuses_ptr = event.statuses;
  for ( int i=0; i<n_requests; ++i ) {
    int request_id = *requests_ptr;
    trace->complete_request( request_id, statuses_ptr, cpu_time, wall_time );
    statuses_ptr++;
    requests_ptr++;
  }

  // Return OK
  return 0;
}

int cb_MPI_Test(const dumpi_test *prm, 
                uint16_t thread, 
                const dumpi_time *cpu, 
                const dumpi_time *wall, 
                const dumpi_perfinfo *perf, 
                void *uarg) 
{
  // Check that event data is OK
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_test event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Check whether a request was actually completed. If not, we can exit early.
  // FIXME: We want to allow user's the option of explicitly representing 
  // unmatched tests in the event graph, but it's not a high priority for now
  if ( event.flag == 0) {
    return 0;
  }
  
  // Otherwise, complete the request
  int request_id = event.request;
  dumpi_status* status_ptr = event.status;
  trace->complete_request( request_id, status_ptr, cpu_time, wall_time );

  // Return OK
  return 0;
}

int cb_MPI_Testany(const dumpi_testany *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg) 
{
  // Check that event data is OK
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_testany event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Check whether a request was actually completed. If not, we can exit early.
  // FIXME: We want to allow user's the option of explicitly representing 
  // unmatched tests in the event graph, but it's not a high priority for now
  if ( event.flag == 0) {
    return 0;
  }

  // Otherwise, complete the request
  int request_id = *( event.requests + event.index );
  dumpi_status* status_ptr = event.status;
  trace->complete_request( request_id, status_ptr, cpu_time, wall_time );

  // Return OK
  return 0;
}

int cb_MPI_Testsome(const dumpi_testsome *prm, 
                    uint16_t thread, 
                    const dumpi_time *cpu, 
                    const dumpi_time *wall, 
                    const dumpi_perfinfo *perf, 
                    void *uarg) 
{
  // Check that event data is OK
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_testsome event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Check whether a request was actually completed. If not, we can exit early.
  // FIXME: We want to allow user's the option of explicitly representing 
  // unmatched tests in the event graph, but it's not a high priority for now
  if ( event.outcount == 0) {
    return 0;
  }
  
  // Otherwise, determine indices of requests that actually completed
  int n_requests_completed = event.outcount;
  int* indices_ptr = event.indices;
  std::set<int> request_indices;
  for ( int i=0; i<n_requests_completed; ++i ) {
    request_indices.insert( *indices_ptr );
    indices_ptr++;
  }
  
  // Iterate over the requests, handling them if they were completed
  int n_requests_tested = event.count;
  dumpi_request* requests_ptr = event.requests;
  dumpi_status* statuses_ptr = event.statuses;
  for ( int i=0; i<n_requests_tested; ++i ) {
    int request_id = *requests_ptr;
    auto search = request_indices.find( i );
    if (search != request_indices.end()) {
      trace->complete_request( request_id, statuses_ptr, cpu_time, wall_time); 
      statuses_ptr++;
    }
    requests_ptr++;
  }

  // Return OK
  return 0;
}

int cb_MPI_Testall(const dumpi_testall *prm, 
                   uint16_t thread, 
                   const dumpi_time *cpu, 
                   const dumpi_time *wall, 
                   const dumpi_perfinfo *perf, 
                   void *uarg) 
{
  // Check that event data is OK
  validate_dumpi_event(prm, cpu, wall, perf);
  Trace* trace = (Trace*) uarg;
  dumpi_testall event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Check whether a request was actually completed. If not, we can exit early.
  // FIXME: We want to allow user's the option of explicitly representing 
  // unmatched tests in the event graph, but it's not a high priority for now
  if ( event.flag == 0) {
    return 0;
  }
  
  // Otherwise, complete all requests
  int n_requests = event.count;
  dumpi_request* requests_ptr = event.requests;
  dumpi_status* statuses_ptr = event.statuses;
  for ( int i=0; i<n_requests; ++i ) {
    int request_id = *requests_ptr;
    trace->complete_request( request_id, statuses_ptr, cpu_time, wall_time );
    statuses_ptr++;
    requests_ptr++;
  }

  // Return OK
  return 0;
}
