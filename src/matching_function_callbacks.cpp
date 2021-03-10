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
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type; 
  long req_addr, event_num;
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_wait event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Update the call index for MPI_Wait
  trace->update_call_idx( "MPI_Wait" );

  // Complete the request
  long request_id = event.request;
  
  trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
  if(call_type != 2){
     std::cerr << "Misaligned Pluto Output in Wait " << call_type << " " << event_num << std::endl;
  }  

  dumpi_status* status_ptr = event.status;
  trace->complete_request( req_addr, 
                           status_ptr, 
                           cpu_time, 
                           wall_time, 
                           perf,
                           "MPI_Wait" );

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
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type; 
  long req_addr, event_num;
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_waitany event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Determine which request was completed
  int request_idx = event.index;

  trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
  if(call_type != 2){
     std::cerr << "Misaligned Pluto Output in Waitany " << call_type << " " << event_num << std::endl;
  }    

  long request_id = *(event.requests + request_idx);
  dumpi_status* status_ptr = event.status;

  // Update the call index for MPI_Waitany
  trace->update_call_idx( "MPI_Waitany" );

  // Complete the request
  trace->complete_request( req_addr, 
                           status_ptr, 
                           cpu_time, 
                           wall_time,
                           perf,
                           "MPI_Waitany" );
  
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
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type; 
  long req_addr, event_num;

  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_waitsome event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Update the call index for MPI_Waitsome
  trace->update_call_idx( "MPI_Waitsome" );
  
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
    long request_id = *requests_ptr;
    auto search = request_indices.find( i );
    if (search != request_indices.end()) {

      trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
      if(call_type != 2){
        std::cerr << "Misaligned Pluto Output in Waitsome " << call_type << " " << event_num << std::endl;
      }
 
      trace->complete_request( req_addr, 
                               statuses_ptr, 
                               cpu_time, 
                               wall_time,
                               perf,
                               "MPI_Waitsome" ); 
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
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type; 
  long req_addr, event_num;

  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_waitall event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Update the call index for MPI_Waitsome
  trace->update_call_idx( "MPI_Waitall" );

  // All requests must have been completed
  int n_requests = event.count;
  dumpi_request* requests_ptr = event.requests;
  dumpi_status* statuses_ptr = event.statuses;
  for ( int i=0; i<n_requests; ++i ) {
    long request_id = *requests_ptr;

    trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
    if(call_type != 2){
      std::cerr << "Misaligned Pluto Output in Waitall " << call_type << " " << event_num << std::endl;
    }


    trace->complete_request( req_addr, 
                             statuses_ptr, 
                             cpu_time, 
                             wall_time,
                             perf,
                             "MPI_Waitall" );
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
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type; 
  long req_addr, event_num;
 
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_test event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Check whether a request was actually completed. If not, we can exit early.
  // FIXME: We want to allow user's the option of explicitly representing 
  // unmatched tests in the event graph, but it's not a high priority for now
  if ( event.flag == 0) {
    return 0;
  }
  
  // Update the call index for MPI_Test
  trace->update_call_idx( "MPI_Test" );
  
  // Otherwise, complete the request
  long request_id = event.request;
  trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
  if(call_type != 2){
      std::cerr << "Misaligned Pluto Output in Test " << call_type << " " << event_num << std::endl;
  }

 
  dumpi_status* status_ptr = event.status;
  trace->complete_request( req_addr, 
                           status_ptr, 
                           cpu_time, 
                           wall_time,
                           perf,
                           "MPI_Test" );

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
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type; 
  long req_addr, event_num;
 
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_testany event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Check whether a request was actually completed. If not, we can exit early.
  // FIXME: We want to allow user's the option of explicitly representing 
  // unmatched tests in the event graph, but it's not a high priority for now
  if ( event.flag == 0) {
    return 0;
  }
  
  // Update the call index for MPI_Testany
  trace->update_call_idx( "MPI_Testany" );

  // Otherwise, complete the request
  long request_id = *( event.requests + event.index );

  trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
  if(call_type != 2){
    std::cerr << "Misaligned Pluto Output in Testany " << call_type << " " << event_num << std::endl;
  }

 

  dumpi_status* status_ptr = event.status;
  trace->complete_request( req_addr, 
                           status_ptr, 
                           cpu_time, 
                           wall_time,
                           perf,
                           "MPI_Testany" );

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
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type; 
  long req_addr, event_num;
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_testsome event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;
  
  // Check whether a request was actually completed. If not, we can exit early.
  // FIXME: We want to allow user's the option of explicitly representing 
  // unmatched tests in the event graph, but it's not a high priority for now
  if ( event.outcount == 0) {
    return 0;
  }
  
  // Update the call index for MPI_Testsome
  trace->update_call_idx( "MPI_Testsome" );
  
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
    long request_id = *requests_ptr;
    auto search = request_indices.find( i );
    if (search != request_indices.end()) {

      trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
      if(call_type != 2){
        std::cerr << "Misaligned Pluto Output in Testsome " << call_type << " " << event_num << std::endl;
      }

      trace->complete_request( req_addr, 
                               statuses_ptr, 
                               cpu_time, 
                               wall_time,
                               perf,
                               "MPI_Testsome" ); 
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
  Trace* trace = (Trace*) uarg;
  // Check that event data is OK 
  bool papi = trace->get_papi_flag();
  int msg_type, call_type; 
  long req_addr, event_num;
  if(papi){
    validate_dumpi_event(prm, cpu, wall, perf);
  }
  else{
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_testall event = *prm;
  dumpi_time cpu_time = *cpu;
  dumpi_time wall_time = *wall;

  // Check whether a request was actually completed. If not, we can exit early.
  // FIXME: We want to allow user's the option of explicitly representing 
  // unmatched tests in the event graph, but it's not a high priority for now
  if ( event.flag == 0) {
    return 0;
  }
  
  // Update the call index for MPI_Testall
  trace->update_call_idx( "MPI_Testall" );
  
  // Otherwise, complete all requests
  int n_requests = event.count;
  dumpi_request* requests_ptr = event.requests;
  dumpi_status* statuses_ptr = event.statuses;
  for ( int i=0; i<n_requests; ++i ) {
    long request_id = *requests_ptr;
    
    trace->get_pluto_entry(msg_type, req_addr, call_type, event_num);
    if(call_type != 2){
      std::cerr << "Misaligned Pluto Output in Testall " << call_type << " " << event_num << std::endl;
    }

    
    trace->complete_request( req_addr, 
                             statuses_ptr, 
                             cpu_time, 
                             wall_time,
                             perf,
                             "MPI_Testall" );
    statuses_ptr++;
    requests_ptr++;
  }

  // Return OK
  return 0;
}
