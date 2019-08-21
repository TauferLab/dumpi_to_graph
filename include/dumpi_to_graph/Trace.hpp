#ifndef TRACE_H
#define TRACE_H

#include <string>
#include <vector>
#include <set>
#include <unordered_map>

// Boost 
#include "boost/mpi.hpp"

// DUMPI-specific 
#include "dumpi/common/argtypes.h" 

// Internal
#include "Configuration.hpp"
#include "Channel.hpp"

class Trace
{
public:
  Trace( d2g_Configuration config, 
         int trace_rank,
         int dumpi_to_graph_rank ) : 
    config( config ),
    trace_rank( trace_rank ),
    dumpi_to_graph_rank( dumpi_to_graph_rank )
    {}
 
  int get_trace_rank();
  int get_dumpi_to_graph_rank();

  // Event registration overloads for MPI startup and shutdown:
  // - MPI_Init
  // - MPI_Init_thread
  // - MPI_Finalize
  void register_event( const dumpi_init init_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );
  void register_event( const dumpi_init_thread init_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );
  void register_event( const dumpi_finalize finalize_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );

  // Event registration overloads for blocking point-to-point communication
  void register_event( const dumpi_send send_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );
  void register_event( const dumpi_recv recv_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );
  
  // Event registration overloads for non-blocking point-to-point communication
  void register_event( const dumpi_isend isend_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );
  void register_event( const dumpi_irecv irecv_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );


private:
  d2g_Configuration config;
  int trace_rank;
  int dumpi_to_graph_rank;
  int curr_vertex_id = 0;
  
  std::unordered_map<Channel, std::vector<int>, ChannelHash> channel_to_send_seq;

};

#endif // TRACE_H
