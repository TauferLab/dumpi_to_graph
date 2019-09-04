#ifndef TRACE_H
#define TRACE_H

#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <cstdint> // uint8_t for event types

// Boost 
#include "boost/mpi.hpp"

// DUMPI
#include "dumpi/common/argtypes.h" 

// Internal
#include "Configuration.hpp"
#include "Channel.hpp"
#include "Event.hpp"

class Trace
{
public:
  Trace( Configuration config, 
         std::string trace_dir,
         int trace_rank,
         int dumpi_to_graph_rank ) : 
    config( config ),
    trace_dir( trace_dir ),
    trace_rank( trace_rank ),
    dumpi_to_graph_rank( dumpi_to_graph_rank )
    {}

  std::string get_trace_dir() const;
  size_t get_next_vertex_id();
  size_t get_curr_vertex_id() const;
  int get_trace_rank() const;
  int get_dumpi_to_graph_rank() const;
  size_t get_initial_vertex_id() const;
  size_t get_final_vertex_id() const;
  size_t get_vertex_id_offset() const;
  std::vector<uint8_t> get_event_seq() const;
  std::unordered_map<size_t,Channel> get_vertex_id_to_channel() const;

  //void update_event_seq( Event* event_ptr );
  void update_event_seq( size_t vertex_id );

  void register_recv( const Channel& channel, size_t recv_vertex_id );
  void register_send( const Channel& channel, size_t send_vertex_id );
  void register_init();
  void register_finalize();

  void register_request( int request_id, const Request& request );

  Channel determine_channel_of_irecv( const Request& request, 
                                      const dumpi_status* status );
  
  Channel determine_channel_of_recv( const dumpi_recv recv );
 

  void cancel_request( int request_id );
  void free_request( int request_id );
  void complete_request( int request_id, 
                         const dumpi_status* status_ptr,
                         const dumpi_time cpu_time,
                         const dumpi_time wall_time );

  void complete_isend_request( Request request );
  void complete_irecv_request( Request request,
                               const dumpi_status* status_ptr,
                               const dumpi_time cpu_time,
                               const dumpi_time wall_time );

  void apply_vertex_id_offset( size_t offset );

  channel_map get_channel_to_recv_seq() const;
  channel_map get_channel_to_send_seq() const;
  std::unordered_map<int,Request> get_id_to_request() const;

  // Convenience functions for printing representations of the trace
  void report_event_seq();
  void report_id_to_request();
  void report_channel_to_send_seq();
  void report_channel_to_recv_seq();

private:

  Configuration config;

  std::string trace_dir;

  // The rank of the MPI process in the traced application that generated the
  // trace file to which this trace object corresponds
  int trace_rank;

  // The rank of the dumpi_to_graph process handling the trace file to which
  // this trace object corresponds
  int dumpi_to_graph_rank;

  // The vertex ID that is currently up for assignment to an event. This is 
  // incrementing every time an event is parsed from the trace file that we 
  // explicitly represent in the event graph with a vertex
  size_t curr_vertex_id = 0;

  // The ID of the initial vertex in this trace's program order. When the trace
  // is initially built up as the trace file is parsed, this will be 0 for all
  // trace ranks. However, once the trace files are all parsed and vertex ID 
  // offsets are determined, this will be equal to the sum of the offsets for 
  // all preceding trace ranks.
  size_t initial_vertex_id = 0;

  // The ID of the final vertex in this trace's program order. 
  size_t final_vertex_id;

  size_t vertex_id_offset;

  // For right now, we just indicate the type of each vertex with a numerical ID
  // 0 := send ( MPI_Send, MPI_Isend, and their specializations )
  // 1 := recv ( MPI_Recv, and the various matching functions )
  // 2 := init ( MPI_Init and MPI_Init_thread )
  // 3 := finalize ( MPI_Finalize )
  // 4 := barrier ( MPI_Barrier )
  std::vector<uint8_t> event_seq;

  // For send and recv vertices, there is an associated channel
  std::unordered_map<size_t, Channel> vertex_id_to_channel;

  std::unordered_map<int, Request> id_to_request;
  std::unordered_map<Channel, std::vector<size_t>, ChannelHash> channel_to_send_seq;
  std::unordered_map<Channel, std::vector<size_t>, ChannelHash> channel_to_recv_seq;

  // A function for completing requests during matching function registration
  void complete_request( int request_id );

};

#endif // TRACE_H
