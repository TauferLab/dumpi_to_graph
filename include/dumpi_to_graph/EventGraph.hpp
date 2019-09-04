#ifndef D2G_EVENT_GRAPH_H
#define D2G_EVENT_GRAPH_H

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>

#include "Channel.hpp"
#include "Trace.hpp"

// Each dumpi_to_graph process will maintain one EventGraph object that contains
// all of the event sequence and channel map data from the possibly several 
// DUMPI trace files it is responsible for
class EventGraph
{
public:
  EventGraph( const Configuration& config,
              const std::unordered_map<int,Trace*> rank_to_trace );
  
  // Function for applying scalar logical timestamps
  void apply_scalar_logical_clock();

  // Convenience functions for printing state
  void report_program_order_edges() const;
  void report_message_order_edges() const;

private:
  void make_program_order_edges();
  void make_message_order_edges();
  void exchange_local_message_matching_data();
  void exchange_remote_message_matching_data();
  void make_collective_edges();

  // Data structures that define the event graph itself 
  // (or at least this dumpi_to_graph process's partial view of the event graph)
  std::vector<size_t> vertex_ids;
  std::vector<std::pair<size_t,size_t>> program_order_edges;
  std::vector<std::pair<size_t,size_t>> message_order_edges;

  // Data structures for event graph labels 
  std::unordered_map<size_t,size_t> logical_timestamps; // FIXME: scalar-only for right now

  // Data directly copied from configuration and traces
  Configuration config;
  std::unordered_map<int,Trace*> rank_to_trace;

  // Data collected by merging corresponding data structures from all traces
  // this dumpi_to_graph is handling
  channel_map channel_to_send_seq;
  channel_map channel_to_recv_seq;
  std::unordered_map<size_t,Channel> vertex_id_to_channel;
  std::unordered_map<size_t,uint8_t> vertex_id_to_event_type;
};


#endif // D2G_EVENT_GRAPH_H
