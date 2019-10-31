#ifndef D2G_EVENT_GRAPH_H
#define D2G_EVENT_GRAPH_H

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>
#include <utility>
#include <functional>

// Boost
#include "boost/functional/hash.hpp"
#include "boost/serialization/access.hpp"

// Internal
#include "Channel.hpp"
#include "Trace.hpp"
#include "CommunicatorManager.hpp" 
#include "CSMPI_Trace.hpp"

struct pair_hash
{
  template<class T1, class T2>
  std::size_t operator() (const std::pair<T1,T2> &pair) const 
  {
    std::size_t hash = 0;
    auto h1 = std::hash<T1>{}(pair.first);
    auto h2 = std::hash<T1>{}(pair.second);
    boost::hash_combine( hash, h1 );
    boost::hash_combine( hash, h2 );
    return hash;
  }
};


// Each dumpi_to_graph process will maintain one EventGraph object that contains
// all of the event sequence and channel map data from the possibly several 
// DUMPI trace files it is responsible for
class EventGraph
{
public:
  EventGraph( const Configuration& config,
              const std::unordered_map<int,Trace*> rank_to_trace,
              const std::unordered_map<int,CSMPI_Trace*> rank_to_csmpi_trace);
  
  // Function for applying scalar logical timestamps
  void apply_scalar_logical_clock();

  // Function to merge the partial edge lists into a single igraph object then
  // write it out to disk
  void merge_and_write();

  // Convenience functions for printing state
  void report_program_order_edges() const;
  void report_message_order_edges() const;


private:
  // Makes sure each dumpi_to_graph process has the same view of any 
  // user-defined communicators 
  CommunicatorManager exchange_user_defined_comm_data(); 
  
  // Makes sure each vertex has a unique ID
  void disambiguate_vertex_ids();

  // Makes sure all channels are specified in terms of global ranks
  void disambiguate_channel_maps();

  // Merges trace data from all trace ranks that this dumpi_to_graph process
  // is handling
  void merge_trace_data(); 

  void make_program_order_edges();
  void make_message_order_edges();
  void make_collective_edges();

  // Helpers for point-to-point message edge construction
  void exchange_local_message_matching_data();
  void exchange_remote_message_matching_data();
  void exchange_message_matching_data_for_communicator( int comm_id );


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

  // Built up from each trace rank's individual view of the communicator data
  CommunicatorManager comm_manager;

  // Data collected by merging corresponding data structures from all traces
  // this dumpi_to_graph is handling
  channel_map channel_to_send_seq;
  channel_map channel_to_recv_seq;
  std::unordered_map<size_t,Channel> vertex_id_to_channel;
  std::unordered_map<size_t,uint8_t> vertex_id_to_event_type;
  std::unordered_map<size_t,double> vertex_id_to_wall_time;
  std::unordered_map<size_t,int> vertex_id_to_pid; 

  //// Communicator properties
  //std::unordered_map<int,size_t> comm_to_size;
  //std::unordered_map<int,int> comm_to_parent;
  //std::unordered_map<int, std::unordered_map<int,std::pair<int,int>>> comm_to_rankcolorkey;
  //std::unordered_map<std::pair<int,int>, int, pair_hash> comm_rank_pair_to_global_rank;
};


#endif // D2G_EVENT_GRAPH_H
