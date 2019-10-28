#ifndef D2G_COMMUNICATOR_MANAGER_H
#define D2G_COMMUNICATOR_MANAGER_H

#include <unordered_map>

// DUMPI
#include "dumpi/common/constants.h"

// Boost
#include "boost/functional/hash.hpp"
#include "boost/serialization/access.hpp"
#include "boost/serialization/unordered_map.hpp"

// Internal 
#include "Channel.hpp"

class CommunicatorManager
{
public:
  CommunicatorManager() {}
  CommunicatorManager( size_t global_comm_size );
  // Used during message matching data exchange to translate from the 
  // communicator-specific ranks observed in the trace data to the global ranks
  // we should actually send to
  int sender_global_rank_to_comm_rank( Channel channel ); 
  int receiver_comm_rank_to_global_rank( Channel channel ); 
  // Function to calculate communicator sizes once all dumpi_to_graph processes
  // have sufficient information to do so
  void calculate_comm_sizes(); 
  // Function to incorporate data from another CommunicatorManager into this one
  void aggregate( const CommunicatorManager& comm_manager );
  // Helper functions
  void update_comm_to_parent( int new_comm_id, int parent_comm_id );
  void associate_rank_with_color( int comm_id, int global_rank, int color );
  void associate_rank_with_key( int comm_id, int global_rank, int key );
  // Accessors for aggregation
  std::unordered_map<int,size_t> get_comm_to_size() const;
  std::unordered_map<int,int> get_comm_to_parent() const;
  std::unordered_map<int,std::unordered_map<int,int>> get_comm_to_rank_to_color() const;
  std::unordered_map<int,std::unordered_map<int,int>> get_comm_to_rank_to_key() const;
  // Also needed for aggregation
  CommunicatorManager operator=( const CommunicatorManager& rhs );
  // Convenience printer
  void print() const;
private:
  std::unordered_map<int,size_t> comm_to_size;
  std::unordered_map<int,int> comm_to_parent;

  std::unordered_map<int, std::unordered_map<int,int>> comm_to_rank_to_color;
  std::unordered_map<int, std::unordered_map<int,int>> comm_to_rank_to_key;
  
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize( Archive& ar, const unsigned int version ) 
  {
    ar & comm_to_size;
    ar & comm_to_parent;
    ar & comm_to_rank_to_color;
    ar & comm_to_rank_to_key;
  }
};

#endif // D2G_COMMUNICATOR_MANAGER_H
