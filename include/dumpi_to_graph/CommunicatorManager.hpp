#ifndef D2G_COMMUNICATOR_MANAGER_H
#define D2G_COMMUNICATOR_MANAGER_H

#include <unordered_map>
#include <map>

// DUMPI
#include "dumpi/common/constants.h"

// Boost
#include "boost/functional/hash.hpp"
#include "boost/serialization/access.hpp"
#include "boost/serialization/unordered_map.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/vector.hpp"

// Internal 
#include "Channel.hpp"
#include "Utilities.hpp"

class CommunicatorManager
{
public:
  CommunicatorManager() {}
  CommunicatorManager( size_t global_comm_size );
 
  // 
  Channel translate_channel( Channel channel, int known_global_rank );

  // translator
  void set_comm_to_translator( std::unordered_map<int,std::unordered_map<std::pair<int,std::vector<int>>,int,rank_seq_hash>> comm_to_translator );
  std::unordered_map<int,std::unordered_map<std::pair<int,std::vector<int>>,int,rank_seq_hash>> get_comm_to_translator() const;

  // 
  std::vector<int> get_color_seq( int global_rank );
  int get_comm_rank( int comm_id, int global_rank );
  
  // Used during message matching data exchange to translate from the 
  // communicator-specific ranks observed in the trace data to the global ranks
  // we should actually send to
  //int sender_comm_rank_to_global_rank( Channel channel ); 
  //int receiver_comm_rank_to_global_rank( Channel channel ); 

  // Function to associate each communicator with its colors
  // Function to calculate communicator sizes once all dumpi_to_graph processes
  // have sufficient information to do so
  void calculate_comm_sizes(); 
  // Function to build mapping from global ranks to comm ranks
  void calculate_global_rank_to_comm_rank_mapping();
  // Function to incorporate data from another CommunicatorManager into this one
  void aggregate( const CommunicatorManager& comm_manager );
  // Helper functions
  void update_comm_to_parent( int new_comm_id, int parent_comm_id );
  void associate_rank_with_color( int comm_id, int global_rank, int color );
  void associate_rank_with_key( int comm_id, int global_rank, int key );
  // Accessors for aggregation
  std::unordered_map<int,size_t> get_comm_to_size() const;
  std::unordered_map<int,int> get_comm_to_parent() const;
  
  //std::unordered_map<int,std::unordered_map<int,int>> get_comm_to_rank_to_color() const;
  std::map<int,std::unordered_map<int,int>> get_comm_to_rank_to_color() const;
  
  
  std::unordered_map<int,std::unordered_map<int,int>> get_comm_to_rank_to_key() const;
  std::unordered_map<int,std::unordered_map<int,int>> get_comm_to_global_rank_to_comm_rank() const;
  // Also needed for aggregation
  CommunicatorManager operator=( const CommunicatorManager& rhs );
  // Convenience printer
  void print() const;


private:
  std::unordered_map<int, std::unordered_map<std::pair<int,std::vector<int>>,int,rank_seq_hash>> comm_to_translator;

  std::unordered_map<int,size_t> comm_to_size;
  std::unordered_map<int,int> comm_to_parent;

  //std::unordered_map<int, std::unordered_map<int,int>> comm_to_rank_to_color;
  std::map<int, std::unordered_map<int,int>> comm_to_rank_to_color;

  std::unordered_map<int, std::unordered_map<int,int>> comm_to_rank_to_key;

  std::unordered_map<int, std::unordered_map<int,int>> comm_to_global_rank_to_comm_rank;

  friend class boost::serialization::access;
  template <typename Archive>
  void serialize( Archive& ar, const unsigned int version ) 
  {
    ar & comm_to_size;
    ar & comm_to_parent;
    ar & comm_to_rank_to_color;
    ar & comm_to_rank_to_key;
    ar & comm_to_global_rank_to_comm_rank;
    ar & comm_to_translator;
  }
};

#endif // D2G_COMMUNICATOR_MANAGER_H
