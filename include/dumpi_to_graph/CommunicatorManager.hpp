#ifndef D2G_COMMUNICATOR_MANAGER_H
#define D2G_COMMUNICATOR_MANAGER_H

#include <unordered_map>

// DUMPI
#include "dumpi/common/constants.h"

class CommunicatorManager
{
public:
  CommunicatorManager() {}
  CommunicatorManager( size_t global_comm_size );
  int get_global_rank( int comm_id, int comm_rank ); 
  // This function broadcasts each dumpi_to_graph process's view of the 
  // communicators of the traced run and ensure that each dumpi_to_graph process
  // has the same complete mapping between global ranks and communicator ranks
  // prior to event graph construction
  void merge();
  // Helper functions
  void update_comm_to_parent( int new_comm_id, int parent_comm_id );
  void associate_rank_with_color( int comm_id, int global_rank, int color );
  void associate_rank_with_key( int comm_id, int global_rank, int key );
private:
  std::unordered_map<int,size_t> comm_to_size;
  std::unordered_map<int,int> comm_to_parent;

  std::unordered_map<int, std::unordered_map<int,int>> comm_to_rank_to_color;
  std::unordered_map<int, std::unordered_map<int,int>> comm_to_rank_to_key;
};


#endif // D2G_COMMUNICATOR_MANAGER_H
