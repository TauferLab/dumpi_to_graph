#include "CommunicatorManager.hpp"

#include "mpi.h"

#include <unordered_map>
#include <stdexcept>
#include <sstream>

// DUMPI
#include "dumpi/common/constants.h" // DUMPI_COMM_WORLD

CommunicatorManager::CommunicatorManager( size_t global_comm_size )
{
  // Set size of global communicator
  comm_to_size.insert( { DUMPI_COMM_WORLD, global_comm_size } );
}

// Unifies view of communicator data across all dumpi_to_graph processes
void CommunicatorManager::merge()
{
  
}


// Updates the mapping between user-defined communicators and their "parent" 
// communicators (e.g., as in MPI_Comm_split)
void CommunicatorManager::update_comm_to_parent( int new_comm_id, int parent_comm_id )
{
  auto search = comm_to_parent.find( new_comm_id );
  if ( search == comm_to_parent.end() ) {
    comm_to_parent.insert( { new_comm_id, parent_comm_id } );
  } else {
    std::stringstream ss;
    ss << "Trying to replace parent communicator of communicator: " << new_comm_id
       << " with communicator: " << parent_comm_id
       << " (old parent communicator: " << search->second << ")" << std::endl;
    throw std::runtime_error( ss.str() );
  }
}

// Associates a rank to a "color" in a user-defined communicator.
// This is MPI's terminology for (roughly) process subgroup and is needed to 
// disambiguate between processes in a communicator that have, e.g., the same 
// rank in their "half" of the "same" communicator
void CommunicatorManager::associate_rank_with_color( int comm_id, 
                                                     int global_rank,
                                                     int color )
{
  auto comm_search = comm_to_rank_to_color.find( comm_id );
  // If communicator has not been encountered yet, create a new rank_to_color
  // mapping and associate it with the communicator ID
  if ( comm_search == comm_to_rank_to_color.end() ) {
    std::unordered_map<int,int> rank_to_color;
    rank_to_color.insert( { global_rank, color } );
    comm_to_rank_to_color.insert( { comm_id, rank_to_color } );
  } 
  // Otherwise just update the rank_to_color mapping
  else {
    auto rank_to_color = comm_search->second;
    auto rank_search = rank_to_color.find( global_rank );
    // Check that rank is not already mapped to a color
    if ( rank_search == rank_to_color.end() ) {
      rank_to_color.insert( { global_rank, color } );
    } else {
      std::stringstream ss;
      ss << "Trying to re-map rank: " << global_rank
         << " in communicator: " << comm_id
         << " from color: " << rank_search->second
         << " to color: " << color << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}

// Associates a rank to a "key" in a user-defined communicator.
void CommunicatorManager::associate_rank_with_key( int comm_id, 
                                                   int global_rank,
                                                   int key )
{
  auto comm_search = comm_to_rank_to_key.find( comm_id );
  // If communicator has not been encountered yet, create a new rank_to_key
  // mapping and associate it with the communicator ID
  if ( comm_search == comm_to_rank_to_key.end() ) {
    std::unordered_map<int,int> rank_to_key;
    rank_to_key.insert( { global_rank, key } );
    comm_to_rank_to_key.insert( { comm_id, rank_to_key } );
  } 
  // Otherwise just update the rank_to_key mapping
  else {
    auto rank_to_key = comm_search->second;
    auto rank_search = rank_to_key.find( global_rank );
    // Check that rank is not already mapped to a key
    if ( rank_search == rank_to_key.end() ) {
      rank_to_key.insert( { global_rank, key } );
    } else {
      std::stringstream ss;
      ss << "Trying to re-map rank: " << global_rank
         << " in communicator: " << comm_id
         << " from key: " << rank_search->second
         << " to key: " << key << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}
