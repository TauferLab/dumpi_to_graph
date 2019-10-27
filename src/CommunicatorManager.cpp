#include "CommunicatorManager.hpp"

#include "mpi.h"

#include <unordered_map>
#include <stdexcept>
#include <sstream>
#include <iostream>

// DUMPI
#include "dumpi/common/constants.h" // DUMPI_COMM_WORLD


CommunicatorManager::CommunicatorManager( size_t global_comm_size )
{
  // Set size of global communicator
  comm_to_size.insert( { DUMPI_COMM_WORLD, global_comm_size } );
}

CommunicatorManager CommunicatorManager::operator=( const CommunicatorManager& rhs )
{
  if ( this == &rhs ) {
    return *this;
  }
  else {
    this->comm_to_size = rhs.get_comm_to_size();
    this->comm_to_parent = rhs.get_comm_to_parent();
    this->comm_to_rank_to_color = rhs.get_comm_to_rank_to_color();
    this->comm_to_rank_to_key = rhs.get_comm_to_rank_to_key();
    return *this;
  }
}
  
void CommunicatorManager::calculate_comm_sizes()
{
  for ( auto kvp : comm_to_rank_to_color ) {
    auto comm_id = kvp.first;
    auto rank_to_color = kvp.second;
    std::set<int> unique_colors;
    for ( auto kvp2 : rank_to_color ) { 
      auto rank = kvp2.first;
      auto color = kvp2.second;
      unique_colors.insert( color );
    }
    int n_colors = unique_colors.size();
    auto parent_comm_id = comm_to_parent.at( comm_id );
    auto parent_size = comm_to_size.at( parent_comm_id );
    int comm_size = parent_size / n_colors;
    comm_to_size.insert( { comm_id, comm_size } );
  } 
}

void CommunicatorManager::aggregate( const CommunicatorManager& rhs_comm_manager )
{
  for ( auto kvp : rhs_comm_manager.get_comm_to_size() ) {
    auto comm_id = kvp.first;
    auto comm_size = kvp.second;
    this->comm_to_size.insert( { comm_id, comm_size } );
  }
  for ( auto kvp : rhs_comm_manager.get_comm_to_parent() ) {
    auto comm_id = kvp.first;
    auto parent_comm_id = kvp.second;
    this->update_comm_to_parent( comm_id, parent_comm_id );
  }
  for ( auto kvp : rhs_comm_manager.get_comm_to_rank_to_color() ) {
    auto comm_id = kvp.first;
    auto rank_to_color = kvp.second;
    for ( auto rank_color : rank_to_color ) {
      auto rank = rank_color.first;
      auto color = rank_color.second;
      this->associate_rank_with_color( comm_id, rank, color );
    }
  }
  for ( auto kvp : rhs_comm_manager.get_comm_to_rank_to_key() ) {
    auto comm_id = kvp.first;
    auto rank_to_key = kvp.second;
    for ( auto rank_key : rank_to_key ) {
      auto rank = rank_key.first;
      auto key = rank_key.second;
      this->associate_rank_with_key( comm_id, rank, key );
    }
  }
}

// Updates the mapping between user-defined communicators and their "parent" 
// communicators (e.g., as in MPI_Comm_split)
void CommunicatorManager::update_comm_to_parent( int new_comm_id, int parent_comm_id )
{
  auto search = comm_to_parent.find( new_comm_id );
  if ( search == comm_to_parent.end() ) {
    comm_to_parent.insert( { new_comm_id, parent_comm_id } );
  } else {
    if ( parent_comm_id != search->second ) {
      std::stringstream ss;
      ss << "Trying to replace parent communicator of communicator: " << new_comm_id
         << " with communicator: " << parent_comm_id
         << " (old parent communicator: " << search->second << ")" << std::endl;
      throw std::runtime_error( ss.str() );
    }
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
      //rank_to_color.insert( { global_rank, color } );
      comm_search->second.insert( { global_rank, color } );
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
      //rank_to_key.insert( { global_rank, key } );
      comm_search->second.insert( { global_rank, key } );
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

// Accessors needed for aggregation
std::unordered_map<int,size_t> CommunicatorManager::get_comm_to_size() const
{
  return this->comm_to_size;
}

std::unordered_map<int,int> CommunicatorManager::get_comm_to_parent() const
{
  return this->comm_to_parent;
}

std::unordered_map<int,std::unordered_map<int,int>> CommunicatorManager::get_comm_to_rank_to_color() const
{
  return this->comm_to_rank_to_color;
}

std::unordered_map<int,std::unordered_map<int,int>> CommunicatorManager::get_comm_to_rank_to_key() const
{
  return this->comm_to_rank_to_key;
}

// Convenience printing function
void CommunicatorManager::print() const
{
  std::vector<int> user_defined_comm_ids;
  for ( auto kvp : comm_to_parent ) {
    user_defined_comm_ids.push_back( kvp.first );
  }

  for ( auto kvp : comm_to_size ) {
    auto comm_id = kvp.first;
    auto comm_size = kvp.second;
    std::cout << "Comm ID: " << comm_id 
              << " - Comm Size: " << comm_size << std::endl;
  }

  for ( auto kvp : comm_to_parent ) {
    auto comm_id = kvp.first;
    auto parent_comm_id = kvp.second;
    std::cout << "Comm ID: " << comm_id 
              << " - Parent Comm ID: " << parent_comm_id << std::endl;
  }

  for ( auto comm_id : user_defined_comm_ids ) {
    std::cout << "Comm ID: " << comm_id << std::endl;
    auto rank_to_color = comm_to_rank_to_color.at( comm_id );    
    auto rank_to_key = comm_to_rank_to_key.at( comm_id );    
    for ( auto kvp : rank_to_color ) {
      auto rank = kvp.first;
      auto color = kvp.second;
      auto key = rank_to_key.at( rank );
      std::cout << "\tRank: " << rank 
                << ", Color: " << color
                << ", Key: " << key << std::endl;
    }
  }

  //for ( auto kvp : comm_to_rank_to_color ) {
  //  auto comm_id = kvp.first;
  //  auto rank_to_color = kvp.second;
  //  std::cout << "Comm ID: " << comm_id << std::endl;
  //  for ( auto rank_color : rank_to_color ) {
  //    auto rank = rank_color.first;
  //    auto color = rank_color.second;
  //    std::cout << "\tRank: " << rank << ", Color: " << color << std::endl;
  //  }
  //}
  //for ( auto kvp : comm_to_rank_to_key ) {
  //  auto comm_id = kvp.first;
  //  auto rank_to_key = kvp.second;
  //  std::cout << "Comm ID: " << comm_id << std::endl;
  //  for ( auto rank_key : rank_to_key ) {
  //    auto rank = rank_key.first;
  //    auto key = rank_key.second;
  //    std::cout << "\tRank: " << rank << ", Key: " << key << std::endl;
  //  }
  //}
}

