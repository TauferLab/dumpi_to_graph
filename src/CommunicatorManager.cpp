#include "CommunicatorManager.hpp"

#include "mpi.h"

#include <unordered_map>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>

#include <functional>
#include <utility> 

// DUMPI
#include "dumpi/common/constants.h" // DUMPI_COMM_WORLD

// Internal 
#include "Channel.hpp"
#include "Utilities.hpp"

int CommunicatorManager::get_comm_rank( int comm_id, int global_rank )
{
  auto global_rank_to_comm_rank = comm_to_global_rank_to_comm_rank.at( comm_id );
  auto comm_rank = global_rank_to_comm_rank.at( global_rank );
  return comm_rank;
}

void CommunicatorManager::set_comm_to_translator( std::unordered_map<int,std::unordered_map<std::pair<int,std::vector<int>>,int,rank_seq_hash>> comm_to_translator ) 
{
  this->comm_to_translator = comm_to_translator;
}

std::unordered_map<int,std::unordered_map<std::pair<int,std::vector<int>>,int,rank_seq_hash>> CommunicatorManager::get_comm_to_translator() const
{
  return this->comm_to_translator;
}


// Translates a channel where either the src or dst ranks (but not both) are
// local to the communicator to a channel where both src and dst are global ranks
Channel CommunicatorManager::translate_channel( Channel channel, int known_global_rank )
{
  // Unpack channel
  auto src = channel.get_src();
  auto dst = channel.get_dst();
  auto tag = channel.get_tag();
  auto comm_id = channel.get_comm();

  // Early exit for global communicator
  if ( comm_id == DUMPI_COMM_WORLD ) {
    return channel;
  }
  
  // Communicator depth
  std::unordered_map<int,int> comm_to_idx;
  for ( auto kvp : comm_to_parent ) {
    int idx = 1;
    auto curr_comm = kvp.first;
    auto parent_comm = kvp.second;
    while ( parent_comm != DUMPI_COMM_WORLD ) {
      curr_comm = parent_comm;
      parent_comm = comm_to_parent.at( curr_comm );
      idx++;
    } 
    comm_to_idx.insert( { kvp.first, idx } );
  }
  auto comm_idx = comm_to_idx.at( comm_id ); 

  //// Get the translation map for this communicator
  //auto translator_search = this->comm_to_translator.find( comm_id );
  //if ( translator_search == this->comm_to_translator.end() ) {
  //  std::stringstream ss;
  //  ss << "Cannot find translator for comm: " << comm_id << std::endl;
  //  throw std::runtime_error( ss.str() );
  //}

  auto translator = this->comm_to_translator.at( comm_id );

  // Case 1: Source rank is global, destination rank is comm-local
  if ( known_global_rank == 0 ) {

    //if ( dst > this->comm_to_size.at( comm_id )-1 ) {
    //  return channel; 
    //}

    auto color_seq = this->get_color_seq( src );
    std::vector<int> color_subseq;
    for ( int i = 0; i < comm_idx; ++i ) {
      color_subseq.push_back( color_seq[i] );
    }
    auto key = std::make_pair( dst, color_subseq );
    
    //// Hack to deal with some dst ranks that are already correct
    //// FIXME: No idea why this is happening
    auto search = translator.find( key );
    if ( search == translator.end() ) {
      //std::stringstream ss;
      //ss << "Cannot translate dst comm rank: " << dst 
      //   << " in comm: " << comm_id
      //   << " to global rank" << std::endl;
      //throw std::runtime_error( ss.str() );
      return channel;
    }

    auto dst_global_rank = translator.at( key );
    Channel translated_channel( src, dst_global_rank, tag, comm_id );
    return translated_channel;
  }
  // Case 2: Destination rank is global, source rank is comm-local
  else if ( known_global_rank == 1 ) {
    auto color_seq = this->get_color_seq( dst );
    std::vector<int> color_subseq;
    for ( int i = 0; i < comm_idx; ++i ) {
      color_subseq.push_back( color_seq[i] );
    }
    auto key = std::make_pair( src, color_subseq );

    auto search = translator.find( key );
    if ( search == translator.end() ) {
      std::stringstream ss;
      ss << "Cannot translate src comm rank: " << src 
         << " in comm: " << comm_id
         << " to global rank" << std::endl;
      throw std::runtime_error( ss.str() );
    }

    auto src_global_rank = translator.at( key );
    Channel translated_channel( src_global_rank, dst, tag, comm_id );
    return translated_channel;
  }
  else {
    throw std::runtime_error("Invalid option for translate_channel");
  }
}


//int CommunicatorManager::sender_comm_rank_to_global_rank( Channel channel )
//{
//  // Unpack channel
//  auto sender_comm_rank     = channel.get_src();
//  auto receiver_global_rank = channel.get_dst();
//  auto comm_id              = channel.get_comm();
//
//  // Need to find out which process group this communication occured in, so
//  // we look up the "color" using the global rank of the sender
//  auto rank_to_color = comm_to_rank_to_color.at( comm_id );
//  auto color         = rank_to_color.at( receiver_global_rank );
//
//  // 
//  std::set<int> colors;
//  for ( auto kvp : rank_to_color ) {
//    auto color = kvp.second;
//    colors.insert( color );
//  }
//  int comm_rank_subrange_idx = 0;
//  int i = 0;
//  for ( auto curr_color : colors ) {
//    if ( curr_color == color ) {
//      comm_rank_subrange_idx = i;
//      break;
//    }
//    i++;
//  }
//  auto comm_size = comm_to_size.at( comm_id );
//  
//  auto global_rank_to_key = comm_to_rank_to_key.at( comm_id );
//
//  std::vector<int> keys;
//  for ( auto kvp : global_rank_to_key ) {
//    auto key = kvp.second;
//    keys.push_back( key );
//  }
//  std::sort( keys.begin(), keys.end() );
//
//  int key_idx = comm_size * comm_rank_subrange_idx + sender_comm_rank;
//  
//  auto key = keys[ key_idx ];
//
//  //std::cout << "Translating comm rank: " << receiver_comm_rank 
//  //          << " in comm: " << comm_id 
//  //          << " color: " << color
//  //          << " subrange size: " << comm_size
//  //          << " subrange idx: " << comm_rank_subrange_idx
//  //          << " key idx: " << key_idx
//  //          << " key: " << key
//  //          << std::endl;
//
//  std::unordered_map<int,int> key_to_global_rank;
//  for ( auto kvp : global_rank_to_key ) {
//    auto curr_key = kvp.second;
//    auto curr_rank = kvp.first;
//    auto curr_color = rank_to_color.at( curr_rank );
//    if ( curr_color == color ) {
//      key_to_global_rank.insert( { curr_key, curr_rank } );
//    }
//  }
//
//  auto sender_global_rank = key_to_global_rank.at( key );
//
//  //std::cout << "Translated comm rank: " << sender_comm_rank 
//  //          << " to global rank: " << sender_global_rank
//  //          << " in comm: " << comm_id
//  //          << std::endl;
//
//  return sender_global_rank; 
//}
//
//int CommunicatorManager::receiver_comm_rank_to_global_rank( Channel channel )
//{
//  // Unpack channel
//  auto sender_global_rank = channel.get_src();
//  auto receiver_comm_rank = channel.get_dst();
//  auto comm_id            = channel.get_comm();
//
//  // We need to determine how to traverse the communicator tree to get to the 
//  // right range of ranks. To do this, we need to first figure out which process
//  // subgroup this communication took place in, not just for this communicator,
//  // but also for all of its ancestors
//  std::vector<int> color_seq;
//  auto curr_comm_id = comm_id;
//  while ( curr_comm_id != DUMPI_COMM_WORLD ) {
//    // Look up color based on the one global rank we already have (the sender)
//    auto rank_to_color = comm_to_rank_to_color.at( curr_comm_id );
//    auto color = rank_to_color.at( sender_global_rank );
//    color_seq.push_back( color );
//    curr_comm_id = comm_to_parent.at( curr_comm_id );
//  }
//  
//
//
//  // Need to find out which process group this communication occured in, so
//  // we look up the "color" using the global rank of the sender
//  auto rank_to_color = comm_to_rank_to_color.at( comm_id );
//  auto color         = rank_to_color.at( sender_global_rank );
//
//
//  // 
//  std::set<int> colors;
//  for ( auto kvp : rank_to_color ) {
//    auto color = kvp.second;
//    colors.insert( color );
//  }
//  int comm_rank_subrange_idx = 0;
//  int i = 0;
//  for ( auto curr_color : colors ) {
//    if ( curr_color == color ) {
//      comm_rank_subrange_idx = i;
//      break;
//    }
//    i++;
//  }
//  auto comm_size = comm_to_size.at( comm_id );
//  
//  auto global_rank_to_key = comm_to_rank_to_key.at( comm_id );
//
//  std::vector<int> keys;
//  for ( auto kvp : global_rank_to_key ) {
//    auto key = kvp.second;
//    keys.push_back( key );
//  }
//  std::sort( keys.begin(), keys.end() );
//
//  int key_idx = comm_size * comm_rank_subrange_idx + receiver_comm_rank;
//  
//  auto key = keys[ key_idx ];
//
//  //std::cout << "Translating comm rank: " << receiver_comm_rank 
//  //          << " in comm: " << comm_id 
//  //          << " color: " << color
//  //          << " subrange size: " << comm_size
//  //          << " subrange idx: " << comm_rank_subrange_idx
//  //          << " key idx: " << key_idx
//  //          << " key: " << key
//  //          << std::endl;
//
//  std::unordered_map<int,int> key_to_global_rank;
//  for ( auto kvp : global_rank_to_key ) {
//    auto curr_key = kvp.second;
//    auto curr_rank = kvp.first;
//    auto curr_color = rank_to_color.at( curr_rank );
//    if ( curr_color == color ) {
//      key_to_global_rank.insert( { curr_key, curr_rank } );
//    }
//  }
//
//  auto global_rank = key_to_global_rank.at( key );
//
//  //std::cout << "Translated comm rank: " << receiver_comm_rank 
//  //          << " to global rank: " << global_rank
//  //          << " in comm: " << comm_id
//  //          << std::endl;
//  
//  return global_rank;
//}

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
    this->comm_to_global_rank_to_comm_rank = rhs.get_comm_to_global_rank_to_comm_rank();
    this->comm_to_translator = rhs.get_comm_to_translator();
    return *this;
  }
}

// Get the sequence of colors indicating which subgroups a global rank is part of
std::vector<int> CommunicatorManager::get_color_seq( int global_rank )
{
  std::vector<int> color_seq;
  for ( auto kvp : comm_to_rank_to_color ) {
    auto rank_to_color = kvp.second;
    auto color = rank_to_color.at( global_rank );
    color_seq.push_back( color );
  }
  return color_seq;
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
  
void CommunicatorManager::calculate_global_rank_to_comm_rank_mapping()
{

  int mpi_rc, my_rank;
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

  // Need this to know how much of the color sequence to match below
  std::unordered_map<int,int> comm_to_idx;
  for ( auto kvp : comm_to_parent ) {
    int idx = 1;
    auto curr_comm = kvp.first;
    auto parent_comm = kvp.second;
    while ( parent_comm != DUMPI_COMM_WORLD ) {
      curr_comm = parent_comm;
      parent_comm = comm_to_parent.at( curr_comm );
      idx++;
    } 
    comm_to_idx.insert( { kvp.first, idx } );
  }

  //if ( my_rank == 12 ) {
  //  for ( auto kvp : comm_to_idx ) {
  //    std::cout << "Comm: " << kvp.first
  //              << " Idx: " << kvp.second
  //              << std::endl;
  //  }
  //}

  //exit(0);


  // Proceed one communicator at a time
  for ( auto kvp : comm_to_rank_to_color ) {
    auto comm_id = kvp.first;

    std::unordered_map<int,int> global_rank_to_comm_rank;
    this->comm_to_global_rank_to_comm_rank.insert( { comm_id, global_rank_to_comm_rank } );

    auto rank_to_color = kvp.second;
    auto rank_to_key = comm_to_rank_to_key.at( comm_id );
    auto comm_size = comm_to_size.at( comm_id );
    for ( auto kvp2 : rank_to_key ) {
      auto global_rank = kvp2.first;
      auto key = kvp2.second;
      auto color_seq = this->get_color_seq( global_rank );
      std::vector<int> color_subseq;
      for ( int i=0; i<comm_to_idx.at( comm_id ); ++i ) {
        color_subseq.push_back( color_seq[i] );
      }
      std::vector<int> keys_in_group;
      keys_in_group.push_back( key );
      // get all the keys for other global ranks that ended up in the same 
      // process group (i.e., global ranks that have the same color sequence 
      // up to the current communicator)
      for ( auto kvp3 : rank_to_key ) {
        auto other_rank = kvp3.first;
        if ( other_rank != global_rank ) {
          auto other_key = kvp3.second;
          auto other_color_seq = this->get_color_seq( other_rank );
          std::vector<int> other_color_subseq;
          for ( int i=0; i<comm_to_idx.at( comm_id ); ++i ) {
            other_color_subseq.push_back( other_color_seq[i] );
          }
          if ( other_color_subseq == color_subseq ) {
            keys_in_group.push_back( other_key );
          }
        }
      }
      std::sort( keys_in_group.begin(), keys_in_group.end() );

      int comm_rank = 0;
      for ( auto current_key : keys_in_group ) {
        // Found comm rank of this global rank
        if ( current_key == key ) {
          //if ( my_rank == 12 ) {
          //  std::cout << "Global Rank: " << global_rank 
          //            << " Comm: " << comm_id
          //            << " Comm Rank: " << comm_rank
          //            << std::endl;
          //}
          this->comm_to_global_rank_to_comm_rank.at( comm_id ).insert( { global_rank, comm_rank } );
        }
        comm_rank++;
      }
    }
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
  for ( auto kvp : rhs_comm_manager.get_comm_to_global_rank_to_comm_rank() ) {
    auto comm_id = kvp.first;
    auto rhs_global_rank_to_comm_rank = kvp.second;
    auto comm_search = this->comm_to_global_rank_to_comm_rank.find( comm_id );
    if ( comm_search == this->comm_to_global_rank_to_comm_rank.end() ) {
      std::unordered_map<int,int> global_rank_to_comm_rank;
      for ( auto kvp2 : rhs_global_rank_to_comm_rank ) {
        auto global_rank = kvp2.first;
        auto comm_rank = kvp2.second;
        global_rank_to_comm_rank.insert( { global_rank, comm_rank } );
      }
      this->comm_to_global_rank_to_comm_rank.insert( { comm_id, global_rank_to_comm_rank } );
    }
    else {
      for ( auto kvp2 : rhs_global_rank_to_comm_rank ) {
        auto global_rank = kvp2.first;
        auto comm_rank = kvp2.second;
        comm_search->second.insert( { global_rank, comm_rank } );
      }
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

std::unordered_map<int,std::unordered_map<int,int>> CommunicatorManager::get_comm_to_global_rank_to_comm_rank() const
{
  return this->comm_to_global_rank_to_comm_rank;
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
    auto rank_to_comm_rank = comm_to_global_rank_to_comm_rank.at( comm_id );
    for ( auto kvp : rank_to_color ) {
      auto rank = kvp.first;
      auto color = kvp.second;
      auto key = rank_to_key.at( rank );
      auto comm_rank = rank_to_comm_rank.at( rank );
      std::cout << "\tGlobal Rank: " << rank 
                << ", Color: " << color
                << ", Key: " << key 
                << ", Comm Rank: " << comm_rank
                << std::endl;
    }
  }
  std::cout << std::endl;
  for ( auto comm_translator : comm_to_translator ) {
    auto comm_id = comm_translator.first;
    auto translator = comm_translator.second;
    std::cout << "Comm: " << comm_id << " Translator:" << std::endl;
    for ( auto kvp : translator ) {
      auto global_rank = kvp.second;
      auto comm_rank = kvp.first.first;
      auto color_seq = kvp.first.second;
      std::cout << "Comm. Rank: " << comm_rank 
                << " Color Subseq: ";
      for ( auto color : color_seq ) {
        std::cout << color << " ";
      }
      std::cout << "--> Global Rank: " << global_rank << std::endl;
    }
  }
}

