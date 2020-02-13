#include "EventGraph.hpp"

#include <mpi.h>

#include <utility> // make_pair
#include <algorithm> // max
#include <cstdint>
#include <climits>
#include <cstdio>
#include <cinttypes>
#include <unordered_map>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <chrono>

// Boost
#include "boost/serialization/unordered_map.hpp" 
#include "boost/mpi.hpp"

// DUMPI
#include "dumpi/common/argtypes.h"
#include "dumpi/common/constants.h"

// Igraph
#include "igraph/igraph.h"

// Internal
#include "Logging.hpp"
#include "Debug.hpp"
#include "CommunicatorManager.hpp" 
#include "Utilities.hpp"
#include "CSMPI_Trace.hpp"
#include "Trace.hpp"

// Super gross workaround for using size_t for scalar logical clock
#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "Could not determine right MPI datatype for size_t"
#endif


// Constructor must do the following:
// - Share user-defined communicator information 
// - Disambiguate all vertex IDs from traces that this dumpi_to_graph process is
//   handling
// - Merge all mappings from channels to vertex ID sequences and mappings from 
//   vertex IDs to channels
// - Construct a vertex set for this dumpi_to_graph process's (partial) view of 
//   the global event graph 
// - Construct a corresponding edge set consisting of all of the program order
//   edges for trace ranks under this process's management, and all message
//   order edges whose sender vertex is in the program order of a trace rank 
//   under this process's management
EventGraph::EventGraph( const Configuration& config,
                        const std::unordered_map<int,Trace*> rank_to_trace,
                        const std::unordered_map<int,CSMPI_Trace*> rank_to_csmpi_trace )
{
  // Establish MPI context
  // Boost is used here for broadcasting std::unordered_maps 
  // (rather than writing our own packing functions)
  int mpi_rc, global_rank, n_procs; 
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &global_rank );
  mpi_rc = MPI_Comm_size( MPI_COMM_WORLD, &n_procs );
  boost::mpi::communicator comm_world;

#ifdef REPORT_PROGRESS
  std::cout << "Rank: " << global_rank 
            << " starting event graph construction" << std::endl;
#endif

  // Set members
  this->config = config;
  this->rank_to_trace = rank_to_trace;
  this->rank_to_csmpi_trace = rank_to_csmpi_trace;
  this->comm_manager = exchange_user_defined_comm_data();

  disambiguate_vertex_ids(); 
  comm_world.barrier();
#ifdef REPORT_PROGRESS
  std::cout << "Rank: " << global_rank 
            << " vertex IDs disambiguated" << std::endl;
#endif


  merge_trace_data();
  comm_world.barrier();
#ifdef REPORT_PROGRESS
  std::cout << "Rank: " << global_rank 
            << " trace data merged" << std::endl;
#endif


#ifdef REPORT_PROGRESS_VERBOSE
  // Sanity check user-defined comm data
  if ( global_rank == REPORTING_RANK ) {
    this->comm_manager.print(); 
  }
#endif

  disambiguate_channel_maps();
  comm_world.barrier();
#ifdef REPORT_PROGRESS
  std::cout << "Rank: " << global_rank 
            << " channel maps disambiguated" << std::endl;
#endif

#ifdef REPORT_PROGRESS_VERBOSE
  for ( auto kvp : channel_to_send_seq ) {
    std::cout << "Rank: " << global_rank
              << " Channel: " << kvp.first
              << " # Sends: " << kvp.second.size()
              << std::endl;
  }
  for ( auto kvp : channel_to_recv_seq ) {
    std::cout << "Rank: " << global_rank
              << " Channel: " << kvp.first
              << " # Recvs: " << kvp.second.size()
              << std::endl;
  }
  comm_world.barrier();
#endif

  // Construct edges
  this->make_program_order_edges();
  comm_world.barrier();

#ifdef PRINT_PROGRESS
  std::cout << "Rank: " << global_rank 
            << " constructed program order edges" << std::endl;
#endif

  this->make_message_order_edges();
  comm_world.barrier();

#ifdef PRINT_PROGRESS
  std::cout << "Rank: " << global_rank 
            << " constructed message order edges" << std::endl;
#endif

  this->make_collective_edges();
#ifdef PRINT_PROGRESS
  std::cout << "Rank: " << global_rank 
            << " constructed collectives edges" << std::endl;
#endif
}

CommunicatorManager EventGraph::exchange_user_defined_comm_data()
{
  // Establish MPI context
  // Boost is used here for broadcasting std::unordered_maps 
  // (rather than writing our own packing functions)
  int mpi_rc, global_rank, n_procs; 
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &global_rank );
  mpi_rc = MPI_Comm_size( MPI_COMM_WORLD, &n_procs );
  boost::mpi::communicator comm_world;

  // Since this dumpi_to_graph process may be managing multiple trace ranks, we
  // need to first merge all of their communicator data locally, then merge 
  // across dumpi_to_graph processes, then finally broadcast the unified view of
  // communicators 
  CommunicatorManager comm_manager;
  // Aggregate from all locally held comm_managers
  for ( auto kvp : rank_to_trace ) {
    auto trace_rank = kvp.first;
    auto trace_ptr = kvp.second;
    auto trace_comm_manager = trace_ptr->get_comm_manager();
    comm_manager.aggregate( trace_comm_manager );
  }
 
#ifdef REPORT_PROGRESS_VERBOSE
  std::cout << "Rank: " << global_rank << " locally-held comm. managers aggregated" << std::endl;
#endif

  // Collect all remotely held comm_managers
  std::vector<CommunicatorManager> remote_comm_managers;
  for ( int rank=0; rank < n_procs; rank++ ) {
    // If it's my rank, broadcast my comm manager
    CommunicatorManager payload;
    if ( global_rank == rank ) {
      payload = comm_manager; 
      boost::mpi::broadcast( comm_world, payload, rank );
    } 
    // If not, receive and aggregate
    else {
      boost::mpi::broadcast( comm_world, payload, rank );
      remote_comm_managers.push_back( payload );     
    }
  }

  // Aggregate the data held in the received comm managers
  for ( auto cm : remote_comm_managers ) {
    comm_manager.aggregate( cm );
  }

#ifdef REPORT_PROGRESS_VERBOSE
  std::cout << "Rank: " << global_rank << " remotely-held comm. managers aggregated" << std::endl;
#endif

  // Now each dumpi_to_graph process has a CommunicatorManager with enough data
  // to calculate the size each communicator. Specifically, the size of a 
  // communicator is the size of its parent divided by the number of distinct 
  // colors in the MPI_Comm_split call that constructed it
  comm_manager.calculate_comm_sizes(); 

#ifdef REPORT_PROGRESS_VERBOSE
  std::cout << "Rank: " << global_rank << " comm sizes calculated" << std::endl;
#endif

  comm_manager.calculate_global_rank_to_comm_rank_mapping();

#ifdef REPORT_PROGRESS_VERBOSE
  std::cout << "Rank: " << global_rank << " global-rank to comm-rank mapping built" << std::endl;
#endif

  comm_world.barrier(); 

  // Maps from a pair of (comm_rank, color_seq) to global_rank
  std::unordered_map<int,std::unordered_map<std::pair<int,std::vector<int>>,int,rank_seq_hash>> comm_to_translator;

  auto comm_to_rank_to_color = comm_manager.get_comm_to_rank_to_color();
  auto comm_to_rank_to_key = comm_manager.get_comm_to_rank_to_key();
  auto comm_to_parent = comm_manager.get_comm_to_parent();

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

  for ( auto kvp : comm_to_rank_to_color ) {

    // Get communicator we're building a translator for
    auto comm_id = kvp.first;

    // Initialize a translator for this communicator
    std::unordered_map< std::pair<int,std::vector<int>>, int, rank_seq_hash> translator;

    // Get the depth of this comm in the comm hierarchy. This will determine 
    // how much of the color sequence of global rank to use in constructing 
    // the key in the translator
    auto comm_idx = comm_to_idx.at( comm_id );
    
    // Get mappings from global ranks to colors and keys for this communicator
    auto global_rank_to_color = kvp.second;
    auto global_rank_to_key = comm_to_rank_to_key.at( comm_id );
    
    // Loop over global ranks
    for ( auto kvp2 : global_rank_to_color ) {
      auto global_rank = kvp2.first;

      // Translate to communicator rank
      auto comm_rank = comm_manager.get_comm_rank( comm_id, global_rank );

      // Get the full color sequence associated with the global rank
      auto color_seq = comm_manager.get_color_seq( global_rank );

      // Take the relevant subsequence of it
      std::vector<int> color_subseq;
      for ( int i = 0; i < comm_idx; ++i ) {
        color_subseq.push_back( color_seq[i] );
      }

      auto key = std::make_pair( comm_rank, color_subseq );
      translator.insert( { key, global_rank } );
    }

    // Associate this translator with its communicator
    comm_to_translator.insert( { comm_id, translator } );
  }

  comm_manager.set_comm_to_translator( comm_to_translator );

  return comm_manager; 
}


// Goes through all data structures that involve a channel and makes sure that 
// each channel's src and dst members refer to global ranks, rather than ranks
// local to a user-defined communicator. Otherwise, message-edge construction
// breaks/deadlocks
void EventGraph::disambiguate_channel_maps()
{
  std::unordered_map<Channel,std::vector<size_t>,ChannelHash> translated_channel_to_send_seq;
  std::unordered_map<Channel,std::vector<size_t>,ChannelHash> translated_channel_to_recv_seq;
  
  for ( auto kvp : this->channel_to_send_seq ) {
    auto channel = kvp.first;
    auto send_seq = kvp.second;
    auto translated_channel = this->comm_manager.translate_channel( channel, 0 );
    auto send_channel_search = translated_channel_to_send_seq.find( translated_channel );
    if ( send_channel_search == translated_channel_to_send_seq.end() ) {
      translated_channel_to_send_seq.insert( { translated_channel, send_seq } );
    } 
    else {
      std::stringstream ss;
      ss << "Trying to overwrite send channel: " << translated_channel << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
  
  for ( auto kvp : this->channel_to_recv_seq ) {
    auto channel = kvp.first;
    auto recv_seq = kvp.second;
    auto translated_channel = this->comm_manager.translate_channel( channel, 1 );
    auto recv_channel_search = translated_channel_to_recv_seq.find( translated_channel );
    if ( recv_channel_search == translated_channel_to_recv_seq.end() ) {
      translated_channel_to_recv_seq.insert( { translated_channel, recv_seq } );
    } 
    else {
      std::stringstream ss;
      ss << "Trying to overwrite recv channel: " << translated_channel << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }

  this->channel_to_send_seq = translated_channel_to_send_seq;
  this->channel_to_recv_seq = translated_channel_to_recv_seq;

  // Finally, update the channels for the vertex_id to channel map
  std::unordered_map<size_t, Channel> vertex_id_to_translated_channel;
  for ( auto kvp : this->vertex_id_to_channel ) {
    auto vertex_id = kvp.first;
    auto channel = kvp.second;
    // Determine if this is a send or a receive. If it's a send, we know the 
    // src member of the channel is the global rank. If it's a recv, we know the
    // dst member of the channel is the global rank.
    auto event_type = this->vertex_id_to_event_type.at( vertex_id );
    if ( event_type == 0 ) {
      auto translated_channel = this->comm_manager.translate_channel( channel, 0 );
      vertex_id_to_translated_channel.insert( { vertex_id, translated_channel } );
    }
    else if ( event_type == 1 ) {
      auto translated_channel = this->comm_manager.translate_channel( channel, 1 );
      vertex_id_to_translated_channel.insert( { vertex_id, translated_channel } );
    }
    else {
      std::stringstream ss;
      ss << "Trying to translate channel for a non-send and non-recv event" << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
  this->vertex_id_to_channel = vertex_id_to_translated_channel;
}





void EventGraph::merge_trace_data()
{
  // Merge data from the various traces this dumpi_to_graph process is handling
  for ( auto kvp : this->rank_to_trace ) {
    // Merge in channel map data
    auto trace_channel_to_send_seq = kvp.second->get_channel_to_send_seq();
    auto trace_channel_to_recv_seq = kvp.second->get_channel_to_recv_seq();
    for ( auto kvp2 : trace_channel_to_send_seq ) {
      auto channel = kvp2.first;
      auto send_seq = kvp2.second;
      this->channel_to_send_seq.insert( { channel, send_seq } );
    }
    for ( auto kvp2 : trace_channel_to_recv_seq ) {
      auto channel = kvp2.first;
      auto recv_seq = kvp2.second;
      this->channel_to_recv_seq.insert( { channel, recv_seq } );
    }
    // Merge in vertex_id_to_channel data
    auto trace_vertex_id_to_channel = kvp.second->get_vertex_id_to_channel();
    for ( auto kvp2 : trace_vertex_id_to_channel ) {
      auto vertex_id = kvp2.first;
      auto channel = kvp2.second;
      this->vertex_id_to_channel.insert( { vertex_id, channel } );
    }
    // Merge in event type, pid, wall-time, and generating function call  data
    auto event_seq = kvp.second->get_event_seq();
    auto wall_time_seq = kvp.second->get_wall_time_seq(); 
    auto fn_call_seq = kvp.second->get_mpi_fn_seq();
    auto perf_counter_seq = kvp.second->get_perf_counter_seq();
    
    size_t n_events = event_seq.size();
    size_t initial_vertex_id = kvp.second->get_initial_vertex_id();
    for ( int i=0; i<n_events; ++i ) {
      uint8_t event_type = event_seq[i];
      double wall_time = wall_time_seq[i];
      size_t vertex_id = initial_vertex_id + i;
      int pid = kvp.first; 
      auto fn_idx_pair = fn_call_seq[i];
      this->vertex_ids.push_back( vertex_id );
      this->vertex_id_to_event_type.insert( { vertex_id, event_type } );
      this->vertex_id_to_wall_time.insert( { vertex_id, wall_time } );
      this->vertex_id_to_pid.insert( { vertex_id, pid } );
      this->vertex_id_to_fn_call.insert( { vertex_id, fn_idx_pair } );
      if(this->config.get_papi_flag()){
        std::string perf_counter = perf_counter_seq[i];
        this->vertex_id_to_papi.insert( {vertex_id, perf_counter } );
      }
    }
  }

  // Merge data from CSMPI traces if available
  if ( this->config.has_csmpi() ) {
#ifdef REPORT_PROGRESS
    int mpi_rc, global_rank;
    mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &global_rank );
    std::cout << "Rank: " << global_rank 
              << " merging CSMPI data, " 
              << "# vertex IDs to assign callstacks to: " << vertex_id_to_fn_call.size() << std::endl;
#endif

    for ( auto kvp : vertex_id_to_fn_call ) {
      auto vertex_id = kvp.first;
      auto fn = kvp.second.first;
      auto call_idx = kvp.second.second;
      auto pid = this->vertex_id_to_pid.at( vertex_id );
      auto search = this->rank_to_csmpi_trace.find( pid );
      auto csmpi_trace_ptr = this->rank_to_csmpi_trace.at( pid );
      auto callstack = csmpi_trace_ptr->lookup_callstack( fn, call_idx );
      this->vertex_id_to_callstack.insert( { vertex_id, callstack } );
    } 
  }
}

void EventGraph::disambiguate_vertex_ids()
{
  // Establish MPI context
  // Boost is used here for broadcasting std::unordered_maps 
  // (rather than writing our own packing functions)
  int mpi_rc, global_rank, n_procs; 
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &global_rank );
  mpi_rc = MPI_Comm_size( MPI_COMM_WORLD, &n_procs );
  boost::mpi::communicator comm_world;
  
  // First determine vertex counts for each trace's event sequence
  std::unordered_map<int,size_t> trace_rank_to_vertex_count;
  for ( auto kvp : this->rank_to_trace ) {
    size_t vertex_count = kvp.second->get_final_vertex_id() + 1;
    trace_rank_to_vertex_count.insert( { kvp.first, vertex_count } );
  }

  // Broadcast vertex counts
  std::unordered_map<int,size_t> all_trace_rank_to_vertex_count;
  for ( int i=0; i<comm_world.size(); ++i) {
    std::unordered_map<int,size_t> payload;
    if ( comm_world.rank() == i ) {
      payload = trace_rank_to_vertex_count;
      boost::mpi::broadcast( comm_world, payload, i );
    } else {
      boost::mpi::broadcast( comm_world, payload, i );
    }
    for ( auto kvp : payload ) {
      all_trace_rank_to_vertex_count.insert( { kvp.first, kvp.second } );
    }
  }

  // Calculate vertex ID offsets
  std::unordered_map<int,size_t> trace_rank_to_offset;
  for ( auto kvp : all_trace_rank_to_vertex_count ) {
    size_t offset = 0;
    for ( int i=0; i<kvp.first; ++i ) {
      offset += all_trace_rank_to_vertex_count.at( i );
    }
    trace_rank_to_offset.insert( { kvp.first, offset } );
  }
  // Apply offsets
  for ( auto kvp : this->rank_to_trace ) {
    size_t offset = trace_rank_to_offset.at( kvp.first );
    kvp.second->apply_vertex_id_offset( offset );
  }
}





//////////////// Top-level functions for edge construction /////////////////////

// Constructs all program order edges 
void EventGraph::make_program_order_edges()
{
  // Loop over the event sequences of the traces and make edges for consecutive
  // pairs of events from the init event to the finalize event
  for ( auto kvp : this->rank_to_trace ) {
    auto event_seq = kvp.second->get_event_seq();
    size_t n_events = event_seq.size();
    size_t initial_vertex_id = kvp.second->get_initial_vertex_id();
    size_t final_vertex_id = kvp.second->get_final_vertex_id();
    for ( int i=1; i<n_events ; ++i ) {
      int src_vertex_id = initial_vertex_id + i - 1;
      int dst_vertex_id = initial_vertex_id + i;
      auto edge = std::make_pair( src_vertex_id, dst_vertex_id );
      this->program_order_edges.push_back( edge );
    }
  }
}

// Top level function invoked in EventGraph constructor for building all message 
// edges
void EventGraph::make_message_order_edges()
{
  this->exchange_local_message_matching_data();
  this->exchange_remote_message_matching_data();
}

// TODO
void EventGraph::make_collective_edges()
{
}

//////////////// Helper functions for message edge construction ////////////////

// Builds message edges between send and recv vertices that are owned by the 
// same dumpi_to_graph process
void EventGraph::exchange_local_message_matching_data()
{
  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  // We need to iterate over this local copy so that we don't invalidate when erasing
  auto channel_to_send_seq_copy = this->channel_to_send_seq;
  for ( auto kvp : channel_to_send_seq_copy ) {
    // The rank of the receiving process in the traced application
    int dst = kvp.first.get_dst();
    // The dumpi_to_graph process managing that rank
    int owning_rank = this->config.lookup_owning_rank( dst );
    // If this dumpi_to_graph process also owns the matching receive sequence,
    // we can determine the set of message edges for this channel without any
    // communication between dumpi_to_graph processes
    if ( owning_rank == rank ) {  
      // Get the matching send and recv sequences for this channel
      auto send_seq = kvp.second;
      auto recv_seq = this->channel_to_recv_seq.at( kvp.first );
#ifdef SANITY_CHECK
      assert( send_seq.size() == recv_seq.size() );
      assert( contains_no_invalid_vertex_ids( send_seq ) );
      assert( contains_no_invalid_vertex_ids( recv_seq ) );
#endif
      // Make the message edges for this channel
      for ( int i=0; i<send_seq.size(); ++i ) {
        auto edge = std::make_pair( send_seq[i], recv_seq[i] );
        this->message_order_edges.push_back( edge );
      }
      // Remove the channel from both maps
      // This simplifies the code for the exchange of message matching data for
      // send and recv sequences held on distinct dumpi_to_graph processes
      this->channel_to_send_seq.erase( kvp.first );
      this->channel_to_recv_seq.erase( kvp.first );
    }
  }
}

// Builds message edges between send and recv vertices that are owened by 
// different dumpi_to_graph processes. 
void EventGraph::exchange_remote_message_matching_data()
{
  // Complete vertex ID exchange for each communicator
  for ( auto kvp : this->comm_manager.get_comm_to_size() ) {
    auto comm_id = kvp.first;
    exchange_message_matching_data_for_communicator( comm_id );
  }
}

// Builds message edges for messages that occurred within a specified communicator  
void EventGraph::exchange_message_matching_data_for_communicator( int current_comm_id )
{
  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // First we exchange all of the recv vertex IDs for messages in the global
  // communicator
  int n_recv_reqs = 0;
  for ( auto kvp : channel_to_send_seq ) {
    auto comm_id = kvp.first.get_comm();
    if ( comm_id == current_comm_id ) {
      n_recv_reqs += 1; 
    }
  }
  // Allocate one receive request per channel that we will receive a list of 
  // recv vertex IDs from
  MPI_Request recv_reqs[ n_recv_reqs ];
  // This will keep track of which channel we're handling
  int recv_req_idx = 0;
  // We are going to receive all of the vertex IDs into a single contiugous array
  // To do this, we need to keep track of how many vertex IDs in total we will 
  // receive, and an offset into the array for each channel since the actual 
  // order of arrival of this data will be unspecified
  std::unordered_map<Channel,int,ChannelHash> channel_to_offset;
  int n_recv_vertex_ids = 0;
  for ( auto kvp : channel_to_send_seq ) {
    auto channel = kvp.first;
    auto n_messages = kvp.second.size();
    auto comm_id = channel.get_comm();
    if ( comm_id == current_comm_id ) {
      channel_to_offset.insert( { channel, n_recv_vertex_ids } );
      n_recv_vertex_ids += n_messages;
    }
  }
  // Allocate the buffer for all of the recv vertex IDs
  int recv_buffer[ n_recv_vertex_ids ];
  // Loop over our send sequences and post a receive for each matching recv
  // sequence at the appropriate offset
  for ( auto kvp : this->channel_to_send_seq ) {
    int comm_id = kvp.first.get_comm();
    if ( comm_id == current_comm_id ) {
      int src = kvp.first.get_dst();
      int tag = kvp.first.get_tag();
      int recv_buffer_offset = channel_to_offset.at( kvp.first );
      int n_elements = kvp.second.size();
      mpi_rc = MPI_Irecv( &recv_buffer[ recv_buffer_offset ],
                          n_elements,
                          MPI_INT, 
                          src,
                          tag,
                          MPI_COMM_WORLD,
                          &recv_reqs[ recv_req_idx ] );
      recv_req_idx++;
    }
  }
  // Loop over our recv sequences and send each to the dumpi_to_graph process
  // with the matching send sequence
  for ( auto kvp : this->channel_to_recv_seq ) {
    int comm_id = kvp.first.get_comm();
    if ( comm_id == current_comm_id ) {
      int dst = kvp.first.get_src(); 
      int tag = kvp.first.get_tag();
      int n_elements = kvp.second.size();
#ifdef SANITY_CHECK
      assert( contains_no_invalid_vertex_ids( kvp.second ) );
#endif
      int send_buffer[ n_elements ];
      for ( int i=0; i<n_elements; ++i ) {
        send_buffer[i] = kvp.second[i];
      }
      mpi_rc = MPI_Send( &send_buffer[0],
                         n_elements,
                         MPI_INT,
                         dst,
                         tag,
                         MPI_COMM_WORLD );
    }
  }
  
  // Complete all of the receives
  mpi_rc = MPI_Waitall( n_recv_reqs, recv_reqs, MPI_STATUSES_IGNORE );

  mpi_rc = MPI_Barrier( MPI_COMM_WORLD );

#ifdef REPORT_PROGRESS
  std::cout << "Rank: " << rank 
            << " completed exchange for communicator: " << current_comm_id 
            << std::endl;
#endif 
  // Now that we have all of the corresponding recv vertex IDs for this 
  // dumpi_to_graph process's sends, we can make the rest of the message edges
  for ( auto kvp : this->channel_to_send_seq ) {
    if ( kvp.first.get_comm() == current_comm_id ) {
      int offset = channel_to_offset.at( kvp.first );
      int n_sends = kvp.second.size();
      for ( int i=0; i<n_sends; ++i ) {
        auto edge = std::make_pair( kvp.second[i], recv_buffer[offset + i] );
        this->message_order_edges.push_back( edge );
      }
    }
  }
#ifdef REPORT_PROGRESS
  std::cout << "Rank: " << rank 
            << " constructed message edges for communicator: " << current_comm_id 
            << std::endl;
#endif 
}



// A function to apply scalar logical timestamps to each vertex. Effectively, 
// we're implementing Lamport's logical clock. 
// Note: Throughout, we use the abbreviation "lts" == "logical timestamp"
// Note: We default to a "ticking policy" that always increments by 1
void EventGraph::apply_scalar_logical_clock()
{
  std::unordered_map<uint8_t, std::string> type_to_name =
  {
    {0, "send"}, {1, "recv"}, {2, "init"}, {3, "finalize"}, {4, "barrier"}
  };

  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  int initial_lts = 0; 
  int tick = 1;
 
  for ( auto kvp : this->rank_to_trace ) {
    auto event_seq = kvp.second->get_event_seq();
    size_t n_vertices = event_seq.size();
    size_t vertex_id_offset = kvp.second->get_vertex_id_offset();
    for ( int i=0; i<n_vertices; ++i ) {
      auto event_type = event_seq[i];
      size_t vertex_id = i + vertex_id_offset;
#ifdef REPORT_PROGRESS_VERBOSE
      std::cout << "Rank: " << rank 
                << " handling trace rank: " << kvp.first
                << " assigning LTS to vertex: " << vertex_id
                << " with event type: " << type_to_name.at( event_type )
                << std::endl;
#endif
      // Case 1: Init events always get the initial logical time stamp
      if ( event_type == 2 ) {
#ifdef PARANOID_INSERTION
        auto vid_search = this->logical_timestamps.find( vertex_id );
        if ( vid_search == this->logical_timestamps.end() ) {
          this->logical_timestamps.insert( { vertex_id, initial_lts } );
        } else {
          std::stringstream oss;
          oss << "Vertex ID: " << vertex_id 
              << " has already been assigned a logical timestamp" 
              << std::endl;
          throw std::runtime_error( oss.str() );
        }
#else
        this->logical_timestamps.insert( { vertex_id, initial_lts } );
#endif
      }
      // Case 2: Vertex represents a send
      // Assign it it's local predecessor's logical timestamp plus the tick, 
      // then send its logical timestamp to its remote successor.
      else if ( event_type == 0 ) {
        // Get local predecessor's lts
        size_t local_pred_vertex_id = vertex_id - 1;
       
        // Check that this send vertex's local predecessor already has a logical 
        // time stamp
        auto search = this->logical_timestamps.find( local_pred_vertex_id );
        if ( search == this->logical_timestamps.end() ) {
          std::ostringstream oss;
          oss << "Rank: " << rank 
              << " send vertex: " << vertex_id
              << " local predecessor: " << local_pred_vertex_id
              << " has no lts" 
              << std::endl;
          throw std::runtime_error( oss.str() ); 
        }
        size_t local_pred_lts = this->logical_timestamps.at( local_pred_vertex_id );

        // Get my lts by incrementing local pred's
        size_t lts = local_pred_lts + tick;
        this->logical_timestamps.insert( { vertex_id, lts } );
        
        // Check that this send vertex has a channel associated with it so that
        // we can propagate its logical time stamp to its remote successor
        // (i.e., a recv vertex owned by some other dumpi_to_graph process)
        auto search2 = this->vertex_id_to_channel.find( vertex_id );
        if ( search2 == this->vertex_id_to_channel.end() ) {
          std::ostringstream oss;
          oss << "Rank: " << rank 
              << " send vertex: " << vertex_id
              << " not mapped to channel" 
              << std::endl;
          throw std::runtime_error( oss.str() ); 
        }
        Channel channel = this->vertex_id_to_channel.at( vertex_id );

        // Propagate the send vertex's logical time stamp to its remote successor
        int dst = channel.get_dst();
        int tag = channel.get_tag();
        int owner = this->config.lookup_owning_rank( dst );
        mpi_rc = MPI_Send( &lts, 1, my_MPI_SIZE_T, dst, tag, MPI_COMM_WORLD );
      }
      // Case 3: Vertex represents a recv
      // Assign it the max of its local and remote predecessors' logical 
      // timestamps plus the tick
      else if ( event_type == 1 ) {
        // Get local predecessor's lts
        size_t local_pred_vertex_id = vertex_id - 1;
        
        // Check that this recv vertex's local predecessor already has a logical
        // time stamp
        auto search = this->logical_timestamps.find( local_pred_vertex_id );
        if ( search == this->logical_timestamps.end() ) {
          std::ostringstream oss;
          oss << "Rank: " << rank 
              << " recv vertex: " << vertex_id
              << " local predecessor has no lts" 
              << std::endl;
          throw std::runtime_error( oss.str() ); 
        }
        size_t local_pred_lts = this->logical_timestamps.at( local_pred_vertex_id );

        // Get remote predecessor's lts
        Channel channel = this->vertex_id_to_channel.at( vertex_id );
        
        int src = channel.get_src();
        int tag = channel.get_tag();
        int owner = this->config.lookup_owning_rank( src );
        size_t remote_pred_lts = 0;
        mpi_rc = MPI_Recv( &remote_pred_lts, 1, my_MPI_SIZE_T, src, tag,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE );

        // Calculate this recv's logical time stamp according to Lamport clock
        // rule
        size_t lts = std::max( local_pred_lts, remote_pred_lts ) + tick;
        this->logical_timestamps.insert( { vertex_id, lts } );

      }
      // Case 4: Everything else. 
      // Currently just vertices representing "barrier" events
      else {
        // Get local predecessor's lts
        size_t local_pred_vertex_id = vertex_id - 1;
        size_t local_pred_lts = this->logical_timestamps.at( local_pred_vertex_id );
        // Get my lts by incrementing local pred's
        size_t lts = local_pred_lts + tick;
        this->logical_timestamps.insert( { vertex_id, lts } );
      }
    } // Loop over vertex sequence for a single trace rank
  } // Loop over all trace ranks managed by this dumpi_to_graph rank
}



void EventGraph::merge()
{
  boost::mpi::communicator world;
  int rank = world.rank();

  // Unique tags for each type of data the root dumpi_to_graph process will 
  // receive from all of the others
  int vertex_ids_tag         = 0;
  int lts_map_tag            = 1;
  int event_type_map_tag     = 2;
  int wall_time_map_tag      = 3;
  int pid_map_tag            = 4;
  int callstack_map_tag      = 5;
  int papi_map_tag           = 433;
  int message_order_edge_tag = 17;
  int program_order_edge_tag = 36;

  // Accumulate all graph data into these
  std::vector<size_t> vertex_ids = this->vertex_ids; 
  std::unordered_map<size_t,size_t> vertex_id_to_lts = this->logical_timestamps;
  std::unordered_map<size_t,uint8_t> vertex_id_to_event_type = this->vertex_id_to_event_type;
  std::unordered_map<size_t,double> vertex_id_to_wall_time = this->vertex_id_to_wall_time;
  std::unordered_map<size_t,int> vertex_id_to_pid = this->vertex_id_to_pid;
  std::unordered_map<size_t,std::string> vertex_id_to_papi = this->vertex_id_to_papi;
  std::unordered_map<size_t,std::string> vertex_id_to_callstack = this->vertex_id_to_callstack;

  std::vector<std::pair<size_t,size_t>> message_order_edges = this->message_order_edges;
  std::vector<std::pair<size_t,size_t>> program_order_edges = this->program_order_edges;

  // Root posts receives for vertices, vertex labels, and edges 
  // All other dumpi_to_graph processes send their partial view of the graph
  if ( rank == 0 ) {
    for ( int i=1; i < world.size(); ++i ) {
      // Receive vertex IDs and vertex labels
      std::vector<size_t> vertex_ids_recv_buffer;
      world.recv( i, vertex_ids_tag, vertex_ids_recv_buffer );
      for ( auto vid : vertex_ids_recv_buffer ) {
        vertex_ids.push_back( vid );
      }

      // Receive vertex label maps
      std::unordered_map<size_t,size_t> lts_map_recv_buffer;
      std::unordered_map<size_t,uint8_t> event_type_map_recv_buffer;
      std::unordered_map<size_t,double> wall_time_map_recv_buffer;
      std::unordered_map<size_t,int> pid_map_recv_buffer;
      std::unordered_map<size_t,std::string> papi_map_recv_buffer;
      std::unordered_map<size_t,std::string> callstack_map_recv_buffer;
      // Receive logical timestamps
      world.recv( i, lts_map_tag, lts_map_recv_buffer );
      for ( auto kvp : lts_map_recv_buffer ) {
        vertex_id_to_lts.insert( kvp );
      }
      // Receive event types
      world.recv( i, event_type_map_tag, event_type_map_recv_buffer );
      for ( auto kvp : event_type_map_recv_buffer ) {
        vertex_id_to_event_type.insert( kvp );
      }
      // Receive wall-time timestamps
      world.recv( i, wall_time_map_tag, wall_time_map_recv_buffer );
      for ( auto kvp : wall_time_map_recv_buffer ) {
        vertex_id_to_wall_time.insert( kvp );
      }
      // Receive process IDs
      world.recv( i, pid_map_tag, pid_map_recv_buffer );
      for ( auto kvp : pid_map_recv_buffer ) {
        vertex_id_to_pid.insert( kvp );
      }

      if ( this->config.has_csmpi() ) {
        // Receive callstacks
        world.recv( i, callstack_map_tag, callstack_map_recv_buffer );
        for ( auto kvp : callstack_map_recv_buffer ) {
          vertex_id_to_callstack.insert( kvp );
        } 
      }

      if ( this->config.get_papi_flag() ) {
        //Receive PAPI counters
        world.recv( i, papi_map_tag, papi_map_recv_buffer );
        for( auto kvp : papi_map_recv_buffer ) {
          vertex_id_to_papi.insert( kvp );
        }
      }

      // A common recv buffer for edges
      std::vector<std::pair<size_t,size_t>> edges_recv_buffer;
      
      // Receive message order edges
      world.recv( i, message_order_edge_tag, edges_recv_buffer );
      for ( auto edge : edges_recv_buffer ) {
        message_order_edges.push_back( edge );
      }

      // Receive program order edges
      world.recv( i, program_order_edge_tag, edges_recv_buffer );
      for ( auto edge : edges_recv_buffer ) {
        program_order_edges.push_back( edge );
      }
    }
    
  }
  else {
    // Send vertex IDs
    world.send( 0, vertex_ids_tag, this->vertex_ids );
    // Send vertex label maps
    world.send( 0, lts_map_tag, this->logical_timestamps );
    world.send( 0, event_type_map_tag, this->vertex_id_to_event_type );
    world.send( 0, wall_time_map_tag, this->vertex_id_to_wall_time );
    world.send( 0, pid_map_tag, this->vertex_id_to_pid );
    if ( this->config.has_csmpi() ) {
      world.send( 0, callstack_map_tag, this->vertex_id_to_callstack );
    }
    fprintf(stderr, "%d flag\n", this->config.get_papi_flag());
    if ( this->config.get_papi_flag() ) {
      world.send(0, papi_map_tag, this->vertex_id_to_papi );
    }
    // Send edges
    world.send( 0, message_order_edge_tag, this->message_order_edges );
    world.send( 0, program_order_edge_tag, this->program_order_edges );
  }

  // Root constructs the igraph representation
  if ( rank == 0 ) {
    std::unordered_map<uint8_t, std::string> type_to_name =
    {
      {0, "send"}, {1, "recv"}, {2, "init"}, {3, "finalize"}, {4, "barrier"}
    };
    int igraph_rc;
    // Turn on the igraph attribute handler 
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    
    // Initialize the graph
    size_t n_vertices = vertex_ids.size();
    igraph_t graph;

    // boolean parameter makes the graph directed
    igraph_rc = igraph_empty( &graph, n_vertices, true ); 

#ifdef REPORT_PROGRESS
    std::cout << "Base igraph object constructed" << std::endl;
#endif

    // Add vertex attributes
    const char event_type_attr_name[16] = "event_type";
    const char lts_attr_name[16]        = "logical_time";
    const char wtime_attr_name[16]      = "wall_time";
    const char pid_attr_name[16]        = "process_id";
    const char callstack_attr_name[16]  = "callstack";
    const char papi_attr_name[16]       = "PAPI_ctrs";

    for ( auto vid : vertex_ids ) {

      // Get vertex labels available directly from DUMPI trace
      uint8_t event_type = vertex_id_to_event_type[vid];
      const std::string event_type_str = type_to_name[event_type];
      size_t lts = vertex_id_to_lts[vid];
      double wtime = vertex_id_to_wall_time[vid];
      int pid = vertex_id_to_pid[vid];

      // Lookup callstack label if available from CSMPI trace
      std::string callstack;
      if ( this->config.has_csmpi() ) {
        callstack = vertex_id_to_callstack[vid];
      }
      std::string papi_counters;
      if ( this->config.get_papi_flag() ) {
        papi_counters = vertex_id_to_papi[vid];
      }

      // Set event type vertex attribute
      igraph_rc = igraph_cattribute_VAS_set( &graph, event_type_attr_name, vid, event_type_str.c_str() );
      // Set logical timestamp vertex attribute
      igraph_rc = igraph_cattribute_VAN_set( &graph, lts_attr_name, vid, lts );
      // Set wall-time timestamp vertex attribute
      igraph_rc = igraph_cattribute_VAN_set( &graph, wtime_attr_name, vid, wtime );
      // Set process ID vertex attribute
      igraph_rc = igraph_cattribute_VAN_set( &graph, pid_attr_name, vid, pid );
      if ( this->config.has_csmpi() ) {
        // Set callstack vertex attribute
        igraph_rc = igraph_cattribute_VAS_set( &graph, callstack_attr_name, vid, callstack.c_str() );
      }
      if ( this->config.get_papi_flag() ) {
        //Set PAPI counters vertex attribute
        igraph_rc = igraph_cattribute_VAS_set( &graph, papi_attr_name, vid, papi_counters.c_str() );
      }
    }
#ifdef REPORT_PROGRESS
    std::cout << "Vertex attributes added" << std::endl;
#endif
    // Add edges
    igraph_vector_t edges;
    size_t n_edges = 2 * ( program_order_edges.size() + message_order_edges.size() );
    igraph_vector_init( &edges, n_edges ); 
    size_t edge_idx = 0;
    for ( auto edge : program_order_edges ) {
      VECTOR(edges)[ edge_idx ] = edge.first;
      edge_idx++;
      VECTOR(edges)[ edge_idx ] = edge.second;
      edge_idx++;
    }
    for ( auto edge : message_order_edges ) {
      VECTOR(edges)[ edge_idx ] = edge.first;
      edge_idx++;
      VECTOR(edges)[ edge_idx ] = edge.second;
      edge_idx++;
    }
    igraph_rc = igraph_add_edges( &graph, &edges, 0 );
 
    // Assign
    this->_graph = graph;

  }
}

void EventGraph::write() const
{
  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if ( rank == 0 ) {
    int igraph_rc;

    std::string trace_dir = this->config.get_trace_dirs()[0];
    std::stringstream ss;
    ss << trace_dir << "/event_graph.graphml"; 
    std::string output_path = ss.str();

    FILE* outfile;
    outfile = fopen( output_path.c_str(), "w" );
    igraph_rc = igraph_write_graph_graphml( &(this->_graph), outfile, false );
    fclose( outfile );
#ifdef REPORT_PROGRESS_MINIMAL
    std::cout << "Wrote event graph to: " << output_path
              << std::endl;
#endif
  }
}



////////////////////////////////////////////////////////////////////////////////
// Convenience functions for printing the vertices and edges of an EventGraph //
////////////////////////////////////////////////////////////////////////////////

void EventGraph::report_program_order_edges() const
{ 
  std::unordered_map<uint8_t, std::string> type_to_name =
  {
    {0, "send"}, {1, "recv"}, {2, "init"}, {3, "finalize"}, {4, "barrier"}
  };
  for ( auto edge : this->program_order_edges ) {
    auto src_vertex_id = edge.first;
    auto dst_vertex_id = edge.second;
    auto src_vertex_type = this->vertex_id_to_event_type.at( src_vertex_id );
    auto dst_vertex_type = this->vertex_id_to_event_type.at( dst_vertex_id );
    auto src_vertex_lts = this->logical_timestamps.at( src_vertex_id );
    auto dst_vertex_lts = this->logical_timestamps.at( dst_vertex_id );
    std::cout << "Program Order Edge: "
              << "ID: " << src_vertex_id 
              << ", Type: " << type_to_name.at( src_vertex_type )
              << ", LTS: " << src_vertex_lts
              << " --> "
              << "ID: " << dst_vertex_id
              << ", Type: " << type_to_name.at( dst_vertex_type )
              << ", LTS: " << dst_vertex_lts
              << std::endl;
  }
}

void EventGraph::report_message_order_edges() const
{
  for ( auto edge : this->message_order_edges ) {
    size_t src_vertex_id = edge.first;
    size_t dst_vertex_id = edge.second;
    size_t src_vertex_lts = this->logical_timestamps.at( src_vertex_id );
    size_t dst_vertex_lts;
    auto search = this->logical_timestamps.find( dst_vertex_id ); 
    if ( search != this->logical_timestamps.end() ) {
      dst_vertex_lts = this->logical_timestamps.at( dst_vertex_id );
    } else {
      dst_vertex_lts = 0;
    }
    auto channel = this->vertex_id_to_channel.at( src_vertex_id );
    std::cout << "Message Order Edge: "
              << "ID: " << src_vertex_id 
              << ", Type: send" 
              << ", LTS: " << src_vertex_lts
              << " --> "
              << "ID: " << dst_vertex_id
              << ", Type: recv"
              << ", LTS: " << dst_vertex_lts
              << ", in Channel: " << channel
              << std::endl;
  }
}



