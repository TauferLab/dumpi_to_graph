#include "EventGraph.hpp"

#include <mpi.h>

#include <utility> // make_pair
#include <algorithm> // max
#include <cstdint>
#include <cinttypes>
#include <unordered_map>
#include <vector>
#include <string>

#include "boost/serialization/unordered_map.hpp" 
#include "boost/mpi.hpp"

#include "Debug.hpp"

// Constructor must do the following:
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
                        const std::unordered_map<int,Trace*> rank_to_trace )
{
  // Set members
  this->config = config;
  this->rank_to_trace = rank_to_trace;

  // First determine the offsets to apply to all of my vertex IDs
  std::unordered_map<int,size_t> trace_rank_to_vertex_count;
  for ( auto kvp : this->rank_to_trace ) {
    size_t vertex_count = kvp.second->get_final_vertex_id() + 1;
    trace_rank_to_vertex_count.insert( { kvp.first, vertex_count } );
  }
  
  boost::mpi::communicator comm_world;

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
  
//#ifdef SANITY_CHECK
//  comm_world.barrier();
//  for ( int i=0; i<comm_world.size(); i++) {
//    if ( comm_world.rank() == i ) {
//      for ( auto kvp : this->rank_to_trace ) {
//        std::cout << "Rank: " << comm_world.rank() 
//                  << " handling trace rank: " << kvp.first << std::endl;
//        kvp.second->report_event_seq();
//        kvp.second->report_channel_to_send_seq();
//        kvp.second->report_channel_to_recv_seq();
//      }
//      std::cout << std::endl;
//    }
//    comm_world.barrier();
//  }
//#endif

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
    // Merge in event type data
    auto event_seq = kvp.second->get_event_seq();
    size_t n_events = event_seq.size();
    size_t initial_vertex_id = kvp.second->get_initial_vertex_id();
    for ( int i=0; i<n_events; ++i ) {
      uint8_t event_type = event_seq[i];
      size_t vertex_id = initial_vertex_id + i;
      this->vertex_id_to_event_type.insert( { vertex_id, event_type } );
    }
  }

//#ifdef SANITY_CHECK
//  comm_world.barrier();
//  std::unordered_map<uint8_t, std::string> type_to_name =
//  {
//    {0, "send"}, {1, "recv"}, {2, "init"}, {3, "finalize"}, {4, "barrier"}
//  };
//  for ( int i=0; i<comm_world.size(); i++) {
//    if ( comm_world.rank() == i ) {
//      std::cout << "Rank: " << comm_world.rank() << std::endl;
//      for ( auto kvp : this->vertex_id_to_event_type ) {
//        std::cout << "Vertex ID: " << kvp.first 
//                  << ", Event Type: " << type_to_name.at( kvp.second ) 
//                  << std::endl;
//      }
//      for ( auto kvp : this->channel_to_send_seq ) {
//        std::cout << "Channel: " << kvp.first << ", Send Vertex IDs: ";
//        for ( auto send : kvp.second ) {
//          std::cout << " " << send;
//        }
//        std::cout << std::endl;
//      }
//      for ( auto kvp : this->channel_to_recv_seq ) {
//        std::cout << "Channel: " << kvp.first << ", Recv Vertex IDs: ";
//        for ( auto recv : kvp.second ) {
//          std::cout << " " << recv;
//        }
//        std::cout << std::endl;
//      }
//      std::cout << std::endl;
//    }
//    comm_world.barrier();
//  }
//#endif

  comm_world.barrier();

  // Construct edges
  this->make_program_order_edges();
  this->make_message_order_edges();
  this->make_collective_edges();

//#ifdef SANITY_CHECK
//  for ( int i=0; i<comm_world.size(); i++) {
//    if ( comm_world.rank() == i ) {
//      comm_world.barrier(); 
//      std::cout << std::endl;
//      std::cout << "Message edges for rank: " << comm_world.rank() << std::endl;
//      this->report_program_order_edges();
//      this->report_message_order_edges();
//      std::cout << std::endl;
//      comm_world.barrier(); 
//    }
//  }
//#endif
}



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

// Helper function for building the message edges between sends and receives
// held on the same dumpi_to_graph process
void EventGraph::exchange_local_message_matching_data()
{
  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  

  // We need to iterate over this local copy so that we don't invalidate when erasing
  auto channel_to_send_seq_copy = this->channel_to_send_seq;

  for ( auto kvp : channel_to_send_seq_copy ) {
    
    //std::cout << "Rank: " << rank << " handling channel: " << kvp.first << std::endl;

    // The rank of the receiving process in the traced application
    int dst = kvp.first.get_dst();
    // The dumpi_to_graph process managing that rank
    int owning_rank = this->config.lookup_owning_rank( dst );

    //std::cout << "Rank that owns destination recv seq is: " << owning_rank << std::endl;

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
      //std::cout << "Size of channel map before erase: " << this->channel_to_send_seq.size() << std::endl;
      this->channel_to_send_seq.erase( kvp.first );
      this->channel_to_recv_seq.erase( kvp.first );
      //std::cout << "Size of channel map after erase: " << this->channel_to_send_seq.size() << std::endl;
    }
  }
}

// Helper function for building the message edges between sends and receives 
// held on distinct dumpi_to_graph processes.
// Assumes that EventGraph::exchange_local_message_matching_data() has already
// executed (i.e., that it has reduced the channel maps to include only channels
// between trace processes that are held on different dumpi_to_graph processes)
void EventGraph::exchange_remote_message_matching_data()
{
  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n_recv_reqs = this->channel_to_send_seq.size();
  MPI_Request recv_reqs[ n_recv_reqs ];
  int recv_req_idx = 0;

  // We will receive all remotely held recv sequences into a single contiguous
  // array. We determine the size of that array and a mapping between the 
  // relevant channels and corresponding offsets into the array here
  std::unordered_map<Channel, int, ChannelHash> channel_to_offset;
  int n_recv_vertex_ids = 0;
  for ( auto kvp : this->channel_to_send_seq ) {
    channel_to_offset.insert( { kvp.first, n_recv_vertex_ids } );
    n_recv_vertex_ids += kvp.second.size();
  }

  // The shared receive buffer that will contain all of the recv vertex IDs from
  // remote dumpi_to_graph processes
  int recv_buffer[ n_recv_vertex_ids ]; 

  // Loop over our send sequences and post a receive for each matching recv
  // sequence at the appropriate offset
  for ( auto kvp : this->channel_to_send_seq ) {
    // The trace rank of the receiver
    int src = kvp.first.get_dst();
    // The dumpi_to_graph process we need to get the matching recv sequence from
    int owner = this->config.lookup_owning_rank( src );
    // The tag of the message
    int tag = kvp.first.get_tag();
    // Offset into the recv buffer
    int recv_buffer_offset = channel_to_offset.at( kvp.first );
    // Number of recv vertex IDs we need to receive
    int n_elements = kvp.second.size();
#ifdef SANITY_CHECK
    // This should never happen if the local message matching data exchange has
    // executed (properly), but no amount of paranoia is excessive in the realm
    // of MPI applications
    assert( owner != rank );
#endif
    // Actually post the receive
    // FIXME: Right now we do all of these exchanges in MPI_COMM_WORLD, but we 
    // should probably mirror the communicators of the traced application. 
    // Not a high priority right now since our traced applications
    // of interest (e.g., Enzo) do everything in MPI_COMM_WORLD themselves.
    mpi_rc = MPI_Irecv( &recv_buffer[ recv_buffer_offset ],
                        n_elements,
                        MPI_INT, 
                        owner,
                        tag,
                        MPI_COMM_WORLD,
                        &recv_reqs[ recv_req_idx ] );
    // Update request index
    recv_req_idx++;
  }
  
  // Loop over our recv sequences and send each to the dumpi_to_graph process
  // with the matching send sequence
  for ( auto kvp : this->channel_to_recv_seq ) {
    // The trace rank of the sender
    int dst = kvp.first.get_src(); 
    // The dumpi_to_graph process we are sending this recv sequence to
    int owner = this->config.lookup_owning_rank( dst );
    // The tag of the message
    int tag = kvp.first.get_tag();
    // Number of recv vertex IDs we are sending
    int n_elements = kvp.second.size();
#ifdef SANITY_CHECK
    assert( owner != rank ); 
    assert( contains_no_invalid_vertex_ids( kvp.second ) );
#endif
    // FIXME: this copy is probably not necessary, but it's working and doesn't
    // seem to be a major source of slowdown 
    int send_buffer[ n_elements ];
    for ( int i=0; i<n_elements; ++i ) {
      send_buffer[i] = kvp.second[i];
    }
    // FIXME: We do a blocking send here because the non-blocking version was
    // delivering corrupted data to the receiver. Pretty sure it has something
    // to do with how MPI_Isend requires you to be very careful about not 
    // overwriting the send buffer (which I didn't think I was doing but...) 
    // Anyway, we use a blocking standard mode send now (which seems to be 
    // working fine on MVAPICH 2.3 at least). There's probably an edge case 
    // where this can deadlock though, so we should eventually implement this
    // as a buffered mode send. 
    // FIXME: Same MPI_COMM_WORLD issue as with the receives above
    mpi_rc = MPI_Send( &send_buffer[0],
                       n_elements,
                       MPI_INT,
                       owner,
                       tag,
                       MPI_COMM_WORLD );
  }

  // Complete all of the receives
  mpi_rc = MPI_Waitall( n_recv_reqs, recv_reqs, MPI_STATUSES_IGNORE );

#ifdef SANITY_CHECK
  assert( validate_remote_recv_seqs( recv_buffer, 
                                     n_recv_vertex_ids,
                                     this->channel_to_send_seq, 
                                     channel_to_offset ) );
#endif
  
  // Now that we have all of the corresponding recv vertex IDs for this 
  // dumpi_to_graph process's sends, we can make the rest of the message edges
  for ( auto kvp : this->channel_to_send_seq ) {
    int offset = channel_to_offset.at( kvp.first );
    int n_sends = kvp.second.size();
    for ( int i=0; i<n_sends; ++i ) {
      auto edge = std::make_pair( kvp.second[i], recv_buffer[offset + i] );
      this->message_order_edges.push_back( edge );
    }
  }
}

void EventGraph::make_collective_edges()
{
  // FIXME: an implementation would be nice 
}


// A function to apply scalar logical timestamps to each vertex. Effectively, 
// we're implementing Lamport's logical clock. 
// Note: Throughout, we use the abbreviation "lts" == "logical timestamp"
// Note: We default to a "ticking policy" that always increments by 1
void EventGraph::apply_scalar_logical_clock()
{
  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  int initial_lts = 0; 
  int tick = 1;
  
  std::cout << "Rank: " << rank << " applying logical clock" << std::endl;
  
  for ( auto kvp : this->rank_to_trace ) {
    auto event_seq = kvp.second->get_event_seq();
    size_t n_vertices = event_seq.size();
    for ( int i=0; i<n_vertices; ++i ) {
      auto event_type = event_seq[i];
      if ( event_type == 2 ) {

      }
    }
  }


  //int n_vertices = this->vertex_ids.size();
  //for ( int i=0; i<n_vertices; ++i) {
  //  int vertex_type = this->vertex_types[ i ];
  //  // Case 1: Vertex represents MPI_Init or MPI_Init_thread. 
  //  // Assign it the initial logical timestamp value
  //  if ( vertex_type == 0 ) {
  //    this->logical_timestamps[ i ] = initial_lts;
  //  }
  //  // Case 2: Vertex represents a send
  //  // Assign it it's local predecessor's logical timestamp plus the tick, then
  //  // send its logical timestamp to its remote successor
  //  else if ( vertex_type == 1 ) {
  //    int local_pred_lts = this->logical_timestamps[ i-1 ];
  //    int lts = local_pred_lts + tick;
  //    this->logical_timestamps[ i ] = lts;
  //    Channel chan = this->vid_to_chan.at( this->vertex_ids[ i ] );
  //    int dst = chan.get_dst();
  //    int tag = chan.get_tag();
  //    // FIXME: Pretty sure this is depending on eager delivery of small 
  //    // messages (i.e., single int) to not deadlock
  //    mpi_rc = MPI_Send( &lts, 1, MPI_INT, dst, tag, MPI_COMM_WORLD );
  //  }
  //  // Case 3: Vertex represents a recv
  //  // Assign it the max of its local and remote predecessors' logical 
  //  // timestamps plus the tick
  //  else if ( vertex_type == 2 ) {
  //    int local_pred_lts = this->logical_timestamps[ i-1 ];
  //    Channel chan = this->vid_to_chan.at( this->vertex_ids[ i ] );
  //    int src = chan.get_src();
  //    int tag = chan.get_tag();
  //    int remote_pred_lts;
  //    mpi_rc = MPI_Recv( &remote_pred_lts, 1, MPI_INT, src, tag, 
  //                       MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  //    int lts = std::max( local_pred_lts, remote_pred_lts ) + tick;
  //    this->logical_timestamps[ i ] = lts;
  //  }
  //}
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
    std::cout << "Program Order Edge: "
              << "ID: " << src_vertex_id 
              << ", Type: " << type_to_name.at( src_vertex_type )
              << " --> "
              << "ID: " << dst_vertex_id
              << ", Type: " << type_to_name.at( dst_vertex_type )
              << std::endl;
  }
}

void EventGraph::report_message_order_edges() const
{
  for ( auto edge : this->message_order_edges ) {
    auto src_vertex_id = edge.first;
    auto dst_vertex_id = edge.second;
    auto channel = this->vertex_id_to_channel.at( src_vertex_id );
    std::cout << "Message Order Edge: "
              << "ID: " << src_vertex_id 
              << ", Type: send" 
              << " --> "
              << "ID: " << dst_vertex_id
              << ", Type: recv" 
              << ", in Channel: " << channel
              << std::endl;
  }
}



