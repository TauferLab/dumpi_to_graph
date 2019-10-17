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
#include <fstream>
#include <iostream>

// Boost
#include "boost/serialization/unordered_map.hpp" 
#include "boost/mpi.hpp"

// Igraph
#include "igraph/igraph.h"

// Internal
#include "Logging.hpp"
#include "Debug.hpp"

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
    // Merge in event type, pid, and wall-time data
    auto event_seq = kvp.second->get_event_seq();
    auto wall_time_seq = kvp.second->get_wall_time_seq(); 
    size_t n_events = event_seq.size();
    size_t initial_vertex_id = kvp.second->get_initial_vertex_id();
    for ( int i=0; i<n_events; ++i ) {
      uint8_t event_type = event_seq[i];
      double wall_time = wall_time_seq[i];
      size_t vertex_id = initial_vertex_id + i;
      int pid = kvp.first; 
      this->vertex_ids.push_back( vertex_id );
      this->vertex_id_to_event_type.insert( { vertex_id, event_type } );
      this->vertex_id_to_wall_time.insert( { vertex_id, wall_time } );
      this->vertex_id_to_pid.insert( { vertex_id, pid } );
    }
  }

  comm_world.barrier();

  // Construct edges
  this->make_program_order_edges();
  this->make_message_order_edges();
  this->make_collective_edges();

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
  std::unordered_map<uint8_t, std::string> type_to_name =
  {
    {0, "send"}, {1, "recv"}, {2, "init"}, {3, "finalize"}, {4, "barrier"}
  };

  int mpi_rc, rank;
  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  int initial_lts = 0; 
  int tick = 1;
 
#ifdef REPORT_PROGRESS
  std::cout << "Rank: " << rank << " applying logical clock" << std::endl;
#endif
  
  for ( auto kvp : this->rank_to_trace ) {
    auto event_seq = kvp.second->get_event_seq();
    size_t n_vertices = event_seq.size();
    size_t vertex_id_offset = kvp.second->get_vertex_id_offset();
#ifdef REPORT_PROGRESS
    std::cout << "Rank: " << rank
              << " handling trace rank: " << kvp.first
              << " with event sequence of length: " << n_vertices
              << " and vertex ID offset: " << vertex_id_offset
              << std::endl;
#endif
    for ( int i=0; i<n_vertices; ++i ) {
      auto event_type = event_seq[i];
      size_t vertex_id = i + vertex_id_offset;

      //std::cout << "Rank: " << rank 
      //          << " handling trace rank: " << kvp.first
      //          << " assigning LTS to vertex: " << vertex_id
      //          << " with event type: " << type_to_name.at( event_type )
      //          << std::endl;

      // Case 1: Init events always get the initial logical time stamp
      if ( event_type == 2 ) {
#ifdef PARANOID_INSERTION
        auto vid_search = this->logical_timestamps.find( vertex_id );
        if ( vid_search == this->logical_timestamps.end() ) {
          this->logical_timestamps.insert( { vertex_id, initial_lts } );
        } else {
          std::stringstream ss;
          ss << "Vertex ID: " << vertex_id 
             << " has already been assigned a logical timestamp" << std::endl;
          throw std::runtime_error( ss.str() );
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
        
        auto search = this->logical_timestamps.find( local_pred_vertex_id );
        if ( search == this->logical_timestamps.end() ) {
          std::cout << "Rank: " << rank << " send vertex: " << vertex_id
                    << " local predecessor: " << local_pred_vertex_id
                    << " has no lts" << std::endl;
          exit(1);
        }
        size_t local_pred_lts = this->logical_timestamps.at( local_pred_vertex_id );
        // Get my lts by incrementing local pred's
        size_t lts = local_pred_lts + tick;
        this->logical_timestamps.insert( { vertex_id, lts } );
        
        auto search2 = this->vertex_id_to_channel.find( vertex_id );
        if ( search2 == this->vertex_id_to_channel.end() ) {
          std::cout << "Rank: " << rank << " send vertex: " << vertex_id
                    << " not mapped to channel" << std::endl;
          exit(1);
        }
        
        Channel channel = this->vertex_id_to_channel.at( vertex_id );
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
        
        auto search = this->logical_timestamps.find( local_pred_vertex_id );
        if ( search == this->logical_timestamps.end() ) {
          std::cout << "Rank: " << rank << " recv vertex: " << vertex_id
                    << " local predecessor has no lts" << std::endl;
          exit(0);
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
        size_t lts = std::max( local_pred_lts, remote_pred_lts ) + tick;
        this->logical_timestamps.insert( { vertex_id, lts } );

      }
      // Case 4: Everything else. 
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



void EventGraph::merge_and_write()
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
  int message_order_edge_tag = 17;
  int program_order_edge_tag = 36;

  // Accumulate all graph data into these
  std::vector<size_t> vertex_ids = this->vertex_ids; 
  std::unordered_map<size_t,size_t> vertex_id_to_lts = this->logical_timestamps;
  std::unordered_map<size_t,uint8_t> vertex_id_to_event_type = this->vertex_id_to_event_type;
  std::unordered_map<size_t,double> vertex_id_to_wall_time = this->vertex_id_to_wall_time;
  std::unordered_map<size_t,int> vertex_id_to_pid = this->vertex_id_to_pid;
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
    for ( auto vid : vertex_ids ) {
      uint8_t event_type = vertex_id_to_event_type[vid];
      const std::string event_type_str = type_to_name[event_type];
      size_t lts = vertex_id_to_lts[vid];
      double wtime = vertex_id_to_wall_time[vid];
      int pid = vertex_id_to_pid[vid];
      // Set event type vertex attribute
      igraph_rc = igraph_cattribute_VAS_set( &graph, event_type_attr_name, vid, event_type_str.c_str() );
      // Set logical timestamp vertex attribute
      igraph_rc = igraph_cattribute_VAN_set( &graph, lts_attr_name, vid, lts );
      // Set wall-time timestamp vertex attribute
      igraph_rc = igraph_cattribute_VAN_set( &graph, wtime_attr_name, vid, wtime );
      // Set process ID vertex attribute
      igraph_rc = igraph_cattribute_VAN_set( &graph, pid_attr_name, vid, pid );
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
  
    std::string trace_dir = this->config.get_trace_dirs()[0];
    std::stringstream ss;
    ss << trace_dir << "/event_graph.graphml"; 
    std::string output_path = ss.str();

    FILE* outfile;
    outfile = fopen( output_path.c_str(), "w" );
#ifdef REPORT_PROGRESS
    std::cout << "Writing event graph..." << std::endl;
#endif
    igraph_rc = igraph_write_graph_graphml( &graph, outfile, false );
    fclose( outfile );
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



