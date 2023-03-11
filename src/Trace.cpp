#include "Trace.hpp"

#include <cinttypes>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>

// DUMPI
#include "dumpi/common/argtypes.h" 
#include "dumpi/common/constants.h" // DUMPI_ANY_SOURCE, DUMPI_ANY_TAG, DUMPI_COMM_WORLD
#include "dumpi/common/types.h" // dumpi_clock 

// Internal
#include "Logging.hpp"
#include "Channel.hpp"
#include "Request.hpp"
#include "Debug.hpp"
#include "CommunicatorManager.hpp"
#include "Utilities.hpp"

Trace::Trace( Configuration config, 
              std::string trace_dir,
              int trace_rank,
              int dumpi_to_graph_rank )
{
  // Set members
  this->config = config;
  this->trace_dir = trace_dir;
  this->trace_rank = trace_rank;
  this->dumpi_to_graph_rank = dumpi_to_graph_rank;
  
  this->pluto_trace.open(this->trace_dir + "/pluto_out" + std::to_string(this->trace_rank)+".txt");
  std::cout << "Using pluto" << this->trace_dir + "/pluto_out" + std::to_string(this->trace_rank) << ".txt" << std::endl;
  this->pluto_trace.ignore(256, '\n');
  // Get # trace ranks 
  int mpi_rc;
  int reduce_recv_buffer;
  int reduce_send_buffer = trace_rank;
  int count = 1;
  int root = 0;
  mpi_rc = MPI_Reduce( &reduce_send_buffer,
                       &reduce_recv_buffer,
                       count,
                       MPI_INT,
                       MPI_MAX,
                       root,
                       MPI_COMM_WORLD );
  int bcast_buffer;
  if ( dumpi_to_graph_rank == 0 ) {
    bcast_buffer = reduce_recv_buffer;
  }
  mpi_rc = MPI_Bcast( &bcast_buffer,
                      count,
                      MPI_INT,
                      root,
                      MPI_COMM_WORLD );
  n_trace_ranks = bcast_buffer + 1;
  
  // Set size of global communicator
  this->comm_manager = CommunicatorManager( n_trace_ranks );
}

// Keep track of which MPI function call (of all the various types encountered 
// in the trace) that we are currently on.
void Trace::update_call_idx( std::string mpi_fn )
{
  auto search = _mpi_fn_to_idx.find( mpi_fn );
  if ( search == _mpi_fn_to_idx.end() ) {
    _mpi_fn_to_idx.insert( { mpi_fn, 0 } );
  } else {
    _mpi_fn_to_idx[ mpi_fn ]++;
  }
}

void Trace::associate_event_with_call( std::string mpi_fn, size_t event_vertex_id )
{
  // Get current call index for this MPI function
  auto call_idx = _mpi_fn_to_idx.at( mpi_fn );
  auto pair = std::make_pair( mpi_fn, call_idx );
  mpi_fn_seq.push_back( pair );
}

CommunicatorManager& Trace::get_comm_manager()
{
  return this->comm_manager;
} 

void Trace::register_comm_split( int parent_comm_id,
                                 int new_comm_id, 
                                 int color, 
                                 int key )
{
  auto global_rank = this->get_trace_rank();
  // Update communicator-hierarchy
  this->comm_manager.update_comm_to_parent( new_comm_id, parent_comm_id );
  // Associate the global rank with its group (color) and rank (key) within
  // the new communicator
  this->comm_manager.associate_rank_with_color( new_comm_id, global_rank, color );
  this->comm_manager.associate_rank_with_key( new_comm_id, global_rank, key );
}

void Trace::register_barrier( size_t event_vertex_id ) 
{
  // Update the sequence of event types
  this->event_seq.push_back(4);
}

void Trace::register_initial_dumpi_timestamp( const dumpi_time& wall_time )
{
  int32_t wall_time_start_sec = wall_time.start.sec;
  int32_t wall_time_stop_sec  = wall_time.stop.sec;
  int32_t wall_time_start_ns = wall_time.start.nsec;
  int32_t wall_time_stop_ns  = wall_time.stop.nsec;

  double start = (double)wall_time_start_sec + (double)wall_time_start_ns / 10e9;
  double stop = (double)wall_time_stop_sec + (double)wall_time_stop_ns / 10e9;
  double midpoint = (start + stop) / 2;
  
  this->initial_timestamp = midpoint;   
  this->wall_time_seq.push_back( 0.0 );
}

void Trace::register_dumpi_timestamp( const dumpi_time& wall_time )
{
  int32_t wall_time_start_sec = wall_time.start.sec;
  int32_t wall_time_stop_sec  = wall_time.stop.sec;
  int32_t wall_time_start_ns = wall_time.start.nsec;
  int32_t wall_time_stop_ns  = wall_time.stop.nsec;

  double start = (double)wall_time_start_sec + (double)wall_time_start_ns / 10e9;
  double stop = (double)wall_time_stop_sec + (double)wall_time_stop_ns / 10e9;
  double midpoint = (start + stop) / 2;
 
  double timestamp = midpoint - this->initial_timestamp;
  this->wall_time_seq.push_back( timestamp );
}

void Trace::register_papi_struct(const dumpi_perfinfo& counters)
{
  this->counter_sets.push_back(stringify_perfinfo(counters));
}

// Helper for updating channel_to_recv_seq and vertex_id_to_channel
void Trace::register_recv( const Channel& channel, size_t recv_vertex_id )
{
  // std::cout << "JACK_ register_recv" << std::endl;
  // First update the sequence of event types
  this->event_seq.push_back(1);
  // First update the mapping from channels to sequences of recv vertex IDs
  auto chan_search = this->channel_to_recv_seq.find( channel );
  if ( chan_search != this->channel_to_recv_seq.end() ) {
    chan_search->second.push_back( recv_vertex_id );
  } else {
    std::vector<size_t> vertex_id_seq = { recv_vertex_id };
    this->channel_to_recv_seq.insert( { channel, vertex_id_seq } );
  }
  // Next update the mapping from vertex IDs to channels
#ifdef PARANOID_INSERTION
  auto vid_search = this->vertex_id_to_channel.find( recv_vertex_id );
  // Vertex has not been assigned a channel yet, all good
  if ( vid_search == this->vertex_id_to_channel.end() ) {
    this->vertex_id_to_channel.insert( { recv_vertex_id, channel } );
  }
  // Vertex already has an associated channel... we messed up somewhere
  else {
    std::stringstream ss;
    ss << "Trying to assign channel: " << channel
       << " to recv vertex ID: " << recv_vertex_id
       << " but that ID is already assigned channel: " << vid_search->second
       << std::endl;
    throw std::runtime_error( ss.str() );
  }
#else
  this->vertex_id_to_channel.insert( { recv_vertex_id, channel } );
#endif
}

// Helper for updating channel_to_send_seq and vertex_id_to_channel
void Trace::register_send( const Channel& channel, size_t send_vertex_id )
{
  // std::cout << "JACK_ register_send con: " << channel << std::endl;
  // First update the sequence of event types
  this->event_seq.push_back(0);
  // Next update the mapping from channels to sequences of send vertex IDs
  auto search = this->channel_to_send_seq.find( channel );
  if ( search != this->channel_to_send_seq.end() ) {
    search->second.push_back( send_vertex_id );
  } else {
    std::vector<size_t> vertex_id_seq = { send_vertex_id };
    this->channel_to_send_seq.insert( { channel, vertex_id_seq } );
  }
  // Next update the mapping from vertex IDs to channels
#ifdef PARANOID_INSERTION
  auto vid_search = this->vertex_id_to_channel.find( send_vertex_id );
  // Vertex has not been assigned a channel yet, all good
  if ( vid_search == this->vertex_id_to_channel.end() ) {
    this->vertex_id_to_channel.insert( { send_vertex_id, channel } );
  }
  // Vertex already has an associated channel... we messed up somewhere
  else {
    std::stringstream ss;
    ss << "Trying to assign channel: " << channel
       << " to send vertex ID: " << send_vertex_id
       << " but that ID is already assigned channel: " << vid_search->second
       << std::endl;
    throw std::runtime_error( ss.str() );
  }
#else
  this->vertex_id_to_channel.insert( { send_vertex_id, channel } );
#endif
}

// Helper for registering init event
void Trace::register_init(std::string init_fn) 
{
  size_t event_vertex_id = this->get_next_vertex_id();
  this->event_seq.push_back(2);
  this->associate_event_with_call( init_fn, event_vertex_id );
}

// Helper for registering finalize event
void Trace::register_finalize()
{
  size_t event_vertex_id = this->get_next_vertex_id();
  this->event_seq.push_back(3);
  this->final_vertex_id = event_vertex_id;
  this->associate_event_with_call( "MPI_Finalize", event_vertex_id );
}





// Helper for updating id_to_request
void Trace::register_request( long request_id, Request& request )
{ 
  auto search = this->id_to_request.find( request_id );
  // std::cout << "JACK::Adding request with ID: " << request_id << std::endl;
  
  // Case 1: Request not already tracked, insert
  if ( search == this->id_to_request.end() ) {
    this->id_to_request.insert( { request_id, request } );
  } 
  // Case 2: Request already tracked. Error. 
  else {
    // TODO: JACK_ add validation to compare the type of the request
    auto prev_request = search->second; 
    
    if(prev_request.get_type() == request.get_type()){ // request is the same
      // std::cout << "type of request: " << prev_request.get_type() << std::endl; // JACK_
      std::stringstream ss;
      ss << "dumpi_to_graph rank: " << this->get_dumpi_to_graph_rank()
        << " trying to map request ID: " << request_id 
        << " to request: " << request
        << " but request ID is already mapped to request: " << prev_request 
        << std::endl;
      throw std::runtime_error( ss.str() );
    }
    else{
      // this->id_to_request.insert( { request_id, request } );
      // JACK_
      request.renew_id();
      std::cout << "JACK::Adding request with ID: " << request.get_id() << " already existed ID: " << request_id << " rank: " << this->get_trace_rank() << std::endl;

      auto inserted = this->id_to_request.insert( { request.get_id(), request } );
      // this->report_id_to_request();
      // std::stringstream ss;
      // std::cout << "JACK_ Inserted: " << inserted.first->second
      //           << ";; Prev: " << prev_request
      //           << ";; actual: " << request << std::endl;
      // throw std::runtime_error( ss.str() );

    }
  }
}




////////////////////////////////////////////////////////////////////////////////

Channel Trace::determine_channel_of_recv( const dumpi_recv recv )
{
  // Destination is just this trace rank and comm is guaranteed to be in recv
  Channel partial_channel( this->get_trace_rank(), recv.comm );
  // If either of the source or tag are available from the recv (i.e., it is not
  // a wildcard receive) then update the channel 
  if ( recv.source != DUMPI_ANY_SOURCE ) {
    partial_channel.set_src( recv.source ); 
  } 
  if ( recv.tag != DUMPI_ANY_TAG ) {
    partial_channel.set_tag( recv.tag );
  }
  // If the channel is complete, just return it, otherwise we need to look at
  // the recv's status
  if ( partial_channel.is_complete() ) {
    return partial_channel;
  } 
  else {
    dumpi_status* status_ptr = recv.status;
    // Status wasn't ignored, we can recover the source and tag
    if ( status_ptr != nullptr && status_ptr != DUMPI_STATUS_IGNORE ) {
      dumpi_status status = *status_ptr;
      partial_channel.set_src( status.source );
      partial_channel.set_tag( status.tag );
      return partial_channel;
    }
    // Status was ignored, can't unambiguously determine message matching
    else {
      std::stringstream ss;
      ss << "Attempting to determine channel of receive, but status was ignored" 
         << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}

Channel Trace::determine_channel_of_irecv( const Request& request, 
                                           const dumpi_status* status )
{
  // Check if we can just look up the channel from the request
  if ( request.get_channel().is_complete() ) {
    return request.get_channel();
  }
  // If not, try to complete the channel from the additional information in the
  // status
  else {
    // Check nullity of status pointer
    if ( status != nullptr ) {
      // Check whether status was ignored
      if ( status != DUMPI_STATUS_IGNORE ) {
        Channel partial_channel = request.get_channel();
        int sender_rank = status->source;
        int message_tag = status->tag;
        Channel complete_channel( sender_rank,
                                  partial_channel.get_dst(),
                                  message_tag,
                                  partial_channel.get_comm() );
        return complete_channel;
      }
      else {
        std::stringstream ss;
        ss << "Attempting to determine channel of receive, but status was ignored" 
           << std::endl;
        throw std::runtime_error( ss.str() );
      }
    }
    else {
      std::stringstream ss;
      ss << "Attempting to determine channel of receive, but status is null" 
         << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Trace::cancel_request( long request_id )
{
  auto search = this->id_to_request.find( request_id );
  if ( search != this->id_to_request.end() ) {
    search->second.cancel();
  } else {
    std::stringstream ss;
    ss << "Rank: " << this->get_dumpi_to_graph_rank() 
       << " trying to cancel request ID: " << request_id
       << " but that request ID is not in id_to_request" << std::endl;
    throw std::runtime_error( ss.str() );
  }
}

void Trace::free_request( long request_id )
{
  auto search = this->id_to_request.find( request_id );
  if ( search != this->id_to_request.end() ) {
    this->id_to_request.erase( search->first );
  } else {
    std::stringstream ss;
    ss << "Rank: " << this->get_dumpi_to_graph_rank() 
       << " trying to free/remove request ID: " << request_id
       << " but that request ID is not in id_to_request" << std::endl;
    throw std::runtime_error( ss.str() );
  }
}

void Trace::complete_request( long request_id,
                              const dumpi_status* status_ptr,
                              const dumpi_time cpu_time,
                              const dumpi_time wall_time,
                              const dumpi_perfinfo *ctrs,
                              std::string matching_fn_call )
{
  // Generally speaking, a class of MPI functions called "matching functions"
  // are called to try to complete requests. Requests may be generated by 
  // non-blocking point-to-point communication functions like MPI_Irecv or 
  // MPI_Isend, non-blocking one-sided communication functions like MPI_Rput or
  // MPI_Rget, and various others. The outcomes of handling these requests are
  // modeled differently in the event graph, so the first thing we have to do
  // is look up what kind of request we are dealing with.
  // TODO: verify insertion
  auto request_search = this->id_to_request.find( request_id );
  int request_type;
  Request request;
  // Case 1: Request is currently tracked, so we can get its type and proceed
  // std::cerr << "JACK_ Request search: " <<  request_search->second.get_channel() << std::endl;
  if ( request_search != this->id_to_request.end() ) {
    request = request_search->second;
    request_type = request.get_type();
  } 
  // Case 2: Request is not currently tracked. It was either not added in the
  // appropriate callback, or the trace is malformed. Either way, we have to 
  // abort. 
  else {
    // std::cout << "JACK:: Request search end: " <<  request_id << " rank: " << this->get_trace_rank() << std::endl;
    std::stringstream ss;
    std::cout << "Tried to handle request corresponding to request ID: " 
       << request_id << " but no such request exists." << " rank: " << this->get_trace_rank()<< std::endl;
    // std::cerr << "Current id_to_request map is:" << std::endl;
    // this->report_id_to_request();
    throw std::runtime_error( ss.str() );
  }

  // Now that we know what kind of request we are handling, call the appropriate
  // specialized request completion handler.
  // Case 1: request_type == 0 --> MPI_Isend
  if ( request_type == 0 ) {
    // std::cout << "JACK_ isendRequest" << std::endl;
    complete_isend_request( request );

  }
  // Case 2: request_type == 1 --> MPI_Irecv
  else if ( request_type == 1 ) {
    // std::cout << "JACK_ irecvRequest" << std::endl;
    complete_irecv_request( request, 
                            status_ptr, 
                            cpu_time, 
                            wall_time,
                            ctrs,
                            matching_fn_call );
  }
  // Case 3+: FIXME: We don't handle persistent communication requests or 
  // one-sided communication just yet
  else {
    // std::cout << "JACK_ case 3+" << std::endl;

    std::stringstream ss;
    ss << "Handling completion of persistent communication requests and "
       << "one-sided communication requests is not implemented yet." << std::endl;
    throw std::runtime_error( ss.str() );
  }

  // No matter what kind of request we are handling, once it is handled it must
  // be removed from the id_to_request map.
  // std::cout << "JACK_ Erasing id:  " << request_id << std::endl;

  this->id_to_request.erase( request_id );
}

void Trace::complete_isend_request( Request request )
{
  // FIXME: Actually, this is not so much a FIXME as it is a WARNING
  // We don't explicitly represent isend request completions with vertices, but
  // we *do* bail on any trace involving cancelled isends because cancelling 
  // isends is pure evil and we cannot guarantee an accurate event graph in the
  // presence of that sort of heresy (seriously, look up what the standard says
  // about cancelling isends and tell me that you can be *absolutely* sure what
  // actually happened given the information in a DUMPI tracefile)
  if ( request.is_cancelled() ) {
    std::stringstream ss;
    ss << "Cancelled isend present in trace. Aborting." << std::endl;
    throw std::runtime_error( ss.str() );
  }
}

void Trace::complete_irecv_request( Request request,
                                    const dumpi_status* status_ptr,
                                    const dumpi_time cpu_time,
                                    const dumpi_time wall_time,
                                    const dumpi_perfinfo *ctrs,
                                    std::string matching_fn_call )
{
  // FIXME: We basically assume that cancelled irecv requests are always 
  // effectively cancelled and thus should not be represented by a recv vertex
  // in the event graph, even though the standard just says that they're 
  // either cancelled or complete normally. SUPER HELPFUL. 
  if ( !request.is_cancelled() ) {
    // Determine the channel in which this receive occurs
    Channel channel = this->determine_channel_of_irecv( request, 
                                                        status_ptr );
    // Create a recv event
    size_t event_vertex_id = this->get_next_vertex_id();

    // Associate this receive event with its channel
    this->register_recv( channel, event_vertex_id );

    // Associate this receive event with a timestamp
    this->register_dumpi_timestamp( wall_time );

    //Associate this receive event with papi struct if necessary
    if(ctrs != nullptr){
      this->register_papi_struct(*ctrs);
    }

    // Associate this receive event with the MPI matching function call that
    // generated it
    this->associate_event_with_call( matching_fn_call, event_vertex_id );
  }
}

void Trace::apply_vertex_id_offset( size_t offset ) 
{
  this->vertex_id_offset = offset;
  this->offset_set = true;
  this->initial_vertex_id += offset;
  this->final_vertex_id += offset;
  // Update the mapping from vertex IDs to channelsS
  std::unordered_map<size_t,Channel> new_vertex_id_to_channel;
  for ( auto kvp : this->vertex_id_to_channel ) {
    size_t old_vertex_id = kvp.first;
    size_t new_vertex_id = old_vertex_id + offset;
    Channel channel = kvp.second;
    new_vertex_id_to_channel.insert( { new_vertex_id, channel } );
  }
  this->vertex_id_to_channel = new_vertex_id_to_channel;
  // Update mapping from channels to sequences of vertex IDs representing sends
  for ( auto kvp : this->channel_to_send_seq ) {
    std::vector<size_t> new_vertex_ids;
    for ( int i=0; i<kvp.second.size(); ++i ) {
      new_vertex_ids.push_back( kvp.second[i] + offset );
    }
    this->channel_to_send_seq[ kvp.first ] = new_vertex_ids;
  }
  // Update mapping from channels to sequences of vertex IDs representing recvs
  for ( auto kvp : this->channel_to_recv_seq ) {
    std::vector<size_t> new_vertex_ids;
    for ( int i=0; i<kvp.second.size(); ++i ) {
      new_vertex_ids.push_back( kvp.second[i] + offset );
    }
    this->channel_to_recv_seq[ kvp.first ] = new_vertex_ids;
  }
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Dumb Getters //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Just gets the current vertex ID that's up for assignment without incrementing
// it
size_t Trace::get_curr_vertex_id() const 
{
  return this->curr_vertex_id;
}

// Returns the next vertex ID to assign to an event and increments the current
// vertex ID
size_t Trace::get_next_vertex_id() 
{
  size_t vertex_id = this->curr_vertex_id;
  this->curr_vertex_id++;
  return vertex_id;
} 

size_t Trace::get_initial_vertex_id() const
{
  return this->initial_vertex_id;
}
  
size_t Trace::get_final_vertex_id() const
{
  return this->final_vertex_id;
}

size_t Trace::get_vertex_id_offset() const
{
  return this->vertex_id_offset;
}

// Returns the directory (effectively, the set of trace files representing a 
// single execution) of which the trace file that this Trace object represents
// is a member
std::string Trace::get_trace_dir() const  
{
  return this->trace_dir;
}

// Returns the MPI rank of the application process that generated the trace file
// that this Trace object represents
int Trace::get_trace_rank() const
{
  return this->trace_rank;
}

// Returns the MPI rank of the dumpi_to_graph process handling the trace file 
// that this Trace object represents
int Trace::get_dumpi_to_graph_rank() const
{
  return this->dumpi_to_graph_rank;
}

// Returns the sequence of event types
std::vector<uint8_t> Trace::get_event_seq() const
{
  return this->event_seq;
}

// Returns the sequence of wall-time time stamps associated with events
std::vector<double> Trace::get_wall_time_seq() const
{
  return this->wall_time_seq;
}

std::vector<std::string> Trace::get_perf_counter_seq() const
{
  return this->counter_sets;
}

// Returns the sequence of (MPI_Function, function_call_index) that generated
// each event. This is used to associate other non-DUMPI trace data (e.g., 
// callstacks) to events in the event sequence
std::vector<std::pair<std::string,size_t>> Trace::get_mpi_fn_seq() const 
{
  return this->mpi_fn_seq;
}

// Returns the mapping from vertex IDs to channels
std::unordered_map<size_t,Channel> Trace::get_vertex_id_to_channel() const
{
  return this->vertex_id_to_channel;
}

// Returns the mapping between channels having this trace process as the 
// receiver and the sequence of vertex IDs representing those receives
channel_map Trace::get_channel_to_recv_seq() const
{
  return this->channel_to_recv_seq;
}

// Returns the mapping between channels having this trace process as the sender
// and the sequence of vertex IDs representing those sends
channel_map Trace::get_channel_to_send_seq() const
{
  return this->channel_to_send_seq;
}

std::unordered_map<long,Request> Trace::get_id_to_request() const
{
  return this->id_to_request;
}


////////////////////////////////////////////////////////////////////////////////
/////////// Convenience functions for printing the state of a Trace ////////////
////////////////////////////////////////////////////////////////////////////////

void Trace::report_event_seq()
{
  std::unordered_map<uint8_t, std::string> type_to_name =
  {
    {0, "send"}, {1, "recv"}, {2, "init"}, {3, "finalize"}, {4, "barrier"}
  };
  std::cout << "Event sequence for trace rank: " 
            << this->get_trace_rank() << std::endl;
  int n_events = this->event_seq.size();
  std::cout << "Number of events: " << n_events << std::endl;
  std::cout << "Current Vertex ID value: " << this->get_curr_vertex_id() << std::endl;

  //std::cout << "event_seq length: " << this->event_seq.size() << std::endl;
  //std::cout << "fn_call_seq length: " << this->_mpi_fn_seq.size() << std::endl;

  if ( this->offset_set ) {
    for ( int i=0; i<n_events; ++i ) {
      std::cout << "Vertex ID: " << i + this->vertex_id_offset
                << ", Type: " << type_to_name.at( this->event_seq[i] ) 
                << ", Generating Function: " << this->mpi_fn_seq[i].first
                << ", Function Call Index: " << this->mpi_fn_seq[i].second
                << std::endl;
    }
  }
  else {
    for ( int i=0; i<n_events; ++i ) {
      std::cout << "Vertex ID: " << i 
                << ", Type: " << type_to_name.at( this->event_seq[i] ) 
                << ", Generating Function: " << this->mpi_fn_seq[i].first
                << ", Function Call Index: " << this->mpi_fn_seq[i].second
                << std::endl;
    }
  }
}

void Trace::report_channel_to_send_seq()
{
  std::cout << "Channel to send sequence map for trace rank: " 
            << this->get_trace_rank() << std::endl;
  for ( auto kvp : this->channel_to_send_seq ) {
    std::cout << "Channel: " << kvp.first << ", Send Vertex IDs: ";
    for ( auto send : kvp.second ) {
      std::cout << " " << send;
    }
    std::cout << std::endl;
  }
}

void Trace::report_channel_to_recv_seq()
{
  std::cout << "Channel to recv sequence map for trace rank: " 
            << this->get_trace_rank() << std::endl;
  for ( auto kvp : this->channel_to_recv_seq ) {
    std::cout << "Channel: " << kvp.first << ", Recv Vertex IDs: ";
    for ( auto recv : kvp.second ) {
      std::cout << " " << recv;
    }
    std::cout << std::endl;
  }
}

void Trace::report_id_to_request()
{
  std::cout << "ID_to_request map for trace rank: " 
            << this->get_trace_rank() <<std::endl;

  for ( auto kvp : this->id_to_request ) {
    std::cout << "Request ID: " << kvp.first 
              << " Request Object: " << kvp.second << std::endl;
  }
}

void Trace::get_pluto_entry(int& msg_type, long& req_addr, int& source_type, long& entrynum){
  this->pluto_trace >> msg_type >> req_addr >> source_type >> entrynum;
}

