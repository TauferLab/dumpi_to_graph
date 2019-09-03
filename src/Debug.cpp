#include "Debug.hpp"

#include <mpi.h>

#include <vector>

#include "Channel.hpp"
#include "Trace.hpp"

bool contains_no_invalid_vertex_ids( std::vector<size_t> vertex_ids )
{
  for ( auto vid : vertex_ids ) {
    if ( vid <= 0 ) {
      return false;
    }
  }
  return true;
}

bool validate_remote_recv_seqs( const int* recv_buffer, 
                                int recv_buffer_len,
                                const channel_map& chan_to_send_seq,
                                const std::unordered_map<Channel,int,ChannelHash>& chan_to_offset )
{
  for ( auto kvp : chan_to_send_seq ) {
    int offset = chan_to_offset.at( kvp.first );
    int n_sends = kvp.second.size();
    std::vector<size_t> recv_seq;
    recv_seq.reserve( n_sends );
    for ( int i=0; i<n_sends; ++i ) {
      int recv_buffer_idx = offset + i;
      assert( recv_buffer_idx < recv_buffer_len );
      recv_seq[i] = recv_buffer[ offset + i ];
    }
    assert( contains_no_invalid_vertex_ids( recv_seq ) );
  }
  return true;
}

bool validate_trace( Trace trace )
{
  // First lets check that the the event sequence begins with an init and ends 
  // with a finalize
  
  // Next, lets check that the channel maps and the mapping from vertex IDs to
  // channels contain vertex sequences with only valid values 
  auto channel_to_send_seq = trace.get_channel_to_send_seq();
  for ( auto kvp : channel_to_send_seq ) {
    assert( contains_no_invalid_vertex_ids( kvp.second ) );
  }
  auto channel_to_recv_seq = trace.get_channel_to_recv_seq();
  for ( auto kvp : channel_to_send_seq ) {
    assert( contains_no_invalid_vertex_ids( kvp.second ) );
  }

  //// Finally, lets check that the id_to_request map is empty
  //if ( !trace.get_id_to_request().empty() ) {
  //  int mpi_rc, rank;
  //  mpi_rc = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  //  trace.report_id_to_request();
  //  //exit(0);
  //}
  ////assert( trace.get_id_to_request().empty() );

  // Return OK
  return true;
}
