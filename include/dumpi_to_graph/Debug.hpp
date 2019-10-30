#ifndef D2G_DEBUG_H
#define D2G_DEBUG_H

#include <cassert>
#include <vector>

#include "Channel.hpp"
#include "Trace.hpp"

// When this is defined, various (possibly slow) safety checks are enabled
#define SAFETY

// When defined, we check insertions into maps that we *should* be able to 
// perform without these checks
#define PARANOID_INSERTION

// When defined, various runtime assertions will be checked at
#define SANITY_CHECK

// Checks that a vertex sequence (e.g., the sequence of sends or recvs mapped to 
// a channel) contains no negative values. While this does not ensure that the 
// vertex IDs are correct (that is largely the duty of the trace parsing 
// callbacks and the member functions of the Trace class) this at least can 
// detect when message buffers are corrupted by MPI wackiness during the 
// exchange of message matching data.
bool contains_no_invalid_vertex_ids( std::vector<size_t> vertex_ids );

bool validate_remote_recv_seqs( const int* recv_buffer, 
                                int recv_buffer_len,
                                const channel_map& chan_to_send_seq,
                                const std::unordered_map<Channel,int,ChannelHash>& chan_to_offset );

bool validate_trace( Trace trace );

#endif // D2G_DEBUG_H
