#ifndef D2G_CHANNEL_H
#define D2G_CHANNEL_H

#include <string>
#include <functional>
#include <vector> 
#include <iostream> 

// DUMPI
#include "dumpi/common/argtypes.h"

// Boost
#include "boost/functional/hash.hpp"
#include "boost/serialization/access.hpp"


class Channel {
public:
  // Constructors 
  Channel() : 
    src(-1), dst(-1), tag(-1), comm(-1) {}
  Channel(int dst, int comm) : 
    src(-1), dst(dst), tag(-1), comm(comm) {} 
  Channel(int dst, int tag, int comm) : 
    src(-1), dst(dst), tag(tag), comm(comm) {} 
  Channel(int src, int dst, int tag, int comm) : 
    src(src), dst(dst), tag(tag), comm(comm) {} 
  Channel( const Channel& rhs ) : 
    src(rhs.src), dst(rhs.dst), tag(rhs.tag), comm(rhs.comm) {}
  Channel( const int sender_rank, const dumpi_send send_event ) :
    src(sender_rank),
    dst(send_event.dest),
    tag(send_event.tag),
    comm(send_event.comm) {}
  Channel( const int sender_rank, const dumpi_isend isend_event ) :
    src(sender_rank),
    dst(isend_event.dest),
    tag(isend_event.tag),
    comm(isend_event.comm) {}
 
  // Accessors
  int get_src() const;
  int get_dst() const;
  int get_tag() const;
  int get_comm() const;

  // Operators
  Channel operator=(const Channel& rhs) 
  {
    if (this == &rhs) {
      return *this;
    } else {
      this->src = rhs.get_src();
      this->dst = rhs.get_dst();
      this->tag = rhs.get_tag();
      this->comm = rhs.get_comm();
      return *this;
    }
  }
  friend std::ostream& operator<< (std::ostream& out, const Channel& chan) 
  {
    out << "(" << chan.get_src() << ", " 
               << chan.get_dst() << ", " 
               << chan.get_tag() << ", " 
               << chan.get_comm() << ")"; 
    return out;
  }
  // Needed for use of Channel as key in std::unordered_map
  bool operator==(const Channel& c) const
  {
    return src  == c.get_src() && 
           dst  == c.get_dst() && 
           tag  == c.get_tag() && 
           comm == c.get_comm(); 
  }
  
  bool is_complete() 
  {
    if (src != -1 && dst != -1 && tag != -1 && comm != -1) {
      return true;
    } else {
      return false; 
    }
  }

  bool is_wildcard_src()
  {
    if (src == -1) {
      return true;
    } else {
      return false;
    }
  }

  bool is_wildcard_tag()
  {
    if (tag == -1) {
      return true;
    } else {
      return false;
    }
  }

  bool is_wildcard_src_and_tag()
  {
    if ( this->is_wildcard_src() && this->is_wildcard_tag() ) {
      return true;
    } else {
      return false;
    }
  }

private:
  // The source rank 
  int src;
  // The destination rank
  int dst;
  // The message tag
  int tag;
  // The communicator in which this channel exists
  int comm;

  // Needed so we can serialize containers of Channels and send them during
  // message matching 
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize( Archive& ar, const unsigned int version ) 
  {
    ar & src;
    ar & dst;
    ar & tag;
    ar & comm;
  }

};

// Hash function for Channels, needed for use of Channel as key in 
// std::unordered_map
struct ChannelHash {
  std::size_t operator() (Channel channel) const
  {
    std::size_t hash = 0;
    boost::hash_combine( hash, boost::hash_value( channel.get_src() ) );
    boost::hash_combine( hash, boost::hash_value( channel.get_dst() ) );
    boost::hash_combine( hash, boost::hash_value( channel.get_tag() ) );
    boost::hash_combine( hash, boost::hash_value( channel.get_comm() ) );
    return hash; 
  }
};


#endif // D2G_CHANNEL_H   
