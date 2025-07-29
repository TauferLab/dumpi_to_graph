#ifndef D2G_COL_CHANNEL
#define D2G_COL_CHANNEL

#include <string>
#include <functional>
#include <vector> 
#include <iostream> 
#include <unordered_map>

// DUMPI
#include "dumpi/common/argtypes.h"

// Boost
#include "boost/functional/hash.hpp"
#include "boost/serialization/access.hpp"

class CollectiveChannel{
    public:
    CollectiveChannel() = default;
    CollectiveChannel(int r, dumpi_comm c, int t): root(r), comm(c), type(t){}

    //accessors
    int get_root() const;
    dumpi_comm get_comm() const;
    int get_type() const;

    //setters
    void set_root(int);
    void set_comm(dumpi_comm);
    void set_type(int);

    bool operator==(const CollectiveChannel& c) const
    {
    return root == c.get_root() && comm == c.get_comm() && type == c.get_type();
    }

    private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int v){
      ar & root;
      ar & comm;
      ar & type;
    }
    int root;
    dumpi_comm comm;
    // Types: 
    // 1. MPI_Reduce
    // 2. MPI_BCast
    // 3. MPI_ALLReduce
    // 4. MPI_Alltoall
    int type;
};

struct CollectiveChannelHash {
  std::size_t operator() (CollectiveChannel channel) const
  {
    std::size_t hash = 0;
    boost::hash_combine( hash, boost::hash_value( channel.get_root() ) );
    boost::hash_combine( hash, boost::hash_value( channel.get_comm() ) );
    boost::hash_combine( hash, boost::hash_value( channel.get_type() ));
    return hash; 
  }
};

using collective_channel_map = std::unordered_map<CollectiveChannel, std::vector<size_t>, CollectiveChannelHash>;

#endif