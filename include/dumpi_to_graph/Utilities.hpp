#ifndef D2G_UTILITIES
#define D2G_UTILITIES

#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>

// DUMPI
#include "dumpi/common/argtypes.h" 

// Boost
#include <boost/functional/hash.hpp> 

int trace_file_to_rank( const std::string trace_file );

// FIXME: Currently does not check PAPI performance counter data
// FIXME: Tons of duplicate code in the explicit specializations below. 
// Structurally identical specializations are delimited by "///////.../////////"
// A reimplementation using type traits / concepts etc. would be highly welcome
template<typename T>
void validate_dumpi_event( T event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr,
                           const dumpi_perfinfo* perf_counter_ptr )
{
  if ( event_ptr     == nullptr ||
       cpu_time_ptr  == nullptr ||
       wall_time_ptr == nullptr 
     ) {
    std::stringstream ss;
    ss << "corrupted event data - exiting" << std::endl;
    throw std::runtime_error( ss.str() );
  }
}

struct rank_seq_hash 
{
template <typename RankType, typename ElemType>
  std::size_t operator () (const std::pair<RankType,std::vector<ElemType>> &p) const 
  {
    std::size_t hash = 0;
    // hash the vector elements
    ElemType current_value;
    for (std::size_t i=0; i<p.second.size(); ++i) {
      current_value = p.second[i];
      boost::hash_combine( hash, boost::hash_value( current_value ) );
    }
    // hash the rank
    auto rank_hash = boost::hash_value( p.first );
    boost::hash_combine( hash, rank_hash );
    return hash;
  }
};


#endif // D2G_UTILITIES
