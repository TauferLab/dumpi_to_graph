#ifndef D2G_UTILITIES
#define D2G_UTILITIES

#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>

// DUMPI
#include "dumpi/common/argtypes.h" 
#include "dumpi/common/types.h"

// Boost
#include <boost/functional/hash.hpp> 


std::string get_csmpi_trace_file( std::string trace_dir, int trace_rank );

int trace_file_to_rank( const std::string trace_file );

// FIXME: Currently does not check PAPI performance counter data
// FIXME: Tons of duplicate code in the explicit specializations below. 
// Structurally identical specializations are delimited by "///////.../////////"
// A reimplementation using type traits / concepts etc. would be highly welcome
template<typename T>
void validate_dumpi_event( T event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr )
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

template<typename T>
void validate_dumpi_event( T event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr,
                           const dumpi_perfinfo* perf_counter_ptr )
{
  if ( event_ptr     == nullptr ||
       cpu_time_ptr  == nullptr ||
       wall_time_ptr == nullptr ||
       perf_counter_ptr == nullptr
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

struct pair_hash
{
  template<typename T1, typename T2>
  std::size_t operator() (const std::pair<T1,T2> &pair) const 
  {
    std::size_t hash = 0;
    auto h1 = boost::hash_value(pair.first);
    auto h2 = boost::hash_value(pair.second);
    boost::hash_combine( hash, h1 );
    boost::hash_combine( hash, h2 );
    return hash;
  }
};

std::string stringify_perfinfo(dumpi_perfinfo counters);

#endif // D2G_UTILITIES
