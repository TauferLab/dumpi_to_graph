#ifndef TRACE_H
#define TRACE_H

#include <string>
#include <vector>
#include <set>
#include <unordered_map>

// Boost 
#include "boost/mpi.hpp"

// DUMPI-specific 
#include "dumpi/common/argtypes.h" 

// Internal
#include "Configuration.hpp"

class Trace
{
public:
  Trace( d2g_Configuration config, 
         int trace_rank,
         int dumpi_to_graph_rank ) : 
    config( config ),
    trace_rank( trace_rank ),
    dumpi_to_graph_rank( dumpi_to_graph_rank )
    {}
 
  int get_trace_rank();
  int get_dumpi_to_graph_rank();

  void register_event( const dumpi_init init_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );
  
  void register_event( const dumpi_init_thread init_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );
  
  void register_event( const dumpi_finalize finalize_event, 
                       const dumpi_time cpu_time,
                       const dumpi_time wall_time );

private:
  d2g_Configuration config;
  int trace_rank;
  int dumpi_to_graph_rank;

};

#endif // TRACE_H
