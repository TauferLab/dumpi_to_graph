#ifndef D2G_UTILITIES
#define D2G_UTILITIES

#include <string>

// DUMPI-specific
#include "dumpi/common/argtypes.h" 

int trace_file_to_rank( const std::string trace_file );

//FIXME: Currently does not check PAPI performance counter data
template<typename T>
bool validate_dumpi_event( T dumpi_event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr,
                           const dumpi_perfinfo* perf_counter_ptr )
{
  if ( dumpi_event_ptr != nullptr &&
       cpu_time_ptr    != nullptr &&
       wall_time_ptr   != nullptr ) {
    return true;
  } else {
    return false;
  }
}

#endif // D2G_UTILITIES
