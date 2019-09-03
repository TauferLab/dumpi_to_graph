#ifndef D2G_UTILITIES
#define D2G_UTILITIES

#include <string>
#include <sstream>
#include <stdexcept>

// DUMPI-specific
#include "dumpi/common/argtypes.h" 

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

////////////////////////////////////////////////////////////////////////////////

template<>
void validate_dumpi_event( dumpi_waitsome* event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr,
                           const dumpi_perfinfo* perf_counter_ptr )
{
  if ( event_ptr     == nullptr ||
       cpu_time_ptr  == nullptr ||
       wall_time_ptr == nullptr 
     ) {
    std::stringstream ss;
    ss << "corrupted waitsome data - exiting" << std::endl;
    throw std::runtime_error( ss.str() );
  } 
  else {
    if ( event_ptr->requests == nullptr || event_ptr->indices  == nullptr ) {
      std::stringstream ss;
      ss << "waitsome requests/indices are null" << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}

template<>
void validate_dumpi_event( dumpi_testsome* event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr,
                           const dumpi_perfinfo* perf_counter_ptr )
{
  if ( event_ptr     == nullptr ||
       cpu_time_ptr  == nullptr ||
       wall_time_ptr == nullptr 
     ) {
    std::stringstream ss;
    ss << "corrupted testsome data - exiting" << std::endl;
    throw std::runtime_error( ss.str() );
  } 
  else {
    if ( event_ptr->requests == nullptr || event_ptr->indices  == nullptr ) {
      std::stringstream ss;
      ss << "testsome requests/indices are null" << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

template<>
void validate_dumpi_event( dumpi_waitall* event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr,
                           const dumpi_perfinfo* perf_counter_ptr )
{
  if ( event_ptr     == nullptr ||
       cpu_time_ptr  == nullptr ||
       wall_time_ptr == nullptr 
     ) {
    std::stringstream ss;
    ss << "corrupted waitall data - exiting" << std::endl;
    throw std::runtime_error( ss.str() );
  } 
  else {
    if ( event_ptr->requests == nullptr ) {
      std::stringstream ss;
      ss << "waitall requests are null" << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}

template<>
void validate_dumpi_event( dumpi_waitany* event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr,
                           const dumpi_perfinfo* perf_counter_ptr )
{
  if ( event_ptr     == nullptr ||
       cpu_time_ptr  == nullptr ||
       wall_time_ptr == nullptr 
     ) {
    std::stringstream ss;
    ss << "corrupted waitany data - exiting" << std::endl;
    throw std::runtime_error( ss.str() );
  } 
  else {
    if ( event_ptr->requests == nullptr ) {
      std::stringstream ss;
      ss << "waitany requests are null" << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}

template<>
void validate_dumpi_event( dumpi_testall* event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr,
                           const dumpi_perfinfo* perf_counter_ptr )
{
  if ( event_ptr     == nullptr ||
       cpu_time_ptr  == nullptr ||
       wall_time_ptr == nullptr 
     ) {
    std::stringstream ss;
    ss << "corrupted testall data - exiting" << std::endl;
    throw std::runtime_error( ss.str() );
  } 
  else {
    if ( event_ptr->requests == nullptr ) {
      std::stringstream ss;
      ss << "testall requests are null" << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}

template<>
void validate_dumpi_event( dumpi_testany* event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr,
                           const dumpi_perfinfo* perf_counter_ptr )
{
  if ( event_ptr     == nullptr ||
       cpu_time_ptr  == nullptr ||
       wall_time_ptr == nullptr 
     ) {
    std::stringstream ss;
    ss << "corrupted testany data - exiting" << std::endl;
    throw std::runtime_error( ss.str() );
  } 
  else {
    if ( event_ptr->requests == nullptr ) {
      std::stringstream ss;
      ss << "testany requests are null" << std::endl;
      throw std::runtime_error( ss.str() );
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

#endif // D2G_UTILITIES
