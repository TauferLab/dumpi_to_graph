#include "Utilities.hpp"

#include <sstream> 
#include <stdexcept>
#include <string>
#include <regex>


int trace_file_to_rank( const std::string trace_file )
{
  std::regex rgx(".*-(\\d+)\\.bin");
  std::smatch match;
  if ( std::regex_search( trace_file.begin(), trace_file.end(), match, rgx ) ) {
    return std::atoi( match.str(1).c_str() );
  } else {
    std::stringstream ss;
    ss << "Could not extract rank from tracefile: " << trace_file << std::endl;
    throw std::runtime_error( ss.str() );
  }
}

std::string get_csmpi_trace_file( std::string trace_dir, int trace_rank )
{
  std::string csmpi_trace_file = trace_dir;
  csmpi_trace_file += "/csmpi/rank_";
  csmpi_trace_file += std::to_string( trace_rank );
  csmpi_trace_file += ".csmpi";
  return csmpi_trace_file;
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
       wall_time_ptr == nullptr || 
       perf_counter_ptr == nullptr 
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
       wall_time_ptr == nullptr || 
       perf_counter_ptr == nullptr 
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
       wall_time_ptr == nullptr ||
       perf_counter_ptr == nullptr 
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
       wall_time_ptr == nullptr || 
       perf_counter_ptr == nullptr 
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
       wall_time_ptr == nullptr || 
       perf_counter_ptr == nullptr 
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
       wall_time_ptr == nullptr ||
       perf_counter_ptr == nullptr 
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

template<>
void validate_dumpi_event( dumpi_waitsome* event_ptr,
                           const dumpi_time* cpu_time_ptr,
                           const dumpi_time* wall_time_ptr)
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
                           const dumpi_time* wall_time_ptr )
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
                           const dumpi_time* wall_time_ptr )
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
                           const dumpi_time* wall_time_ptr )
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
                           const dumpi_time* wall_time_ptr )
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
                           const dumpi_time* wall_time_ptr )
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

///////////////////////////////////////////////////////////////////////////////

std::string stringify_perfinfo(dumpi_perfinfo counters){
  std::ostringstream oss; 
  for (size_t i = 0; i < counters.count; i++){
    oss << counters.counter_tag[i] << " " << counters.outvalue[i];
  }
  oss << std::endl;
  return oss.str();
}
