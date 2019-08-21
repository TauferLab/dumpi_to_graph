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
