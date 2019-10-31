#include "CSMPI_Callstack.hpp"

#include <string>
#include <sstream>
#include <vector>

// Returns a string representation of the call stack so that it can be added as
// a vertex label to the event graph
std::string CSMPI_Callstack::str() const
{
  std::ostringstream oss;
  //oss << std::hex;
  auto n_addresses = _addresses.size();
  for ( size_t i=0; i < n_addresses; ++i ) {
    oss << _addresses[i];
    if ( i != n_addresses - 1 ) {
      oss << ", ";
    }
  }
  oss << std::endl;
  return oss.str();
}
