#ifndef D2G_CSMPI_TRACE_H
#define D2G_CSMPI_TRACE_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// Internal
#include "CSMPI_Callstack.hpp"
#include "Utilities.hpp"

class CSMPI_Trace
{
public:
  CSMPI_Trace( std::string csmpi_trace_file, int trace_file_rank );
  std::string lookup_callstack( std::string fn, int call_idx ) const;
private:
  std::string trace_file;
  int trace_rank;
  std::unordered_map<std::pair<std::string,int>,CSMPI_Callstack,pair_hash> fn_idx_pair_to_callstack;
};

#endif 
