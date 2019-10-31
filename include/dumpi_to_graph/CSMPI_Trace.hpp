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

class CSMPI_Trace
{
public:
  CSMPI_Trace( std::string csmpi_trace_file, int trace_file_rank );
private:
  std::string trace_file;
  int trace_rank;
  std::unordered_map<std::string,
                     std::vector<std::pair<int,CSMPI_Callstack>>> fn_to_callstack_seq;
};

#endif 
