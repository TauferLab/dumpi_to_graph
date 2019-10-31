#ifndef D2G_CSMPI_CALLSTACK_H
#define D2G_CSMPI_CALLSTACK_H

#include <string>
#include <vector>

class CSMPI_Callstack
{
public:
  CSMPI_Callstack( int idx, std::vector<std::string> addresses ) :
    _idx(idx), _addresses(addresses) {}
  std::string str() const;
private:
  int _idx;
  std::vector<std::string> _addresses;
};


#endif
