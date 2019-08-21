#ifndef GLOB_H
#define GLOB_H

#include <vector>
#include <string>

// Wrapper function for C glob function
std::vector<std::string> glob(const std::string& pattern);

#endif // GLOB_H
