#include "Glob.hpp"

// glob_t, glob, globfree, GLOB_TILDE
#include <glob.h> 
// memset
#include <string.h>

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept> 
#include <iostream> 

// glob implementation taken from: 
// https://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
std::vector<std::string> glob(const std::string& pattern)
{
  glob_t glob_result;
  memset(&glob_result, 0, sizeof(glob_result));
  int rc;
  rc = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
  // Handle case where glob finds no matching files
  if (rc != 0) {
    globfree(&glob_result);
    std::stringstream ss;
    ss << "glob() failed with return value = " << rc << std::endl;
    throw std::runtime_error(ss.str());
  }
  // Otherwise, store the matching files
  std::vector<std::string> filenames;
  for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
    filenames.push_back(std::string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  return filenames; 
} 
