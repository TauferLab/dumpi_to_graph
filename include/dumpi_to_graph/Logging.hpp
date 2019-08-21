#ifndef D2G_LOGGING_H
#define D2G_LOGGING_H

// Include this here because we almost always want to restrict logging to a
// single process, or at the very least tag the output of logging with which
// process it is coming from. 
#include <mpi.h>
// These should be fairly obvious...
#include <iostream>
#include <sstream>
#include <fstream>

#define REPORT_PROGRESS
#define REPORTING_RANK 0

#endif // D2G_LOGGING
