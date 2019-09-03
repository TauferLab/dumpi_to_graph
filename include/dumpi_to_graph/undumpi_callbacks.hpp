#ifndef UNDUMPI_CALLBACKS_H
#define UNDUMPI_CALLBACKS_H

#include "dumpi/libundumpi/libundumpi.h"
#include "Configuration.hpp" 

void set_callbacks(libundumpi_callbacks* callbacks, 
                   Configuration config);

#endif // UNDUMPI_CALLBACKS_H

