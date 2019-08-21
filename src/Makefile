MPICXX=mpicxx
RM=rm -f

CPPFLAGS= -Wall -Werror -pedantic 
CXXFLAGS= -std=c++11 -O3

INCLUDES= -I. 
INCLUDES+= -I../include
INCLUDES+= -I/g/g17/chapp1/repos/dumpi_to_graph/include/dumpi_to_graph/
INCLUDES+= -I/g/g17/chapp1/repos/dumpi_to_graph/include/external/
INCLUDES+= -I/g/g17/chapp1/repos/sst-dumpi/
INCLUDES+= -I/g/g17/chapp1/repos/boost_1_67_0/

LDFLAGS= 
LD_LIBRARY_FLAGS=

LIBS= -lundumpi 
LIBS+= -lboost_mpi 
LIBS+= -lboost_serialization 
LIBS+= -lboost_filesystem
LIBS+= -lboost_system

# Application 
SOURCES= dumpi_to_graph.cpp 
# CLI and configuration 
SOURCES+= Configuration.cpp 
SOURCES+= Glob.cpp 
# Callbacks for parsing event stream
SOURCES+= undumpi_callbacks.cpp
SOURCES+= init_callbacks.cpp
SOURCES+= finalize_callbacks.cpp
# Trace and event graph representation
SOURCES+= Trace.cpp
# Utilities
SOURCES+= Utilities.cpp

OBJECTS=$(SOURCES:.cpp=.o)

TARGET=dumpi_to_graph.exe

all: $(TARGET)

%.o: %.cpp
	$(MPICXX) -c $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $< -o $@

$(TARGET): $(OBJECTS)
	$(MPICXX) $(LDFLAGS) -o $@ $(OBJECTS) $(LD_LIBRARY_FLAGS) $(LIBS)


clean:
	$(RM) *.o
	$(RM) *.exe
	$(RM) *.out
	$(RM) *.err
	$(RM) *.core
