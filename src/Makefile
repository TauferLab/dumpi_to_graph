MPICXX=mpicxx
RM=rm -f

CPPFLAGS= -Wall #-Werror -pedantic 
CXXFLAGS= -std=c++11 -O3

INCLUDES= -I. 
INCLUDES+= -I../include
#INCLUDES+= -I/g/g17/chapp1/repos/dumpi_to_graph/include/dumpi_to_graph/
#INCLUDES+= -I/g/g17/chapp1/repos/dumpi_to_graph/include/external/
INCLUDES+= -I../include/dumpi_to_graph/
INCLUDES+= -I../include/external/
INCLUDES+= -I/g/g17/chapp1/repos/sst-dumpi/
#INCLUDES+= -I/g/g17/chapp1/repos/boost_1_67_0/
#INCLUDES+= -I/usr/tce/packages/boost/boost-1.69.0-mvapich2-2.3-gcc-8.1.0/include
INCLUDES+= -I/usr/tce/packages/boost/boost-1.69.0-mvapich2-2.3-gcc-4.9.3/include
INCLUDES+= -I/g/g17/chapp1/repos/igraph-0.7.1/build/include/

LDFLAGS= 
LD_LIBRARY_FLAGS= -L/usr/tce/packages/boost/boost-1.69.0-mvapich2-2.3-gcc-4.9.3/lib -L/g/g17/chapp1/repos/sst-dumpi/build_catalyst/lib -L/g/g17/chapp1/repos/igraph-0.7.1/build/lib/

LIBS= -lundumpi 
LIBS+= -lboost_mpi 
LIBS+= -lboost_serialization 
LIBS+= -lboost_filesystem
LIBS+= -lboost_system
LIBS+= -ligraph

# Application 
SOURCES= dumpi_to_graph.cpp 
# CLI and configuration 
SOURCES+= Configuration.cpp 
SOURCES+= Glob.cpp 
# Callbacks for parsing event stream
SOURCES+= undumpi_callbacks.cpp
SOURCES+= init_callbacks.cpp
SOURCES+= finalize_callbacks.cpp
SOURCES+= blocking_p2p_callbacks.cpp
SOURCES+= nonblocking_p2p_callbacks.cpp
SOURCES+= matching_function_callbacks.cpp
SOURCES+= request_mutating_callbacks.cpp
SOURCES+= blocking_collective_callbacks.cpp
SOURCES+= communicator_management_callbacks.cpp
# Trace and event graph representation
SOURCES+= Trace.cpp
SOURCES+= Channel.cpp
SOURCES+= Request.cpp
SOURCES+= EventGraph.cpp
# Utilities
SOURCES+= Utilities.cpp
SOURCES+= Debug.cpp

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
