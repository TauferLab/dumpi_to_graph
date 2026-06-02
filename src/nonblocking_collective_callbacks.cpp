#include "nonblocking_collective_callbacks.hpp"

// DUMPI
#include "dumpi/common/argtypes.h"

// Internal
#include "Trace.hpp"
#include "Utilities.hpp"
#include "CollectiveChannel.hpp"

// MPI_Ireduce: register a pending collective at the call site. The vertex is
// minted later, when the matching MPI_Wait/Test completes request==event.request
// (handled in matching_function_callbacks.cpp via complete_collective_request).
// No vertex and no Pluto entry are consumed here (mirrors how cb_MPI_Alltoall
// already skips Pluto).
int cb_MPI_Ireduce(const dumpi_ireduce *prm,
                   uint16_t thread,
                   const dumpi_time *cpu,
                   const dumpi_time *wall,
                   const dumpi_perfinfo *perf,
                   void *uarg)
{
  Trace* trace = (Trace*) uarg;
  bool papi = trace->get_papi_flag();
  if (papi) {
    validate_dumpi_event(prm, cpu, wall, perf);
  } else {
    validate_dumpi_event(prm, cpu, wall);
  }
  dumpi_ireduce event = *prm;

  // Type 1 == reduce (fan-in to root), matching the blocking MPI_Reduce code
  CollectiveChannel channel( event.root, event.comm, 1 );
  trace->register_collective_request( event.request, channel, "MPI_Ireduce" );

  return 0;
}

// Helper macro for the remaining nonblocking collectives. Each registers a
// pending collective keyed by the DUMPI request id, with the CollectiveChannel
// type code matching the blocking sibling (so make_collective_edges draws the
// same fan-in/fan-out topology). All-to-all collectives use root -1.
#define D2G_NBC_CB(FNNAME, DUMPITYPE, ROOTEXPR, TYPECODE)                       \
int FNNAME(const DUMPITYPE *prm, uint16_t thread, const dumpi_time *cpu,        \
           const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg) {    \
  Trace* trace = (Trace*) uarg;                                                 \
  bool papi = trace->get_papi_flag();                                          \
  if (papi) { validate_dumpi_event(prm, cpu, wall, perf); }                     \
  else { validate_dumpi_event(prm, cpu, wall); }                               \
  DUMPITYPE event = *prm;                                                       \
  CollectiveChannel channel( (ROOTEXPR), event.comm, (TYPECODE) );             \
  trace->register_collective_request( event.request, channel, #FNNAME + 3 );   \
  return 0;                                                                     \
}
// NOTE: #FNNAME + 3 strips the leading "cb_" to yield e.g. "MPI_Ibcast".
D2G_NBC_CB(cb_MPI_Iallreduce, dumpi_iallreduce, -1,          3)
D2G_NBC_CB(cb_MPI_Ibcast,     dumpi_ibcast,     event.root,  2)
D2G_NBC_CB(cb_MPI_Igather,    dumpi_igather,    event.root,  5)
D2G_NBC_CB(cb_MPI_Igatherv,   dumpi_igatherv,   event.root,  6)
D2G_NBC_CB(cb_MPI_Iscatter,   dumpi_iscatter,   event.root,  7)
D2G_NBC_CB(cb_MPI_Iscatterv,  dumpi_iscatterv,  event.root,  8)
D2G_NBC_CB(cb_MPI_Ialltoall,  dumpi_ialltoall,  -1,          4)
D2G_NBC_CB(cb_MPI_Ialltoallv, dumpi_ialltoallv, -1,          4)
#undef D2G_NBC_CB
