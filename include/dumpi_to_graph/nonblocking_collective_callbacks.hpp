#ifndef D2G_NONBLOCKING_COLLECTIVE_CALLBACKS_H
#define D2G_NONBLOCKING_COLLECTIVE_CALLBACKS_H

// DUMPI
#include "dumpi/common/argtypes.h"

// Callbacks for nonblocking collectives. Unlike blocking collectives, these do
// not create a vertex at the call site; they register a pending collective
// (keyed by the DUMPI request id) that is completed when the matching
// MPI_Wait/Waitall/Test is encountered. comm/root/type are taken directly from
// the DUMPI event struct, so no Pluto coupling is required.

int cb_MPI_Ireduce(const dumpi_ireduce *prm, uint16_t thread, const dumpi_time *cpu, const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg);
int cb_MPI_Iallreduce(const dumpi_iallreduce *prm, uint16_t thread, const dumpi_time *cpu, const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg);
int cb_MPI_Ibcast(const dumpi_ibcast *prm, uint16_t thread, const dumpi_time *cpu, const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg);
int cb_MPI_Igather(const dumpi_igather *prm, uint16_t thread, const dumpi_time *cpu, const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg);
int cb_MPI_Igatherv(const dumpi_igatherv *prm, uint16_t thread, const dumpi_time *cpu, const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg);
int cb_MPI_Iscatter(const dumpi_iscatter *prm, uint16_t thread, const dumpi_time *cpu, const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg);
int cb_MPI_Iscatterv(const dumpi_iscatterv *prm, uint16_t thread, const dumpi_time *cpu, const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg);
int cb_MPI_Ialltoall(const dumpi_ialltoall *prm, uint16_t thread, const dumpi_time *cpu, const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg);
int cb_MPI_Ialltoallv(const dumpi_ialltoallv *prm, uint16_t thread, const dumpi_time *cpu, const dumpi_time *wall, const dumpi_perfinfo *perf, void *uarg);

#endif // D2G_NONBLOCKING_COLLECTIVE_CALLBACKS_H
