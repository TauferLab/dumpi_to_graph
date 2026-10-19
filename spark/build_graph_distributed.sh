#!/usr/bin/env bash
# Distributed replacement for ANACIN-X's workflow_scripts/build_graph.sh.
#
# 1. dumpi_to_graph (MPI) parses the traces and each rank writes its own part
#    of the event graph as JSON Lines under <trace_dir>/event_graph_jsonl/.
#    The config must contain "output_format": "jsonl".
# 2. Spark sorts and formats the parts and writes <trace_dir>/event_graph.graphml,
#    byte-identical to the single-rank merge() output.
#
# Usage: build_graph_distributed.sh <n_procs> <dumpi_to_graph_bin> <config.json> <trace_dir>
# Environment:
#   MPIRUN          launcher (default: mpirun)
#   SPARK_SUBMIT    spark-submit to use (default: spark-submit on PATH)
#   SPARK_ARGS      extra spark-submit arguments, e.g. "--master local[8] --driver-memory 2g"
set -euo pipefail

n_procs=$1
dumpi_to_graph_bin=$2
dumpi_to_graph_config=$3
trace_dir=$4
here=$(cd "$(dirname "$0")" && pwd)

if ! grep -q '"output_format"[[:space:]]*:[[:space:]]*"jsonl"' "${dumpi_to_graph_config}"; then
  echo "error: ${dumpi_to_graph_config} must set \"output_format\": \"jsonl\"" >&2
  exit 1
fi

${MPIRUN:-mpirun} -np "${n_procs}" "${dumpi_to_graph_bin}" "${dumpi_to_graph_config}" "${trace_dir}"
# shellcheck disable=SC2086
${SPARK_SUBMIT:-spark-submit} ${SPARK_ARGS:-} "${here}/jsonl_to_graphml.py" \
  "${trace_dir}/event_graph_jsonl" -o "${trace_dir}/event_graph.graphml"
