#!/usr/bin/env bash
# Builds test/jsonl_format_test.cpp against igraph, runs spark/jsonl_to_graphml.py
# on the partitions it writes, and checks the result is byte-identical to
# igraph's own GraphML output.
#
# Usage: IGRAPH_ROOT=/path/to/igraph-0.7.1/install test/check_jsonl_format.sh
# Needs: a C++ compiler, igraph, and pyspark (python3 -c "import pyspark").
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
igraph_root=${IGRAPH_ROOT:?set IGRAPH_ROOT to the igraph install prefix}
python=${PYTHON:-python3}
work=$(mktemp -d)
trap 'rm -rf "${work}"' EXIT

${CXX:-c++} -std=c++11 -O1 \
  -I"${repo}/include/dumpi_to_graph" -I"${igraph_root}/include" \
  "${repo}/test/jsonl_format_test.cpp" \
  -L"${igraph_root}/lib" -Wl,-rpath,"${igraph_root}/lib" -ligraph \
  -o "${work}/jsonl_format_test"

"${work}/jsonl_format_test" "${work}/out"
"${python}" "${repo}/spark/jsonl_to_graphml.py" "${work}/out/event_graph_jsonl" \
  -o "${work}/out/actual.graphml" --partitions 3 > "${work}/spark.log" 2>&1 \
  || { cat "${work}/spark.log"; exit 1; }

if cmp "${work}/out/expected.graphml" "${work}/out/actual.graphml"; then
  echo "PASS: Spark-assembled GraphML is byte-identical to igraph's"
else
  diff <(cat -A "${work}/out/expected.graphml") <(cat -A "${work}/out/actual.graphml") | head -40
  echo "FAIL"
  exit 1
fi
