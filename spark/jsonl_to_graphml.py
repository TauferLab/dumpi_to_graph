#!/usr/bin/env python3
"""Assemble a dumpi_to_graph event graph from JSON Lines partitions with Spark.

dumpi_to_graph run with "output_format": "jsonl" writes, per dumpi_to_graph
rank, its partial view of the event graph to

    <trace_dir>/event_graph_jsonl/_metadata.json
    <trace_dir>/event_graph_jsonl/vertices/part-NNNNN.jsonl
    <trace_dir>/event_graph_jsonl/edges/part-NNNNN.jsonl

instead of gathering everything on rank 0. This job is the "reduce" side: it
reads the partitions in parallel, sorts vertices by id and edges by the order
rank 0's merge() would have added them to igraph, formats each element as
GraphML in parallel, and streams the sorted pieces into one file. The result
is byte-identical to the event_graph.graphml written by the original
merge() + write() path (igraph 0.7.1 GraphML writer), so downstream ANACIN-X
tooling does not change. No step holds the whole graph in one process.

Usage:
    spark-submit spark/jsonl_to_graphml.py <trace_dir>/event_graph_jsonl \
        [-o <trace_dir>/event_graph.graphml] [--partitions N] [--keep-parts]

    (plain `python3 spark/jsonl_to_graphml.py ...` also works and runs Spark
    in local mode.)

The temporary part files are written next to the output file, so the output
directory must be on a filesystem that the Spark executors and the driver can
all see (true for local mode and for the shared filesystems HPC clusters use).

The module also exposes load_vertices()/load_edges() for reading the
partitions as DataFrames, e.g. to look up a few vertices without building the
GraphML file at all.
"""

import argparse
import json
import math
import os
import shutil
import sys
import uuid

FORMAT_NAME = "dumpi_to_graph.event_graph.jsonl"
SUPPORTED_VERSIONS = (1,)

# Strings in the partitions are raw bytes from the traces (callstacks, PAPI
# counters). Decoding with surrogateescape and encoding back the same way
# keeps any non-UTF-8 byte exactly as igraph would have written it.
ENCODING = "utf-8"
ERRORS = "surrogateescape"

GRAPHML_HEADER = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<graphml xmlns="http://graphml.graphdrawing.org/xmlns"\n'
    '         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n'
    '         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns\n'
    '         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">\n'
    "<!-- Created by igraph -->\n"
)
GRAPHML_FOOTER = "  </graph>\n</graphml>\n"


# --------------------------------------------------------------------------
# Formatting helpers. These reproduce igraph 0.7.1's igraph_write_graph_graphml
# (src/foreign-graphml.c) exactly.
# --------------------------------------------------------------------------

def xml_escape(s):
    """igraph_i_xml_escape: escapes & < > " ' and nothing else."""
    return (s.replace("&", "&amp;")
             .replace("<", "&lt;")
             .replace(">", "&gt;")
             .replace('"', "&quot;")
             .replace("'", "&apos;"))


def to_double(value):
    """JSON value -> the double igraph stored (non-finite values are strings)."""
    if isinstance(value, str):
        return float(value)  # "nan", "inf", "-inf"
    return float(value)


def format_double(x):
    """C printf("%g"); Python's %-formatting follows the same rules."""
    return "%g" % x


def key_lines(attributes):
    lines = []
    for attr in attributes:
        name = xml_escape(attr["name"])
        lines.append('  <key id="%s" for="node" attr.name="%s" attr.type="%s"/>\n'
                     % (name, name, attr["type"]))
    return "".join(lines)


def format_vertex(record, attributes):
    """One <node> element, without its final newline (_write_parts adds it)."""
    parts = ['    <node id="n%d">\n' % record["id"]]
    for attr in attributes:
        name = attr["name"]
        value = record.get(name)
        if attr["type"] == "double":
            x = to_double(value if value is not None else 0)
            if math.isnan(x):
                continue  # igraph leaves NaN numeric attributes out
            parts.append('      <data key="%s">%s</data>\n'
                         % (xml_escape(name), format_double(x)))
        else:
            parts.append('      <data key="%s">%s</data>\n'
                         % (xml_escape(name), xml_escape(value or "")))
    parts.append("    </node>")
    return "".join(parts)


def format_edge(record):
    return ('    <edge source="n%d" target="n%d">\n    </edge>'
            % (record["src"], record["dst"]))


# --------------------------------------------------------------------------
# Reading the partitions
# --------------------------------------------------------------------------

def read_metadata(input_dir):
    path = os.path.join(input_dir, "_metadata.json")
    if not os.path.exists(path):
        raise SystemExit(
            "No _metadata.json in %s. Was dumpi_to_graph run with "
            '"output_format": "jsonl", and did it finish?' % input_dir)
    with open(path, encoding=ENCODING) as f:
        meta = json.load(f)
    if meta.get("format") != FORMAT_NAME or meta.get("version") not in SUPPORTED_VERSIONS:
        raise SystemExit("%s is not a supported event graph partition index "
                         "(format=%r, version=%r)"
                         % (path, meta.get("format"), meta.get("version")))
    return meta


def partition_paths(input_dir, meta, kind):
    """Paths listed in the index (never a directory glob, so stale part files
    from an earlier run with more ranks are ignored)."""
    paths = [os.path.join(input_dir, p[kind]) for p in meta["partitions"]]
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise SystemExit("Missing partition files: %s" % ", ".join(missing))
    return [os.path.abspath(p) for p in paths]


# Spark's own text reader and writer go through Java strings, which would
# replace non-UTF-8 bytes. The partitions are therefore read (and the GraphML
# pieces written) by plain Python file I/O inside the tasks, split into byte
# ranges the same way Hadoop splits text files.
SPLIT_BYTES = 64 * 1024 * 1024


def _splits(paths, split_bytes):
    for path in paths:
        size = os.path.getsize(path)
        for start in range(0, max(size, 1), split_bytes):
            yield (path, start, min(start + split_bytes, size))


def _read_split(split):
    """Lines that start inside [start, end) of the file; a line that crosses
    the end belongs to this split, one that crosses the start does not."""
    path, start, end = split
    with open(path, "rb") as f:
        if start > 0:
            f.seek(start - 1)
            f.readline()  # finish the line that started before this split
        while f.tell() < end:
            raw = f.readline()
            if not raw:
                break
            raw = raw.rstrip(b"\n")
            if raw:
                yield json.loads(raw.decode(ENCODING, ERRORS))


def read_records(sc, paths, split_bytes=SPLIT_BYTES):
    """RDD of parsed JSON records from the given JSON Lines files."""
    splits = list(_splits(paths, split_bytes))
    return (sc.parallelize(splits, max(1, len(splits)))
              .flatMap(_read_split))


def load_vertices(spark, input_dir):
    """Vertices as a DataFrame (one row per vertex), for on-demand queries.
    (Spark's JSON reader replaces any non-UTF-8 bytes in string attributes;
    use read_records() when those bytes matter.)"""
    meta = read_metadata(input_dir)
    return spark.read.json(partition_paths(input_dir, meta, "vertices"))


def load_edges(spark, input_dir):
    """Edges as a DataFrame with columns src, dst, order (0 program,
    1 message), part, seq."""
    meta = read_metadata(input_dir)
    return spark.read.json(partition_paths(input_dir, meta, "edges"))


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------

def _encode(s):
    return s.encode(ENCODING, ERRORS)


def _write_parts(rdd, part_dir):
    """Write each partition of an RDD of bytes to part_dir/part-NNNNN, one
    element per line. Renamed into place, so a retried task cannot leave a
    half-written part behind."""
    def write(index, items):
        path = os.path.join(part_dir, "part-%05d" % index)
        tmp = "%s.%s.tmp" % (path, uuid.uuid4().hex[:8])
        with open(tmp, "wb") as f:
            for item in items:
                f.write(item)
                f.write(b"\n")
        os.replace(tmp, path)
        yield index
    os.makedirs(part_dir, exist_ok=True)
    rdd.mapPartitionsWithIndex(write).count()


def _concat(out, part_dir):
    names = sorted(n for n in os.listdir(part_dir)
                   if n.startswith("part-") and not n.endswith(".tmp"))
    for name in names:
        with open(os.path.join(part_dir, name), "rb") as f:
            shutil.copyfileobj(f, out, 16 * 1024 * 1024)


def build_graphml(spark, input_dir, output_path, n_partitions=None, keep_parts=False):
    sc = spark.sparkContext
    meta = read_metadata(input_dir)
    attributes = meta["vertex_attributes"]
    n_vertices = int(meta["n_vertices"])
    n_edges = int(meta["n_edges"])
    if n_partitions is None:
        # Roughly 1M vertices per output partition, at least one per core
        n_partitions = max(sc.defaultParallelism, (n_vertices + n_edges) // 1000000 + 1)

    vertices = read_records(sc, partition_paths(input_dir, meta, "vertices"))
    edges = read_records(sc, partition_paths(input_dir, meta, "edges"))

    # --- validate: vertex ids must be exactly 0..N-1 (igraph uses them as
    # indices), and every edge must point at an existing vertex
    vertex_ids = vertices.map(lambda r: r["id"]).cache()
    count = vertex_ids.count()
    if count != n_vertices:
        raise SystemExit("Index says %d vertices but partitions hold %d"
                         % (n_vertices, count))
    if count:
        lo, hi = vertex_ids.min(), vertex_ids.max()
        distinct = vertex_ids.distinct().count()
        if lo != 0 or hi != count - 1 or distinct != count:
            raise SystemExit("Vertex ids are not a contiguous 0..N-1 range "
                             "(min=%d, max=%d, distinct=%d, count=%d)"
                             % (lo, hi, distinct, count))
    vertex_ids.unpersist()

    bad_edges = edges.filter(
        lambda r: not (0 <= r["src"] < n_vertices and 0 <= r["dst"] < n_vertices)).count()
    if bad_edges:
        raise SystemExit("%d edges reference vertices outside 0..%d"
                         % (bad_edges, n_vertices - 1))

    # --- map: format every element in parallel; sort: restore igraph's order
    attrs = attributes  # captured by the closures below
    vertex_xml = (vertices
                  .map(lambda r: (r["id"], format_vertex(r, attrs)))
                  .sortByKey(numPartitions=n_partitions)
                  .map(lambda kv: _encode(kv[1])))
    edge_xml = (edges
                .map(lambda r: ((r["order"], r["part"], r["seq"]), format_edge(r)))
                .sortByKey(numPartitions=n_partitions)
                .map(lambda kv: _encode(kv[1])))

    output_path = os.path.abspath(output_path)
    out_dir = os.path.dirname(output_path)
    os.makedirs(out_dir, exist_ok=True)
    work_dir = os.path.join(out_dir, "_graphml_parts_" + uuid.uuid4().hex[:8])
    vertex_parts = os.path.join(work_dir, "vertices")
    edge_parts = os.path.join(work_dir, "edges")
    _write_parts(vertex_xml, vertex_parts)
    _write_parts(edge_xml, edge_parts)

    # --- reduce: stream header + sorted parts + footer into one file
    tmp_output = output_path + ".tmp"
    with open(tmp_output, "wb") as out:
        out.write(_encode(GRAPHML_HEADER))
        if n_vertices:
            out.write(_encode(key_lines(attributes)))
        out.write(_encode('  <graph id="G" edgedefault="%s">\n'
                          % ("directed" if meta.get("directed", True) else "undirected")))
        _concat(out, vertex_parts)
        _concat(out, edge_parts)
        out.write(_encode(GRAPHML_FOOTER))
    os.replace(tmp_output, output_path)

    if keep_parts:
        print("Kept GraphML part files in %s" % work_dir)
    else:
        shutil.rmtree(work_dir, ignore_errors=True)
    print("Wrote event graph to: %s (%d vertices, %d edges)"
          % (output_path, n_vertices, n_edges))
    return output_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("input_dir",
                        help="directory written by dumpi_to_graph with "
                             'output_format "jsonl" (<trace_dir>/event_graph_jsonl)')
    parser.add_argument("-o", "--output",
                        help="GraphML file to write (default: event_graph.graphml "
                             "next to input_dir, i.e. in the trace directory)")
    parser.add_argument("--partitions", type=int, default=None,
                        help="number of sorted output partitions")
    parser.add_argument("--keep-parts", action="store_true",
                        help="keep the intermediate GraphML part files")
    args = parser.parse_args(argv)

    input_dir = os.path.abspath(args.input_dir.rstrip("/"))
    output = args.output or os.path.join(os.path.dirname(input_dir), "event_graph.graphml")

    from pyspark.sql import SparkSession
    spark = SparkSession.builder.appName("dumpi_to_graph jsonl_to_graphml").getOrCreate()
    try:
        build_graphml(spark, input_dir, output, args.partitions, args.keep_parts)
    finally:
        spark.stop()


if __name__ == "__main__":
    sys.exit(main())
