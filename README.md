## dumpi\_to\_graph

### Summary:
This tool converts sets of DUMPI execution traces to a graph representation of the execution. 

### Usage:
To run `dumpi_to_graph` on your execution trace data, do the following:

`mpirun -n ${n_processes} dumpi_to_graph.exe ${dumpi_to_graph_config.json} ${trace_directory}`

Depending on your environment, you may use your scheduler's `mpirun`-like equivalent, e.g., `srun`. 
Since `dumpi_to_graph` was developed on a system scheduled by SLURM, most examples in `scripts` use `srun`.


### Configuration Options:
`dumpi_to_graph` takes as its first argument the path to a JSON configuration file that specifies what kind of event graph model to generate.
An example configuration file is in `config/default.json` The following options can be modified:

* Which MPI functions to represent as vertices in the event graph
* Which happens-before orders to represent as edges in the event graph
* What vertex labels to assign
* What edge labels to assign
* Whether to represent unmatched tests at all
* If so, wether to merge consecutive unmatched tests into a single vertex 
* Whether to represent each receive separately or to merge receives completed by the same matching function into a single vertex

### Distributed graph construction (JSON Lines + Spark):
By default every `dumpi_to_graph` rank sends its part of the event graph to rank 0, which builds the whole graph in memory with igraph and writes `event_graph.graphml`.
For large traces rank 0 can run out of memory. Setting `"output_format": "jsonl"` in the configuration file avoids that gather:

1. Each rank writes its own part of the graph as JSON Lines partitions, and rank 0 writes a small index:
    ```
    <trace_dir>/event_graph_jsonl/_metadata.json            # schema + vertex-id range, edge counts, process ids per partition
    <trace_dir>/event_graph_jsonl/vertices/part-NNNNN.jsonl  # {"id":12,"event_type":"recv","logical_time":4,"wall_time":7.88e-05,"process_id":0,"mpi_function":"MPI_Recv"}
    <trace_dir>/event_graph_jsonl/edges/part-NNNNN.jsonl     # {"src":3,"dst":4,"order":0,"part":0,"seq":3}  (order 0 = program, 1 = message)
    ```
2. `spark/jsonl_to_graphml.py` reads the partitions in parallel, sorts vertices by id and edges into igraph's order, formats them as GraphML in parallel, and streams the pieces into `<trace_dir>/event_graph.graphml`.
   The file is byte-identical to the one the default path writes, so nothing downstream changes.

```
mpirun -n ${n_processes} dumpi_to_graph.exe config/dumpi_only_jsonl.json ${trace_directory}
spark-submit spark/jsonl_to_graphml.py ${trace_directory}/event_graph_jsonl
```
or both steps at once with `spark/build_graph_distributed.sh ${n_processes} dumpi_to_graph.exe config/dumpi_only_jsonl.json ${trace_directory}`.

Notes:
* Requires `pyspark` (`pip install -r spark/requirements.txt`) and Java. Without `spark-submit`, `python3 spark/jsonl_to_graphml.py ...` runs Spark in local mode.
* The trace directory must be on a filesystem the Spark executors and driver can all see (local disk in local mode, or the cluster's shared filesystem).
* The partitions can also be queried directly without building the GraphML file, e.g. `load_vertices(spark, dir)` / `load_edges(spark, dir)` in `spark/jsonl_to_graphml.py`, or by using `_metadata.json` to find the partition that holds a vertex id or process.
* `test/check_jsonl_format.sh` checks the JSON Lines encoding and the Spark job against igraph's GraphML writer on awkward values (XML special characters, raw bytes, NaN/inf).

### Dependencies:
* MPI
    * Ideally, an implementation from the 21st century 
    * Whatever application you want to trace is probably more demanding feature-wise than `dumpi_to_graph` is. 
* SST-DUMPI (https://github.com/sstsimulator/sst-dumpi)
* Boost >=1.70.0 (https://www.boost.org/users/history/version_1_70_0.html)
* JSON for Modern C++ (https://github.com/nlohmann/json)
* For `"output_format": "jsonl"` only: Apache Spark / PySpark >= 3.3 and Java

### Known Issues:
* SST-DUMPI may not build properly for your favorite MPI implementation
    * Known compatible implementations include:
        * MVAPICH2/2.3

### Copyright and License:

Copyright (c) 2021, Global Computing Lab

ANACIN-X is distributed under terms of the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) with LLVM Exceptions.

See [LICENSE](https://github.com/TauferLab/dumpi_to_graph/blob/master/LICENSE) for more details.
