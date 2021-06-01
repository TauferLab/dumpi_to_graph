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

### Dependencies:
* MPI
    * Ideally, an implementation from the 21st century 
    * Whatever application you want to trace is probably more demanding feature-wise than `dumpi_to_graph` is. 
* SST-DUMPI (https://github.com/sstsimulator/sst-dumpi)
* Boost >=1.70.0 (https://www.boost.org/users/history/version_1_70_0.html)
* JSON for Modern C++ (https://github.com/nlohmann/json)

### Known Issues:
* SST-DUMPI may not build properly for your favorite MPI implementation
    * Known compatible implementations include:
        * MVAPICH2/2.3

## Copyright and License:

Copyright (c) 2021, Global Computing Lab

ANACIN-X is distributed under terms of the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) with LLVM Exceptions.

See [LICENSE](https://github.com/TauferLab/ANACIN-X/blob/documentation/LICENSE) for more details.
