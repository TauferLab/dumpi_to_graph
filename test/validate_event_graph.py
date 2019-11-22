#!/usr/bin/env python3

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -o validate_event_graph-%j.out
#SBATCH -e validate_event_graph-%j.err

import argparse
import igraph

import pprint

from functools import wraps
from time import time
# A simple decorator for timing function calls
def timer(f):
    @wraps(f)
    def wrapper(*args, **kwargs):                                                  
        start = time()                                                             
        result = f(*args, **kwargs)                                                
        end = time()                                                               
        print("{} - Elapsed time: {}".format(f, end-start))
        return result                                                              
    return wrapper

# Reads in a single graph file via igraph
@timer
def read_graph( graph_path ):
    print("Reading in graph: {}".format(graph_path))
    graph = igraph.read( graph_path )
    return graph


@timer
def validate( graph, check_logical_timestamps, clock_condition ):
    # Check properties of init vertices
    print("Checking init vertex properties...")
    init_vertices = graph.vs.select( event_type_eq="init" )
    for v in init_vertices:
        preds = v.predecessors()
        succs = v.successors()
        assert( len(preds) == 0 ) 
        assert( len(succs) == 1 )
    # Check properties of finalize vertices
    print("Checking finalize vertex properties...")
    finalize_vertices = graph.vs.select( event_type_eq="finalize" )
    for v in finalize_vertices:
        preds = v.predecessors()
        succs = v.successors()
        assert( len(preds) == 1 ) 
        assert( len(succs) == 0 )
    # Check properties of barrier vertices
    print("Checking barrier vertex properties...")
    barrier_vertices = graph.vs.select( event_type_eq="barrier" )
    for v in barrier_vertices:
        preds = v.predecessors()
        succs = v.successors()
        assert( len(preds) == 1 ) 
        assert( len(succs) == 1 )
    # Check properties of recv vertices
    print("Checking recv vertex properties...")
    recv_vertices = graph.vs.select( event_type_eq="recv" )
    for v in recv_vertices:
        preds = v.predecessors()
        succs = v.successors()
        # Check that recv has two predecessors and one successor
        assert( len(preds) == 2 )
        assert( len(succs) == 1 )
        # Check that the one predecessor is in the same program order as this
        # recv (i.e., the local predecessor) and the other is a send in a 
        # different program order (i.e., the remote predecessor or sender)
        recv_pid = v["process_id"]
        assert( ( preds[0]["process_id"] != recv_pid and 
                  preds[0]["event_type"] == "send" and 
                  preds[1]["process_id"] == recv_pid ) 
                  or
                ( preds[1]["process_id"] != recv_pid and
                  preds[1]["event_type"] == "send" and 
                  preds[0]["process_id"] == recv_pid ) 
              )
        # Check that the recv's successor is in the same program order as it
        assert( succs[0]["process_id"] == recv_pid )
    # Check properties of send vertices
    print("Checking send vertex properties...")
    send_vertices = graph.vs.select( event_type_eq="send" )
    for v in send_vertices:
        preds = v.predecessors()
        succs = v.successors()
        # Check that send has one predecessors and two successors
        assert( len(preds) == 1 )
        assert( len(succs) == 2 )
        # Check that the one successor is in the same program order as this
        # send (i.e., the local successor) and the other is a send in a 
        # different program order (i.e., the remote successor or sender)
        send_pid = v["process_id"]
        assert( ( succs[0]["process_id"] != send_pid and 
                  succs[0]["event_type"] == "recv" and 
                  succs[1]["process_id"] == send_pid ) 
                  or
                ( succs[1]["process_id"] != send_pid and
                  succs[1]["event_type"] == "recv" and 
                  succs[0]["process_id"] == send_pid ) 
              )
        # Check that the send's predecessor is in the same program order as it
        assert( preds[0]["process_id"] == send_pid )
    # Check logical timestamp properties if desired
    if check_logical_timestamps:
        print("Checking logical timestamp properties...")
        for v in graph.vs[:]:
            lts = v["logical_time"]
            preds = v.predecessors()
            if len(preds) == 0:
                assert( lts == 0 )
            elif len(preds) == 1:
                pred_lts = preds[0]["logical_time"]
                if clock_condition == "lamport":
                    assert( lts == pred_lts + 1 )
                elif clock_condition == "scalar":
                    assert( lts > pred_lts )
            elif len(preds) == 2:
                pred_lts = [ p["logical_time"] for p in preds ]
                assert( lts == max(pred_lts) + 1 )

def main( graph_path, check_logical_timestamps, clock_condition ):
    graph = read_graph( graph_path )    
    validate( graph, check_logical_timestamps, clock_condition ) 




if __name__ == "__main__":
    desc = "Checks that an event graph is constructed properly"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("graph",
                        help="A GraphML file")
    parser.add_argument("-l", "--check_logical_timestamps", action="store_true", default=False,
                        help="Check that logical timestamps obey some given logical clock conditions")
    parser.add_argument("-c", "--clock_condition", type=str, default="scalar",
                        help="Specify which logical clock conditions to check. Options: \"scalar\", \"lamport\"")

    args = parser.parse_args()

    main( args.graph, args.check_logical_timestamps, args.clock_condition ) 
