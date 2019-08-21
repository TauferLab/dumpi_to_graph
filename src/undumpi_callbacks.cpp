#include "undumpi_callbacks.hpp"

#include "Configuration.hpp" 
#include "init_callbacks.hpp"
#include "finalize_callbacks.hpp"
//#include "blocking_p2p_callbacks.hpp"
//#include "nonblocking_p2p_callbacks.hpp"
//#include "persistent_p2p_callbacks.hpp"
//#include "request_mutating_callbacks.hpp"
//#include "matching_function_callbacks.hpp" 
//#include "blocking_collectives_callbacks.hpp"
//#include "communicator_management_callbacks.hpp"

void set_callbacks(libundumpi_callbacks* callbacks, 
                   d2g_Configuration config)
{
  assert(callbacks != NULL); 

  // First, we always toggle on the callbacks for MPI_Init, MPI_Init_thread,
  // and MPI_Finalize. This is necessary to allow the event graph construction
  // to fail gracefully if given insufficient data to construct the graph 
  // according to the desired configuration
  callbacks->on_init        = cb_MPI_Init;
  callbacks->on_init_thread = cb_MPI_Init_thread;
  callbacks->on_finalize    = cb_MPI_Finalize;

  // Toggle on callbacks for all other functions specified in configuration
  for ( auto fn : config.get_mpi_functions() ) {
    //// Blocking point-to-point communication functions
    //if      (fn == "MPI_Send")         { callbacks->on_send         = cb_MPI_Send;         }
    //else if (fn == "MPI_Recv")         { callbacks->on_recv         = cb_MPI_Recv;         }
    //// Nonblocking point-to-point communication functions
    //else if (fn == "MPI_Isend")        { callbacks->on_isend        = cb_MPI_Isend;        }
    //else if (fn == "MPI_Irecv")        { callbacks->on_irecv        = cb_MPI_Irecv;        }
    //// Persistent point-to-point communication functions
    //else if (fn == "MPI_Recv_init")    { callbacks->on_recv_init    = cb_MPI_Recv_init;    }
    //else if (fn == "MPI_Send_init")    { callbacks->on_send_init    = cb_MPI_Send_init;    }
    ////// Request-mutating functions
    ////else if (fn == "MPI_Cancel")       { callbacks->on_cancel       = cb_MPI_Cancel;       }
    ////else if (fn == "MPI_Request_free") { callbacks->on_request_free = cb_MPI_Request_free; }
    ////else if (fn == "MPI_Start")        { callbacks->on_start        = cb_MPI_Start;        }
    ////else if (fn == "MPI_Startall")     { callbacks->on_startall     = cb_MPI_Startall;     }
    //// Blocking matching functions
    //else if (fn == "MPI_Wait")         { callbacks->on_wait         = cb_MPI_Wait;         }
    //else if (fn == "MPI_Waitany")      { callbacks->on_waitany      = cb_MPI_Waitany;      }
    //else if (fn == "MPI_Waitsome")     { callbacks->on_waitsome     = cb_MPI_Waitsome;     }
    //else if (fn == "MPI_Waitall")      { callbacks->on_waitall      = cb_MPI_Waitall;      }
    //// Non-blocking matching functions
    //else if (fn == "MPI_Test")         { callbacks->on_test         = cb_MPI_Test;         }
    //else if (fn == "MPI_Testsome")     { callbacks->on_testsome     = cb_MPI_Testsome;     }
    //else if (fn == "MPI_Testany")      { callbacks->on_testany      = cb_MPI_Testany;      }
    //else if (fn == "MPI_Testall")      { callbacks->on_testall      = cb_MPI_Testall;      }
    //// Blocking collectives
    //else if (fn == "MPI_Barrier")      { callbacks->on_barrier      = cb_MPI_Barrier;      }
    //// Communicator management
    //else if (fn == "MPI_Comm_rank")    { callbacks->on_comm_rank    = cb_MPI_Comm_rank;    }
    //else if (fn == "MPI_Comm_size")    { callbacks->on_comm_size    = cb_MPI_Comm_size;    }
    //else if (fn == "MPI_Comm_split")   { callbacks->on_comm_split   = cb_MPI_Comm_split;   }
  }
}




