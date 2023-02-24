#include "Request.hpp"

// Internal
#include "Channel.hpp"
// Request types:
// 0 - MPI_Isend
// 1 - MPI_Irecv
// 2 - MPI_Send_init
// 3 - MPI_Recv_init


Request::Request( int type, long id )
{
  this->type = type;
  this->id = id;
  std::cout << "Request ID (1): " << id << std::endl;
  if ( type == 0 || type == 1 ) {
    this->active = true;
  } else {
    this->active = false;
  }
  // If the full channel associated with this request is not available at the 
  // time of construction (e.g., when constructing the request for a wildcard
  // MPI_Irecv) default contruct a channel anyway
  this->channel = Channel(); 
}

Request::Request( int type, long id, const Channel& channel )
{
  this->type = type;
  this->id = id;
  // TODO: JACK_ create id_fk
  this->id_fk = id+type; // JACK_ MODIFY
  
  if ( type == 0 || type == 1 ) {
    this->active = true;
  } else {
    this->active = false;
  }
  this->channel = channel;
  // std::cout << "JACK:: Request ID_fk: " << this->id_fk << " ID: " << id << std::endl;
}

int Request::get_type() const
{
  return this->type;
}

long Request::get_id() const
{
  return this->id;
}

long Request::get_id_fk() const
{
  return this->id_fk;
}

Channel Request::get_channel() const
{
  return this->channel;
}

void Request::set_channel( const Channel& channel )
{
  this->channel = channel;
}

bool Request::is_active() const
{
  return this->active;
}

bool Request::is_cancelled() const
{
  return this->cancelled;
}

void Request::activate() 
{
  this->active = true;
}

void Request::deactivate()
{
  this->active = false;
}

void Request::cancel()
{
  this->cancelled = true;
}
