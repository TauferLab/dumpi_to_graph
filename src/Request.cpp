#include "Request.hpp"

// Internal
#include "Channel.hpp"

Request::Request( int type, long id )
{
  this->type = type;
  this->id = id;
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
  if ( type == 0 || type == 1 ) {
    this->active = true;
  } else {
    this->active = false;
  }
  this->channel = channel;
}

int Request::get_type() const
{
  return this->type;
}

long Request::get_id() const
{
  return this->id;
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
