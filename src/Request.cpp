#include "Request.hpp"

// Internal
#include "Channel.hpp"

int Request::get_id() const
{
  return this->id;
}

const Channel& NonblockingCommunicationRequest::get_channel() const
{
  return this->channel;
}

bool NonblockingCommunicationRequest::is_cancelled() const
{
  return this->cancelled;
}

void NonblockingCommunicationRequest::cancel() 
{
  this->cancelled = true;
}

//#include <string>

//// Helper function for outputting request type
//std::string request_properties::RequestType_to_string(RequestType type)
//{
//  std::string s; 
//  switch( (int)type ) {
//    case 0: s = "irecv"; break;
//    case 1: s = "isend"; break;
//    case 2: s = "recv_init"; break;
//    case 3: s = "send_init"; break;
//  }
//  return s;
//}
//
//// Helper function for outputting request state
//std::string request_properties::RequestState_to_string(RequestState state)
//{
//  std::string s; 
//  switch( (int)state ) {
//    case 0: s = "active"; break;
//    case 1: s = "inactive"; break;
//    case 2: s = "cancelled"; break;
//  }
//  return s;
//}
