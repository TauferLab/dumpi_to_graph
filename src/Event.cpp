#include "Event.hpp"

#include <cinttypes>
#include <string>
#include <vector>

// Internal 
#include "Debug.hpp"

std::vector<std::string> Event::get_vertex_labels() const
{
  return this->vertex_labels;
}

size_t Event::get_vertex_id() const
{
  return this->vertex_id;
}

void Event::set_logical_time(size_t logical_time) 
{
  this->logical_time = logical_time;
}

size_t Event::get_logical_time() const 
{
  return this->logical_time;
}
  
int32_t Event::get_cpu_time_start() const
{
  return this->cpu_time.start.nsec;
}

int32_t Event::get_cpu_time_stop() const
{
  return this->cpu_time.stop.nsec;
}

int32_t Event::get_cpu_time_elapsed() const
{
  return (this->cpu_time.stop.nsec - this->cpu_time.start.nsec); 
}

double Event::get_cpu_time_midpoint() const
{
  double mid = (double)(this->cpu_time.stop.nsec - this->cpu_time.start.nsec)/2; 
  return (double)(this->cpu_time.start.nsec) + mid;
}

int32_t Event::get_wall_time_start() const
{
  return this->wall_time.start.nsec;
}

int32_t Event::get_wall_time_stop() const
{
  return this->wall_time.stop.nsec;
}

int32_t Event::get_wall_time_elapsed() const
{
  return (this->wall_time.stop.nsec - this->wall_time.start.nsec); 
}

double Event::get_wall_time_midpoint() const
{
  double mid = (double)(this->wall_time.stop.nsec - this->wall_time.start.nsec)/2; 
  return (double)(this->wall_time.start.nsec) + mid;
}

////////////////////////////////////////////////////////////////////////////////

size_t InitEvent::get_program_order_successor() const
{
  return this->program_order_successor;
}

void InitEvent::set_program_order_successor(size_t successor)
{
  this->program_order_successor = successor;
}

////////////////////////////////////////////////////////////////////////////////

size_t FinalizeEvent::get_program_order_predecessor() const
{
  return this->program_order_predecessor;
}

////////////////////////////////////////////////////////////////////////////////

size_t InteriorEvent::get_program_order_successor() const
{
  return this->program_order_successor;
}

void InteriorEvent::set_program_order_successor(size_t successor)
{
  this->program_order_successor = successor;
}

size_t InteriorEvent::get_program_order_predecessor() const
{
  return this->program_order_predecessor;
}

////////////////////////////////////////////////////////////////////////////////

void SendEvent::set_message_order_successor(size_t successor)
{
  this->message_order_successor = successor;
}

size_t SendEvent::get_message_order_successor() const
{
  return this->message_order_successor;
}

////////////////////////////////////////////////////////////////////////////////

void RecvEvent::set_message_order_predecessors(std::vector<size_t> predecessors)
{
  this->message_order_predecessors = predecessors; 
}

void RecvEvent::set_message_order_predecessor(size_t predecessor)
{
  this->message_order_predecessors.push_back( predecessor );
}

std::vector<size_t> RecvEvent::get_message_order_predecessors() const
{
#ifdef SAFETY
  assert( this->message_order_predecessors.size() > 0 );
#endif
  return this->message_order_predecessors;
}

size_t RecvEvent::get_message_order_predecessor() const
{
#ifdef SAFETY
  assert( this->message_order_predecessors.size() > 0 );
#endif
  return this->message_order_predecessors[0];
}


