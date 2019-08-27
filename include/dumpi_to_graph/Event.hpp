#ifndef D2G_EVENT_H
#define D2G_EVENT_H

#include <cinttypes>
#include <string>
#include <vector>
#include <unordered_map>

// Internal
#include "Channel.hpp"
#include "Request.hpp"

class Event
{
public:
  Event(size_t vertex_id, 
        const dumpi_time cpu_time, 
        const dumpi_time wall_time) :
    vertex_id(vertex_id), cpu_time(cpu_time), wall_time(wall_time) 
    {}
  virtual ~Event() {}
  size_t get_vertex_id() const;
  size_t get_logical_time() const;
  int32_t get_cpu_time_start() const;
  int32_t get_cpu_time_stop() const;
  int32_t get_cpu_time_elapsed() const;
  double get_cpu_time_midpoint() const;
  int32_t get_wall_time_start() const;
  int32_t get_wall_time_stop() const;
  int32_t get_wall_time_elapsed() const;
  double get_wall_time_midpoint() const;
  std::vector<std::string> get_vertex_labels() const;
  void set_logical_time(const size_t logical_time); 
private:
  size_t vertex_id;
  // FIXME: Logical timestamp should be generalized to:
  // - Scalar values of arbitrary width (e.g., via GMP)
  // - Vector and matrix valued timestamps
  size_t logical_time;
  const dumpi_time cpu_time;
  const dumpi_time wall_time;
  std::vector<std::string> vertex_labels; 
};

class InitEvent : public Event
{
public:
  InitEvent(size_t vertex_id, 
            const dumpi_time cpu_time, 
            const dumpi_time wall_time) :
    Event(vertex_id, cpu_time, wall_time) 
    {}
  void set_program_order_successor(size_t successor);
  size_t get_program_order_successor() const; 
private:
  size_t program_order_successor;
};

class FinalizeEvent : public Event
{
public:
  FinalizeEvent(size_t vertex_id,
                const dumpi_time cpu_time, 
                const dumpi_time wall_time,
                const size_t predecessor) :
    Event(vertex_id, cpu_time, wall_time), 
    program_order_predecessor(predecessor) 
    {}
  size_t get_program_order_predecessor() const;
private:
  size_t program_order_predecessor;
};

class InteriorEvent : public Event
{
public:
  InteriorEvent(size_t vertex_id,
                const dumpi_time cpu_time,
                const dumpi_time wall_time,
                const size_t predecessor) :
    Event(vertex_id, cpu_time, wall_time),
    program_order_predecessor(predecessor) 
    {}
  virtual ~InteriorEvent() {}
  size_t get_program_order_predecessor() const;
  size_t get_program_order_successor() const; 
  void set_program_order_successor(size_t successor);
private:
  size_t program_order_predecessor;
  size_t program_order_successor;
};

class SendEvent : public InteriorEvent
{
public:
  SendEvent(size_t vertex_id,
            const dumpi_time cpu_time,
            const dumpi_time wall_time,
            const size_t predecessor) :
    InteriorEvent(vertex_id, cpu_time, wall_time, predecessor)
    {}
  void set_message_order_successor(size_t successor);
  size_t get_message_order_successor() const;
private:
  size_t message_order_successor;
  Channel channel;
  Request request;
};

class RecvEvent : public InteriorEvent
{
public:
  RecvEvent(size_t vertex_id,
            const dumpi_time cpu_time,
            const dumpi_time wall_time,
            const size_t predecessor) :
    InteriorEvent(vertex_id, cpu_time, wall_time, predecessor)
    {}
  void set_message_order_predecessors(std::vector<size_t> predecessors);
  void set_message_order_predecessor(size_t predecessor);
  std::vector<size_t> get_message_order_predecessors() const;
  size_t get_message_order_predecessor() const;
private:
  std::vector<size_t> message_order_predecessors;
  Channel channel;
  Request request;
};

//class CollectiveEvent : public InteriorEvent
//{
//public: 
//  CollectiveEvent(const size_t vertex_id,
//                  const dumpi_time cpu_time,
//                  const dumpi_time wall_time,
//                  const size_t predecessor) :
//    InteriorEvent(vertex_id, cpu_time, wall_time, predecessor)
//    {}
//  void add_peer(size_t peer_rank, size_t peer_vertex_id);
//  std::unordered_map get_peers() const;
//private:
//  std::unordered_map<size_t,size_t> peers;
//}


#endif // D2G_EVENT_H
