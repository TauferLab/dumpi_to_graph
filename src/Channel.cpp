#include "Channel.hpp"

int Channel::get_src() const
{
  return this->src;
}

int Channel::get_dst() const
{
  return this->dst;
}

int Channel::get_tag() const
{
  return this->tag;
}

int Channel::get_comm() const
{
  return this->comm;
}

void Channel::set_src( int src )
{
  this->src = src;
}
void Channel::set_dst( int dst )
{
  this->dst = dst;
}
void Channel::set_tag( int tag )
{
  this->tag = tag;
}
void Channel::set_comm( int comm )
{
  this->comm = comm;
}
