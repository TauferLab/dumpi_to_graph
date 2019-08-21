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
