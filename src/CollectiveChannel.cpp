#include "CollectiveChannel.hpp"

int CollectiveChannel::get_root() const {return this->root;}

dumpi_comm CollectiveChannel::get_comm() const {return this->comm;}

int CollectiveChannel::get_type() const {return this->type;}

std::vector<int> CollectiveChannel::get_global_ranks() const {return this->global_ranks;}

void CollectiveChannel::set_root(int r){this->root = r;}

void CollectiveChannel::set_comm(dumpi_comm d){this->comm = d;}

void CollectiveChannel::set_type(int t){this->type = t;}

void CollectiveChannel::set_global_ranks(std::vector<int> v){this->global_ranks = v;}