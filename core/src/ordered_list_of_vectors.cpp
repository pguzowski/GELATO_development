#include "ordered_list_of_vectors.hpp"

#include <utility>
#include <algorithm>

size_t decaygen::ordered_list_of_vectors::size() const {
  return map.size();
}

const decaygen::vector3& decaygen::ordered_list_of_vectors::first() const {
  // throws exception if empty
  if(map.empty()) {
    throw std::runtime_error("empty list of vectors");
  }
  return map.begin()->second;
}

const decaygen::vector3& decaygen::ordered_list_of_vectors::last() const {
  // throws exception if empty
  if(map.empty()) {
    throw std::runtime_error("empty list of vectors");
  }
  return map.rbegin()->second;
}

decaygen::ordered_list_of_vectors& decaygen::ordered_list_of_vectors::add(double pos, const decaygen::vector3& vec) {
  map.emplace(pos, vec);
  return *this;
}

decaygen::ordered_list_of_vectors& decaygen::ordered_list_of_vectors::add(double pos, decaygen::vector3&& vec) {
  map.emplace(pos,std::move(vec));
  return *this;
}

