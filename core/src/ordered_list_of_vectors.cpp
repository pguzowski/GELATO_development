#include "dkgen/core/ordered_list_of_vectors.hpp"

#include <utility>
#include <algorithm>

size_t dkgen::core::ordered_list_of_vectors::size() const {
  return map.size();
}

const dkgen::core::vector3& dkgen::core::ordered_list_of_vectors::first() const {
  // throws exception if empty
  if(map.empty()) {
    throw std::runtime_error("empty list of vectors");
  }
  return map.begin()->second;
}

const dkgen::core::vector3& dkgen::core::ordered_list_of_vectors::last() const {
  // throws exception if empty
  if(map.empty()) {
    throw std::runtime_error("empty list of vectors");
  }
  return map.rbegin()->second;
}

dkgen::core::ordered_list_of_vectors& dkgen::core::ordered_list_of_vectors::add(double pos, const dkgen::core::vector3& vec) {
  map.emplace(pos, vec);
  return *this;
}

dkgen::core::ordered_list_of_vectors& dkgen::core::ordered_list_of_vectors::add(double pos, dkgen::core::vector3&& vec) {
  map.emplace(pos,std::move(vec));
  return *this;
}

