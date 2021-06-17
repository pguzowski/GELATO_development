#include "GELATO/core/ordered_list_of_vectors.hpp"

#include <utility>
#include <algorithm>

size_t GELATO::core::ordered_list_of_vectors::size() const {
  return map.size();
}

const GELATO::core::vector3& GELATO::core::ordered_list_of_vectors::first() const {
  // throws exception if empty
  if(map.empty()) {
    throw std::runtime_error("empty list of vectors");
  }
  auto f = std::min_element(map.begin(), map.end(), [](auto& a, auto& b){ return a.first < b.first; });
  return f->second;
}

const GELATO::core::vector3& GELATO::core::ordered_list_of_vectors::last() const {
  // throws exception if empty
  if(map.empty()) {
    throw std::runtime_error("empty list of vectors");
  }
  auto f = std::min_element(map.begin(), map.end(), [](auto& a, auto& b){ return a.first > b.first; });
  return f->second;
}

GELATO::core::ordered_list_of_vectors& GELATO::core::ordered_list_of_vectors::add(double pos, const vector3& vec) {
  map.emplace_back(pos, vec);
  return *this;
}

GELATO::core::ordered_list_of_vectors& GELATO::core::ordered_list_of_vectors::add(double pos, vector3&& vec) {
  map.emplace_back(pos,std::move(vec));
  return *this;
}

