#ifndef __ordered_list_of_vectors_hpp__
#define __ordered_list_of_vectors_hpp__

#include "vectors.hpp"
//#include <algorithm>
#include <map>

namespace decaygen {
  class ordered_list_of_vectors {
    public:
      ordered_list_of_vectors() = default;
      ~ordered_list_of_vectors() = default;
      ordered_list_of_vectors(const ordered_list_of_vectors&) = default;
      ordered_list_of_vectors(ordered_list_of_vectors&&) = default;
      ordered_list_of_vectors& operator=(const ordered_list_of_vectors&) = default;
      ordered_list_of_vectors& operator=(ordered_list_of_vectors&&) = default;
      const vector3& first() const;
      const vector3& last() const;
      size_t size() const;
      ordered_list_of_vectors& add(double position, const vector3& vector);
      ordered_list_of_vectors& add(double position, vector3&& vector);

      template<class UnaryFunction>
      ordered_list_of_vectors& apply_transformation(UnaryFunction transform) {
        for(auto& kv : map) {
          transform(kv.second);
        }
        return *this;
      };
    private:
      std::map<double,vector3> map;
  };
}

#endif
