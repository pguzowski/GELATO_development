#ifndef __dkgen_core_ordered_list_of_vectors_hpp__
#define __dkgen_core_ordered_list_of_vectors_hpp__

#include "dkgen/core/vectors.hpp"

#include <vector>

namespace dkgen {
  namespace core {
    class ordered_list_of_vectors {
      public:
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
          }
      private:
        std::vector<std::pair<double,vector3>> map;
    };
  }
}

#endif
