#ifndef __GELATO_core_ordered_list_of_vectors_hpp__
#define __GELATO_core_ordered_list_of_vectors_hpp__

#include "GELATO/core/vectors.hpp"

#include <vector>

namespace GELATO {
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
            for(auto& kv : list_of_vectors) {
              transform(kv.second);
            }
            return *this;
          }
      private:
        std::vector<std::pair<double,vector3>> list_of_vectors;
    };
  }
}

#endif
