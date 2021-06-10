#ifndef __dkgen_physics_detail_utils_hpp__
#define __dkgen_physics_detail_utils_hpp__

namespace dkgen {
  namespace physics {
    namespace detail {
      namespace utils {
        inline double sqrt_kallen_lambda(double a, double b, double c) {
          return std::sqrt(a*a + b*b + c*c -2*a*b -2*a*c -2*b*c);
        };
        auto sqrtkl = sqrt_kallen_lambda; // alias with shorter name
      }
    }
  }
}

#endif //__dkgen_physics_detail_utils_hpp__
