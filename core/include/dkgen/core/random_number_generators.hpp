#ifndef __dkgen_core_random_number_generators_hpp__
#define __dkgen_core_random_number_generators_hpp__

#include <functional>

namespace dkgen {
  namespace core  {
    // a function that returns a random double between 0 and 1
    using random_uniform_0_1_generator = std::function<double()>;
  }
}

#endif
