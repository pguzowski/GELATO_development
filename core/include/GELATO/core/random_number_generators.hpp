#ifndef __GELATO_core_random_number_generators_hpp__
#define __GELATO_core_random_number_generators_hpp__

#include <functional>

namespace GELATO {
  namespace core  {
    // an end-user-defined function that returns a random double between 0 and 1
    using random_uniform_0_1_generator = std::function<double()>;
  }
}

#endif
