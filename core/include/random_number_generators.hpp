#ifndef __random_number_generators_hpp__
#define __random_number_generators_hpp__

#include <functional>

namespace decaygen {
  // a function that returns a random double between 0 and 1
  typedef std::function<double()> random_uniform_0_1_generator;
}

#endif
