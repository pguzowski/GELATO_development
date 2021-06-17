#ifndef __GELATO_physics_detail_integration_hpp__
#define __GELATO_physics_detail_integration_hpp__

#include <gsl/gsl_integration.h>

namespace GELATO {
  namespace physics {
    namespace detail {
      namespace integration {
        // allows passing lambda into gsl_integrate
        template<typename F>
        class gsl_function_wrapper : public gsl_function {
          public:
            gsl_function_wrapper(const F& func) : _func(func) {
              function = &gsl_function_wrapper::invoke;
              params=this;
            }
          private:
            const F& _func;
            static double invoke(double x, void *params) {
              return static_cast<gsl_function_wrapper*>(params)->_func(x);
            }
        };
        
        template<typename F>
        double integrate(F fn, double low, double high, double tolerance = 1e-5) {
          double result, error;
          size_t neval;
          gsl_function_wrapper<decltype(fn)> wrapped_fn(fn);
          gsl_function *func = static_cast<gsl_function*>(&wrapped_fn);   
          gsl_integration_qng(func, low, high, tolerance, tolerance, &result, &error, &neval);
          return result;
        }

        using integrate_alt_workspace = gsl_integration_workspace;
        template<typename F>
        double integrate_alt(F fn, double low, double high,/*
            integrate_alt_workspace* ws, size_t limit,*/ double tolerance = 1e-5){ 
          double result, error;
          gsl_function_wrapper<decltype(fn)> wrapped_fn(fn);
          gsl_function *func = static_cast<gsl_function*>(&wrapped_fn);
          size_t limit = 10000;
          auto ws = gsl_integration_workspace_alloc(limit);
          gsl_integration_qags(func, low, high, tolerance, tolerance, limit,  ws, &result, &error);
          gsl_integration_workspace_free(ws);
          return result;
        }
      }
    }
  }
}

#endif // __GELATO_physics_detail_integration_hpp__

