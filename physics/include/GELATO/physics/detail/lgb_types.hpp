#ifndef __GELATO_physics_detail_lgb_types_hpp__
#define __GELATO_physics_detail_lgb_types_hpp__

#include "GELATO/physics/detail/integration.hpp"

#include "GELATO/core/config.hpp"
#include <cmath>
#include <complex>

#include <iostream>

namespace GELATO {
  namespace physics {
    namespace detail {
      namespace lgb {

        constexpr int gauge_boson_pdg = 32;
          
        constexpr bool self_conjugate = true;
        constexpr bool final_state = true;
        constexpr bool NOT_final_state = !final_state;

        enum class decay_modes { nu_nu, e_e, mu_mu, tau_tau };
        enum class production_modes { k_mu, k_e, pi_mu, pi_e, pi0, eta };
        enum class leptophilic_gauges { Le_minus_Ltau, Le_minus_Lmu, Lmu_minus_Ltau };


        struct model_parameters {
          double boson_mass;
          double gauge_coupling;
          leptophilic_gauges gauge;

          // have to specify these because of the other constructor
          model_parameters() = default;
          ~model_parameters() = default;
          model_parameters(const model_parameters&) = default;
          model_parameters(model_parameters&&) = default;
          model_parameters& operator=(const model_parameters&) = default;
          model_parameters& operator=(model_parameters&&) = default;

          // copy from GELATO::physics::leptophilic_gauge_boson::model_parameters
          // without defining that struct in this header file or includes
          template<typename T> model_parameters(const T& p)
            : boson_mass{p.boson_mass}, gauge_coupling{p.gauge_coupling}, gauge{p.gauge} { }
        };
        // so we don't have to constantly recalculate these
        struct derived_params {
          const double coup2;
          double kinetic_mixing;
          const double mB2;
          double mla2;
          double mlb2;
          const double elec_charge;
          const double gFermi2;
          const double elec_mass;
          const double muon_mass;
          const double tau_mass;
          const double pion_pm_mass;
          const double kaon_pm_mass;
          const double pi3_inv;
          const model_parameters& raw_params;
          const GELATO::core::config& config;
          derived_params(const model_parameters& p,
              const GELATO::core::config& c) :
            coup2{std::pow(p.gauge_coupling,2)},
            mB2{std::pow(p.boson_mass,2)},
            elec_charge{c.physical_params().elec_charge},
            gFermi2{std::pow(c.physical_params().gFermi,2)},
            elec_mass{c.physical_params().find_particle("elec").mass},
            muon_mass{c.physical_params().find_particle("muon").mass},
            tau_mass{c.physical_params().find_particle("tau").mass},
            pion_pm_mass{c.physical_params().find_particle("pion_pm").mass},
            kaon_pm_mass{c.physical_params().find_particle("kaon_pm").mass},
            pi3_inv{1./std::pow(M_PI,3)},
            raw_params(p), config(c) {
              switch(p.gauge) {
                case leptophilic_gauges::Le_minus_Ltau:
                  mla2 = std::pow(elec_mass,2);
                  mlb2 = std::pow(tau_mass,2);
                  break;
                case leptophilic_gauges::Le_minus_Lmu:
                  mla2 = std::pow(elec_mass,2);
                  mlb2 = std::pow(muon_mass,2);
                  break;
                case leptophilic_gauges::Lmu_minus_Ltau:
                  mla2 = std::pow(muon_mass,2);
                  mlb2 = std::pow(tau_mass,2);
                default:
                  break;
              }
              //std::cerr << mla2 << " "<< mlb2 << " "<<mB2<< std::endl;
              auto integrand_re = [this](double x) {
                using complex = std::complex<double>;
                const double xbit = x*(1.-x);
                const double abit = mla2 - xbit * mB2;
                const double bbit = mlb2 - xbit * mB2;
                const double i = (abit <= 0. && bbit <= 0.) ? xbit * (std::log(-bbit) - std::log(-abit)) :
                  (abit > 0. && bbit > 0.) ? xbit * (std::log(bbit) - std::log(abit)) :
                  std::real(xbit * (std::log(complex{bbit,0.}) - std::log(complex{abit,0.})));
                //std::cerr << x <<" re "<<i<<" "<<abit<<" "<<bbit<<'\n';
                return i;
              };
              
              auto integrand_im = [this](double x) {
                using complex = std::complex<double>;
                const double xbit = x*(1.-x);
                const double abit = mla2 - xbit * mB2;
                const double bbit = mlb2 - xbit * mB2;
                const double i = (abit <= 0. && bbit <= 0.) ? 0. :
                  (abit > 0. && bbit > 0.) ? 0. :
                  std::imag(xbit * (std::log(complex{bbit,0.}) - std::log(complex{abit,0.})));
                //std::cerr << x <<" im "<<i<<" "<<abit<<" "<<bbit<<'\n';
                return i;
              };
              
              const double int_re = integration::integrate_alt(integrand_re, 0., 1.);
              const double int_im = integration::integrate_alt(integrand_im, 0., 1.);

              kinetic_mixing = -elec_charge * p.gauge_coupling / 2 / M_PI / M_PI * std::sqrt(int_re*int_re + int_im*int_im);
              //std::cerr << "kinetic mixing: "<<kinetic_mixing<<'\n';
            }
        };
      }
    }
  }
}
#endif
