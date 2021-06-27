#ifndef __GELATO_physics_leptophilic_gauge_boson_hpp__
#define __GELATO_physics_leptophilic_gauge_boson_hpp__

#include "GELATO/physics/detail/lgb.hpp"

#include "GELATO/core/driver.hpp"
#include "GELATO/core/particle.hpp"

/*
 *
 * API for the end user to use. All implementation in detail:: namespace
 *
 */


namespace GELATO {
  namespace physics {
    namespace leptophilic_gauge_boson {

      using decay_modes = detail::lgb::decay_modes; // { nu_nu, e_e, mu_mu, tau_tau }
      using production_modes = detail::lgb::production_modes; // { k_mu2, k_mu3, k_e2, k_e3, k0_mu, k0_e, pi_mu, pi_e, mu_e }
      using leptophilic_gauges = detail::lgb::leptophilic_gauges; // { Le_minus_Ltau, Le_minus_Lmu, Lmu_minus_Ltau }
      
      const std::vector<decay_modes> all_decay_modes{
        decay_modes::e_e,
        decay_modes::mu_mu,
        decay_modes::tau_tau,
      };
      const std::vector<production_modes> all_production_modes{
        production_modes::k_e,
        production_modes::k_mu,
        production_modes::pi_e,
        production_modes::pi_mu,
        production_modes::pi0,
        production_modes::eta,
      };
      const std::vector<production_modes> all_kaon_pm_production_modes{
        production_modes::k_mu, production_modes::k_e, };
      const std::vector<production_modes> all_pion_pm_production_modes{
        production_modes::pi_mu, production_modes::pi_e, };
      
      struct model_parameters {
        double boson_mass;
        double gauge_coupling;
        leptophilic_gauges gauge;
      };
      


      // Generate the particle (decay table) map for this model
      // model_parameters: the model parameters (mass, mixing angles)
      // conf: configuration of the main program
      // decay_modes_to_use: list of the HNL decay modes to simulate
      //    e.g. only N -> mu pi
      // production_modes_to_use: list of the HNL production modes to simulate
      //    e.g. only K -> mu N
      // ensure_observation_possible: whether to crash if no HNL production & decay possible
      //    given the model parameters that have been given
      //    (set to 'false' if you want to for example calculate branching ratios
      //       for various mass points in a loop)

      inline GELATO::core::driver::particle_map create_particle_content(
          const model_parameters& params,
          const GELATO::core::config& conf,
          // all "visible" final states
          const std::vector<decay_modes>& decay_modes_to_use = all_decay_modes,
          const std::vector<production_modes>& production_modes_to_use = all_production_modes,
          bool ensure_observation_possible = true) {

        return detail::lgb::create_particle_content(params, conf,
            decay_modes_to_use, production_modes_to_use, ensure_observation_possible);

      }
      
      
      
      

      // calculates total branching ratio of charged kaon to the requested production modes of boson, so that
      // weights can be applied to final event for given model parameters
      inline double kaon_pm_lgb_branching_ratio(const model_parameters& params, const GELATO::core::config& conf,
          const std::vector<production_modes>& production_modes_to_use) {
        return detail::lgb::kaon_pm_lgb_branching_ratio(params, conf, production_modes_to_use);
      }


      // calculates total branching ratio of charged pion to the requested production modes of boson, so that
      // weights can be applied to final event for given model parameters
      inline double pion_lgb_branching_ratio(const model_parameters& params, const GELATO::core::config& conf,
          const std::vector<production_modes>& production_modes_to_use) {
        return detail::lgb::pion_lgb_branching_ratio(params, conf, production_modes_to_use);
      }

      
      // calculates total branching ratio of neutral pion to the requested production modes of boson, so that
      // weights can be applied to final event for given model parameters
      inline double pi0_lgb_branching_ratio(const model_parameters& params, const GELATO::core::config& conf,
          const std::vector<production_modes>& production_modes_to_use) {
        return detail::lgb::pi0_lgb_branching_ratio(params, conf, production_modes_to_use);
      }
      
      
      // calculates total branching ratio of eta to the requested decay production of boson, so that
      // weights can be applied to final event for given model parameters
      inline double eta_lgb_branching_ratio(const model_parameters& params, const GELATO::core::config& conf,
          const std::vector<production_modes>& production_modes_to_use) {
        return detail::lgb::eta_lgb_branching_ratio(params, conf, production_modes_to_use);
      }
      
    }
  }
}

#endif
