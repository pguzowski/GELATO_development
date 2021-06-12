#ifndef __dkgen_physics_heavy_neutral_leptons_hpp__
#define __dkgen_physics_heavy_neutral_leptons_hpp__

#include "dkgen/physics/detail/hnl.hpp"

#include "dkgen/core/driver.hpp"
#include "dkgen/core/particle.hpp"

/*
 *
 * API for the end user to use. All implementation in detail:: namespace
 *
 */


namespace dkgen {
  namespace physics {
    namespace heavy_neutral_leptons {

      using decay_modes = detail::hnl::decay_modes; // { nu_nu_nu, e_e_nu, e_mu_nu, mu_e_nu, mu_mu_nu, e_pi, mu_pi, pi0_nu }
      using production_modes = detail::hnl::production_modes; // { k_mu2, k_mu3, k_e2, k_e3, k0_mu, k0_e, pi_mu, pi_e, mu_e }
      
      const std::vector<decay_modes> all_decay_modes{
        decay_modes::nu_nu_nu,
        decay_modes::e_e_nu,
        decay_modes::e_mu_nu,
        decay_modes::mu_e_nu,
        decay_modes::mu_mu_nu,
        decay_modes::e_pi,
        decay_modes::mu_pi,
        decay_modes::pi0_nu
      };
      const std::vector<production_modes> all_production_modes{
        production_modes::k_mu2,
        production_modes::k_mu3,
        production_modes::k_e2,
        production_modes::k_e3,
        production_modes::k0_mu,
        production_modes::k0_e,
        production_modes::pi_mu,
        production_modes::pi_e,
        production_modes::mu_e
      };
      const std::vector<production_modes> all_kaon_pm_production_modes{
        production_modes::k_mu2, production_modes::k_mu3, production_modes::k_e2, production_modes::k_e3, };
      const std::vector<production_modes> all_pion_production_modes{
        production_modes::pi_mu, production_modes::pi_e, };
      const std::vector<production_modes> all_kaon_0L_production_modes{
        production_modes::k0_mu, production_modes::k0_e, };
      const std::vector<production_modes> all_muon_production_modes{
        production_modes::mu_e };
      
      struct model_parameters {
        double HNL_mass;
        double U_e4, U_m4, U_t4;
        bool is_majorana;
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

      inline dkgen::core::driver::particle_map create_particle_content(
          const model_parameters& params,
          const dkgen::core::config& conf,
          // all "visible" final states
          const std::vector<decay_modes>& decay_modes_to_use = all_decay_modes,
          const std::vector<production_modes>& production_modes_to_use = all_production_modes,
          bool ensure_observation_possible = true) {
        return detail::hnl::create_particle_content(params, conf, decay_modes_to_use,
            production_modes_to_use, ensure_observation_possible);
      }
      
      
      
      

      // calculates total branching ratio of neutral kaon to the requested production modes of HNL, so that
      // weights can be applied to final event for given model parameters
      inline double kaon_pm_hnl_branching_ratio(const model_parameters& params, const dkgen::core::config& conf,
          const std::vector<production_modes>& production_modes_to_use) {
        return detail::hnl::kaon_pm_hnl_branching_ratio(params, conf, production_modes_to_use);
      }

      // calculates total branching ratio of neutral kaon to the requested production modes of HNL, so that
      // weights can be applied to final event for given model parameters
      inline double kaon_0L_hnl_branching_ratio(const model_parameters& params, const dkgen::core::config& conf,
          const std::vector<production_modes>& production_modes_to_use) {
        return detail::hnl::kaon_0L_hnl_branching_ratio(params, conf, production_modes_to_use);
      }

      // calculates total branching ratio of charged pion to the requested production modes of HNL, so that
      // weights can be applied to final event for given model parameters
      inline double pion_hnl_branching_ratio(const model_parameters& params, const dkgen::core::config& conf,
          const std::vector<production_modes>& production_modes_to_use) {
        return detail::hnl::pion_hnl_branching_ratio(params, conf, production_modes_to_use);
      }

      // calculates total branching ratio of pion to the requested decay production of HNL, so that
      // weights can be applied to final event for given model parameters
      inline double muon_hnl_branching_ratio(const model_parameters& params, const dkgen::core::config& conf,
          const std::vector<production_modes>& production_modes_to_use) {
        return detail::hnl::muon_hnl_branching_ratio(params, conf, production_modes_to_use);
      }
    }
  }
}

#endif
