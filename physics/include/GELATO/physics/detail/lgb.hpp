#ifndef __GELATO_physics_detail_lgb_hpp__
#define __GELATO_physics_detail_lgb_hpp__

#include "GELATO/physics/detail/lgb_types.hpp"
#include "GELATO/physics/detail/lgb_decays.hpp"
#include "GELATO/physics/detail/lgb_production.hpp"

#include "GELATO/core/driver.hpp"

namespace GELATO {
  namespace physics {
    namespace detail {
      namespace lgb {

        GELATO::core::driver::particle_map create_particle_content(
            const model_parameters& params,
            const GELATO::core::config& conf,
            const std::vector<decay_modes>& decay_modes_to_use,
            const std::vector<production_modes>& production_modes_to_use,
            bool ensure_observation_possible) {

          if(decay_modes_to_use.empty() || production_modes_to_use.empty()) {
            throw std::runtime_error("No boson production or decay modes requested!");
          }

          GELATO::core::driver::particle_map ret;

          const derived_params dparams(params,conf);

          //////////////// KAON DECAYS ///////////////////////////////////////////////////////////////////////////////

          ret.push_back(make_kaon_definition(dparams,production_modes_to_use));
          const size_t n_kaon_decay_modes = ret.back().get_decay_table().size();

          //////////////// PION DECAYS ///////////////////////////////////////////////////////////////////////////////

          ret.push_back(make_pion_definition(dparams,production_modes_to_use));
          const size_t n_pion_decay_modes = ret.back().get_decay_table().size();

          //////////////// PI0 DECAYS  ///////////////////////////////////////////////////////////////////////////////

          
          ret.push_back(make_pi0_definition(dparams,production_modes_to_use));
          const size_t n_pi0_decay_modes = ret.back().get_decay_table().size();

          //////////////// ETA DECAYS  ///////////////////////////////////////////////////////////////////////////////

          ret.push_back(make_eta_definition(dparams,production_modes_to_use));
          const size_t n_eta_decay_modes = ret.back().get_decay_table().size();
          



          if(ensure_observation_possible &&
              n_kaon_decay_modes + n_pion_decay_modes + n_pi0_decay_modes + n_eta_decay_modes  == 0) {
            throw std::runtime_error("Unable to produce HNLs for these parameters and options");
          }

          //////////////// HNL DECAYS ///////////////////////////////////////////////////////////////////////////////

          ret.push_back(make_boson_definition(dparams, decay_modes_to_use));
          const size_t n_boson_decay_modes = ret.back().get_decay_table().size();
          
          if(ensure_observation_possible && n_boson_decay_modes == 0) {
            throw std::runtime_error("Unable to have decaying HNLs for these parameters and options");
          }
          
          //////////////// REMAINING PARTICLES ////////////////////////////////////////////////////////////////////

          auto& elec    = conf.physical_params().find_particle("elec");
          auto& muon    = conf.physical_params().find_particle("muon");
          auto& tau     = conf.physical_params().find_particle("tau");
          auto& nu      = conf.physical_params().find_particle("nu_e"); // all neutrinos will be saved as nu_e
          auto& gamma   = conf.physical_params().find_particle("gamma"); // all neutrinos will be saved as nu_e
          const int elec_pdg    = elec.pdgcode;
          const int muon_pdg    = muon.pdgcode;
          const int tau_pdg     = tau.pdgcode;
          const int nu_pdg   = nu.pdgcode;
          const int gamma_pdg = gamma.pdgcode;
          const double elec_mass    = elec.mass;
          const double muon_mass    = muon.mass;
          const double tau_mass     = tau.mass;
          const double elec_lt    = elec.lifetime;
          const double muon_lt    = muon.lifetime;
          const double tau_lt     = tau.lifetime;

          ret.push_back(GELATO::core::particle_definition{elec_pdg,elec_mass,elec_lt});
          ret.push_back(GELATO::core::particle_definition{muon_pdg,muon_mass,muon_lt});
          ret.push_back(GELATO::core::particle_definition{tau_pdg,tau_mass,tau_lt});
          ret.push_back(GELATO::core::particle_definition{nu_pdg,0.,-1.});
          ret.push_back(GELATO::core::particle_definition{gamma_pdg,0.,-1.});

          return ret;
        }
      }
    }
  }
}

#endif // __GELATO_physics_detail_lgb_hpp__
