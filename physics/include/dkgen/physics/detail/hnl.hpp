#ifndef __dkgen_physics_detail_hnl_hpp__
#define __dkgen_physics_detail_hnl_hpp__

#include "dkgen/physics/detail/hnl_types.hpp"
#include "dkgen/physics/detail/hnl_decays.hpp"
#include "dkgen/physics/detail/hnl_production.hpp"

namespace dkgen {
  namespace physics {
    namespace detail {
      namespace hnl {

        dkgen::core::driver::particle_map create_particle_content(
            const model_parameters& params,
            const dkgen::core::config& conf,
            const std::vector<decay_modes>& decay_modes_to_use,
            const std::vector<production_modes>& production_modes_to_use,
            bool ensure_observation_possible) {

          if(decay_modes_to_use.empty() || production_modes_to_use.empty()) {
            throw std::runtime_error("No HNL production or decay modes requested!");
          }

          dkgen::core::driver::particle_map ret;

          const derived_params dparams(params,conf);

          //////////////// KAON DECAYS ///////////////////////////////////////////////////////////////////////////////

          ret.push_back(make_kaon_definition(dparams,production_modes_to_use));
          const size_t n_kaon_decay_modes = ret.back().get_decay_table().size();

          //////////////// K0L DECAYS ///////////////////////////////////////////////////////////////////////////////

          ret.push_back(make_K0L_definition(dparams,production_modes_to_use));
          const size_t n_K0L_decay_modes = ret.back().get_decay_table().size();

          //////////////// PION DECAYS ///////////////////////////////////////////////////////////////////////////////

          ret.push_back(make_pion_definition(dparams,production_modes_to_use));
          const size_t n_pion_decay_modes = ret.back().get_decay_table().size();

          //////////////// MUON DECAYS ///////////////////////////////////////////////////////////////////////////////

          ret.push_back(make_muon_definition(dparams,production_modes_to_use));
          const size_t n_muon_decay_modes = ret.back().get_decay_table().size();



          if(ensure_observation_possible &&
              n_kaon_decay_modes + n_K0L_decay_modes + n_pion_decay_modes + n_muon_decay_modes == 0) {
            throw std::runtime_error("Unable to produce HNLs for these parameters and options");
          }

          //////////////// HNL DECAYS ///////////////////////////////////////////////////////////////////////////////

          auto const& hnls = make_hnl_definitions(dparams, decay_modes_to_use);
          
          if(ensure_observation_possible && hnls.first.get_decay_table().empty()
              && hnls.second.get_decay_table().empty()) {
            throw std::runtime_error("Unable to have decaying HNLs for these parameters and options");
          }
          
          ret.push_back(hnls.first);
          ret.push_back(hnls.second);
          
          //////////////// REMAINING PARTICLES ////////////////////////////////////////////////////////////////////

          auto& pion_0  = conf.physical_params().find_particle("pion_0");
          auto& elec    = conf.physical_params().find_particle("elec");
          auto& nu      = conf.physical_params().find_particle("nu_e"); // all neutrinos will be saved as nu_e
          const int pion_0_pdg  = pion_0.pdgcode;
          const int elec_pdg    = elec.pdgcode;
          const int nu_pdg   = nu.pdgcode;
          const double pion_0_mass  = pion_0.mass;
          const double elec_mass    = elec.mass;
          const double pion_0_lt  = pion_0.lifetime;
          const double elec_lt    = elec.lifetime;

          ret.push_back(dkgen::core::particle_definition{pion_0_pdg,pion_0_mass,pion_0_lt,self_conjugate});
          ret.push_back(dkgen::core::particle_definition{elec_pdg,elec_mass,elec_lt});
          ret.push_back(dkgen::core::particle_definition{nu_pdg,0.,-1.});

          return ret;
        }
      }
    }
  }
}

#endif // __dkgen_physics_detail_hnl_hpp__
