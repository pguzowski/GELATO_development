#ifndef __GELATO_physics_detail_lgb_decays_hpp__
#define __GELATO_physics_detail_lgb_decays_hpp__


#include "GELATO/physics/detail/integration.hpp"
#include "GELATO/physics/detail/utils.hpp"
#include "GELATO/physics/detail/lgb_types.hpp"

#include "GELATO/core/particle.hpp"
#include <cmath>
#include <numeric>

#undef DEBUG
#ifdef DEBUG
#include <iostream>
#endif

namespace GELATO {
  namespace physics {
    namespace detail {
      namespace lgb {

        core::particle_definition make_boson_definition(const derived_params& dparams,
            const std::vector<decay_modes>& decay_modes_to_use) {

          auto& elec    = dparams.config.physical_params().find_particle("elec");
          auto& muon    = dparams.config.physical_params().find_particle("muon");
          auto& tau     = dparams.config.physical_params().find_particle("tau");
          const int elec_pdg    = elec.pdgcode;
          const int muon_pdg    = muon.pdgcode;
          const int tau_pdg     = tau.pdgcode;
          const double elec_mass    = elec.mass;
          const double muon_mass    = muon.mass;
          const double tau_mass     = tau.mass;

          const double mass_boson = dparams.raw_params.boson_mass;

          const double hbar = dparams.config.physical_params().hbar;

          const double decay_width_to_nu_nu = std::pow(dparams.raw_params.gauge_coupling,2)
            * dparams.raw_params.boson_mass / 12. / M_PI;

          auto decay_width_to_leptons = [](const double coupling2, const double mass_lep2, const double mass_boson) {
            const double mass_boson2 = std::pow(mass_boson, 2);
            return coupling2 * mass_boson / 12. / M_PI * (1. + 2 * mass_lep2 / mass_boson2)
              * std::sqrt(1. - 4. * mass_lep2 / mass_boson2);
          };

          auto decay_mode_enabled = [&decay_modes_to_use](const decay_modes& mode) {
            return std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), mode) != decay_modes_to_use.end();
          };

          const double kinmix2 = std::pow(dparams.kinetic_mixing * dparams.elec_charge,2);

          const double decay_width_to_e_e = (mass_boson <= 2. * elec_mass) ? 0.
            : decay_width_to_leptons(
                dparams.raw_params.gauge == leptophilic_gauges::Lmu_minus_Ltau ? kinmix2 : std::pow(dparams.raw_params.gauge_coupling, 2),
                std::pow(elec_mass, 2), dparams.raw_params.boson_mass
                );

          const double decay_width_to_mu_mu = (mass_boson <= 2. * muon_mass) ? 0.
            : decay_width_to_leptons(
                dparams.raw_params.gauge == leptophilic_gauges::Le_minus_Ltau ? kinmix2 : std::pow(dparams.raw_params.gauge_coupling, 2),
                std::pow(muon_mass, 2), dparams.raw_params.boson_mass
                );

          const double decay_width_to_tau_tau = (mass_boson <= 2. * tau_mass) ? 0.
            : decay_width_to_leptons(
                dparams.raw_params.gauge == leptophilic_gauges::Lmu_minus_Ltau ? kinmix2 : std::pow(dparams.raw_params.gauge_coupling, 2),
                std::pow(tau_mass, 2), dparams.raw_params.boson_mass
                );


          const double total_decay_width = decay_width_to_nu_nu + decay_width_to_e_e + decay_width_to_mu_mu + decay_width_to_tau_tau;
          const double boson_lifetime = hbar / total_decay_width;


          GELATO::core::particle_definition boson_info{gauge_boson_pdg, dparams.raw_params.boson_mass, boson_lifetime, self_conjugate};

          if(mass_boson > 2. * elec_mass && decay_mode_enabled(decay_modes::e_e)) {
            boson_info.add_decay({decay_width_to_e_e/total_decay_width, {{elec_pdg,true},{-elec_pdg,true}}});
          }
          if(mass_boson > 2. * muon_mass && decay_mode_enabled(decay_modes::mu_mu)) {
            boson_info.add_decay({decay_width_to_mu_mu/total_decay_width, {{muon_pdg,true},{-muon_pdg,true}}});
          }
          if(mass_boson > 2. * tau_mass && decay_mode_enabled(decay_modes::tau_tau)) {
            boson_info.add_decay({decay_width_to_tau_tau/total_decay_width, {{tau_pdg,true},{-tau_pdg,true}}});
          }

          boson_info.finalise_decay_table();

          return boson_info;

        }

      }
    }
  }
}
#undef DEBUG

#endif // __GELATO_physics_detail_lgb_hpp__
