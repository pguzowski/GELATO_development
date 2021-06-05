#include "driver.hpp"
#include "particle.hpp"

#include <cmath>

namespace decaygen {
  namespace higgs_portal_from_kaons {
    decaygen::driver::particle_map create_particle_content(double scalar_mass, double theta,
        const decaygen::config& conf,
        decaygen::driver::particle_map input = {}) {

      auto& pion_pm = conf.physical_params().find_particle("pion_pm");
      auto& pion_0  = conf.physical_params().find_particle("pion_0");
      auto& kaon_pm = conf.physical_params().find_particle("kaon_pm");
      auto& kaon_0L = conf.physical_params().find_particle("kaon_0L");
      auto& elec    = conf.physical_params().find_particle("elec");
      auto& muon    = conf.physical_params().find_particle("muon");

      const int scalar_pdg  = 54; // free pdg code for bsm particles
      const int pion_pm_pdg = pion_pm.pdgcode;
      const int pion_0_pdg  = pion_0.pdgcode;
      const int kaon_pm_pdg = kaon_pm.pdgcode;
      const int kaon_0L_pdg = kaon_0L.pdgcode;
      const int elec_pdg    = elec.pdgcode;
      const int muon_pdg    = muon.pdgcode;

      // GeV system of units
      const double pion_pm_mass = pion_pm.mass;
      const double pion_0_mass  = pion_0.mass;
      const double kaon_pm_mass = kaon_pm.mass;
      const double kaon_0L_mass = kaon_0L.mass;
      const double elec_mass    = elec.mass;
      const double muon_mass    = muon.mass;

      const double hbar = conf.physical_params().hbar;
      const double higgs_vev = conf.physical_params().higgs_vev;

      using decay_type = enum {lepton, pion};
      auto decay_rate = [&scalar_mass, &theta, &higgs_vev](double daughter_mass, decay_type daughter_type ) {
        if(daughter_type == lepton) {
          return 
            theta*theta * daughter_mass*daughter_mass*scalar_mass
            / (8.*M_PI*higgs_vev*higgs_vev)
            * std::pow(1.-4.*daughter_mass*daughter_mass/scalar_mass/scalar_mass,1.5);
        }
        else {
          const double g = 2.*scalar_mass*scalar_mass/9. + 11.*daughter_mass*daughter_mass/9.;
          return
            theta*theta * 3.*g*g
            / (32.*M_PI*higgs_vev*higgs_vev*scalar_mass)
            * std::sqrt(1.-4.*daughter_mass*daughter_mass/scalar_mass/scalar_mass);
        }
      };

      double total_decay_rate = 0.;
      if(scalar_pdg > 2*elec_mass) {
        total_decay_rate += decay_rate(elec_mass, lepton);
      }
      if(scalar_pdg > 2*muon_mass) {
        total_decay_rate += decay_rate(muon_mass, lepton);
      }
      if(scalar_pdg > 2*pion_0_mass) {
        total_decay_rate += decay_rate(pion_0_mass, pion);
      }
      if(scalar_pdg > 2*pion_pm_mass) {
        total_decay_rate += decay_rate(pion_pm_mass, pion);
      }
      const double lifetime = hbar / total_decay_rate;

      // these will always be parent mesons, so no need to worry about lifetime
      input[kaon_pm_pdg] = decaygen::particle{kaon_pm_pdg,kaon_pm_mass,0.}.add_decay({1.,{{pion_pm_pdg,true},{scalar_pdg,false}}});
      input[kaon_0L_pdg] = decaygen::particle{kaon_0L_pdg,kaon_0L_mass,0.}.add_decay({1.,{{pion_0_pdg,true},{scalar_pdg,false}}});

      input[scalar_pdg] = decaygen::particle{scalar_pdg, scalar_mass, lifetime};
      
      if(scalar_mass > 2*elec_mass) {
        input[scalar_pdg].add_decay({decay_rate(elec_mass, lepton)/total_decay_rate,{{elec_pdg,true},{-elec_pdg,true}}});
      }
      if(scalar_mass > 2*muon_mass) {
        input[scalar_pdg].add_decay({decay_rate(muon_mass, lepton)/total_decay_rate,{{muon_pdg,true},{-muon_pdg,true}}});
      }
      if(scalar_mass > 2*pion_0_mass) {
        input[scalar_pdg].add_decay({decay_rate(pion_0_mass, pion)/total_decay_rate,{{pion_0_pdg,true},{pion_0_pdg,true}}});
      }
      if(scalar_mass > 2*pion_pm_mass) {
        input[scalar_pdg].add_decay({decay_rate(pion_pm_mass, pion)/total_decay_rate,{{pion_pm_mass,true},{-pion_pm_mass,true}}});
      }

      return input;
    }

    double kaon_branching_ratio(double scalar_mass, double theta, int kaon_pdg,
        const decaygen::config& conf) {

      auto& kaon_pm = conf.physical_params().find_particle("kaon_pm");
      auto& kaon_0L = conf.physical_params().find_particle("kaon_0L");
      
      const int kaon_pm_pdg = kaon_pm.pdgcode;
      const int kaon_0L_pdg = kaon_0L.pdgcode;
      
      if(kaon_pdg != kaon_0L_pdg && std::abs(kaon_pdg) != kaon_pm_pdg) {
        throw std::runtime_error("Kaon branching ratio requested for non-kaon pdg "+std::to_string(kaon_pdg));
      }

      const double kaon_pm_mass = kaon_pm.mass;
      const double kaon_0L_mass = kaon_0L.mass;
      const double kaon_pm_lifetime = kaon_pm.lifetime;
      const double kaon_0L_lifetime = kaon_0L.lifetime;
      const double pion_pm_mass = conf.physical_params().find_particle("pion_pm").mass;
      const double pion_0_mass = conf.physical_params().find_particle("pion_0").mass;
      
      const double hbar = conf.physical_params().hbar;
      const double higgs_vev = conf.physical_params().higgs_vev;
      const double Vtd = conf.physical_params().CKM_Vtd;
      const double Vts = conf.physical_params().CKM_Vts;
      const double top_mass = conf.physical_params().find_particle("top").mass;
      
      const double kmass = (kaon_pdg == kaon_0L_pdg ? kaon_0L_mass : kaon_pm_mass);
      const double k_lifetime = (kaon_pdg == kaon_0L_pdg ? kaon_0L_lifetime : kaon_pm_lifetime);
      const double pimass = (kaon_pdg == kaon_0L_pdg ? pion_0_mass : pion_pm_mass);
      
      // equation (5) of arXiv:1909.11670v1

      auto kallen_lambda = [](double a, double b, double c)->double {
        return a*a + b*b + c*c -2*a*b -2*a*c -2*b*c;
      };

      return
        theta * theta / (16*M_PI*kmass)
        * std::pow(3 * Vtd * Vts * top_mass*top_mass * kmass*kmass / (32*M_PI*M_PI * higgs_vev*higgs_vev*higgs_vev), 2) 
        * std::sqrt(kallen_lambda(1., scalar_mass*scalar_mass/kmass/kmass, pimass*pimass/kmass/kmass))
        / (hbar / k_lifetime);
    }
  }
}
