#ifndef __GELATO_physics_detail_lgb_production_hpp__
#define __GELATO_physics_detail_lgb_production_hpp__

#include "GELATO/physics/detail/integration.hpp"
#include "GELATO/physics/detail/utils.hpp"
#include "GELATO/physics/detail/lgb_types.hpp"


#include "GELATO/core/config.hpp"
#include "GELATO/core/particle.hpp"
#include <cmath>
#include <numeric>

#undef DEBUG
#ifdef DEBUG
#include <iostream>
#endif

#include <iostream>
namespace GELATO {
  namespace physics {
    namespace detail {
      namespace lgb {

        template<typename T, typename U>
        inline bool production_mode_enabled(T pm, const U& production_modes_to_use) {
          return std::find(production_modes_to_use.begin(), production_modes_to_use.end(), pm) != production_modes_to_use.end();
        }

        inline double L3(double w, double x, double y) {
          const double sqrt_bit = w*w-2.*w*(1.+y)+std::pow(1.-y,2);
          const double S = (y-x)*std::sqrt(sqrt_bit);
          const double log_term = (x+y)*(1.-w)+x*y-y*y;
          //if(log_term + S < 0.) std::cerr << log_term <<" + S "<<S<<std::endl;
          //if(log_term - S < 0.) std::cerr << log_term <<" - S "<<S<<std::endl;
          // mistake in B.7  should be ... (1+^^^w^^^)xy -x(x+(1-w)2)-wy2 ... 2w(w+x-2)+2^^^y^^^(w-1)+y2+2 ...
          const double i = S * ((w+y-1.)*(x+3.*y)/(2.*y*y*y) + (1.-w)*(x+2.*w)/((1.+w)*x*y-x*(x+std::pow(1.-w,2))-w*y*y))
            + (2.*w*(w+x-2.)+2.*y*(w-1.)+y*y+2.)/y * (std::log(log_term + S)-std::log(log_term - S));
          //std::cerr << "w: "<<w<<" x: "<<x<<" y: "<<y<<" S: "<<S<<" sqrt "<<sqrt_bit<<" i: "<<i<<'\n';
          //std::cerr <<"num: "<< log_term + S << " S: "<<S<< " denom: "<<log_term -S<<'\n';
          return i;
        }

        double integral_L3(const double boson_mass2, const double meson_mass2, const double lep_mass2) {
          const double x = boson_mass2 / meson_mass2;
          const double w = lep_mass2 / meson_mass2;
          auto integrand_L3 = [x,w](double y) {
            return L3(w,x,y);
          };
          double ret = integration::integrate_alt(integrand_L3, x, std::pow(1.-std::sqrt(w),2));
          if(ret < 0. ) {
            std::cerr << " -ve integral "<<ret<<'\n';
          }
          return ret;
        }

        double partial_width_meson_to_lep_nu_vec(const derived_params& dparams, const double meson_mass2,
            const double VCKM2, const double lep_mass2, const double f_meson2) {
          const double meson_mass = std::sqrt(meson_mass2);
          return dparams.coup2 * dparams.gFermi2 * f_meson2 * VCKM2 * meson_mass * lep_mass2 * dparams.pi3_inv / 16.
            * integral_L3(dparams.mB2, meson_mass2, lep_mass2);
        }

        auto make_diff_rate_meson_to_lep_nu_vec(const derived_params& dparams, const double meson_mass2,
            const double lep_mass2) {
          const double x = dparams.mB2 / meson_mass2;
          const double w = lep_mass2 / meson_mass2;
          const double tot_integral = integral_L3(dparams.mB2, meson_mass2, lep_mass2);
          auto diff_rate_meson_to_lep_nu_vec = [w,x,tot_integral](double reduced_s12, double reduced_s13) {
            const double reduced_s23 = (1. + w + x - reduced_s12 - reduced_s13);
            const double y = reduced_s23;
            return L3(w,x,y) / tot_integral;
          };
          return diff_rate_meson_to_lep_nu_vec;
        }

        double neutral_meson_lgb_branching_ratio(const derived_params& dparams,
            double sm_BR, double meson_mass2, double boson_mass2) {
          if(boson_mass2 > meson_mass2) return 0.;
          const double kinmix2 = std::pow(dparams.kinetic_mixing,2);
          return sm_BR * 2. * kinmix2 * std::pow(1. - boson_mass2 / meson_mass2, 3);
        }


        template<typename T>
          double sum_decay_widths(const T& map1) {
            return std::accumulate(map1.begin(), map1.end(), 0.,
                [](double total, auto& pair){
                  return total + pair.second;
                });
          } 

        using decay_width_map_t = std::map<production_modes, double>;
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////                      KAONS                            /////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////



        decay_width_map_t get_kaon_pm_decay_widths(
            const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {

          decay_width_map_t decay_rates;

          const double CKM_Vus2 = std::pow(dparams.config.physical_params().CKM_Vus,2);
          const double fK2 = std::pow(dparams.config.physical_params().kaon_decay_constant,2);

          const double boson_mass = dparams.raw_params.boson_mass;
          const double muon_mass = dparams.muon_mass;
          const double elec_mass = dparams.elec_mass;
          const double kaon_pm_mass = dparams.kaon_pm_mass;


          if(kaon_pm_mass > boson_mass + muon_mass
              && dparams.raw_params.gauge != leptophilic_gauges::Le_minus_Ltau
              && production_mode_enabled(production_modes::k_mu,production_modes_to_use)) {
            decay_rates[production_modes::k_mu] = partial_width_meson_to_lep_nu_vec(dparams,
                std::pow(kaon_pm_mass,2), CKM_Vus2, std::pow(muon_mass,2), fK2);
          }
          if(kaon_pm_mass > boson_mass + elec_mass
              && dparams.raw_params.gauge != leptophilic_gauges::Lmu_minus_Ltau
              && production_mode_enabled(production_modes::k_e,production_modes_to_use)) {
            decay_rates[production_modes::k_e] = partial_width_meson_to_lep_nu_vec(dparams,
                std::pow(kaon_pm_mass,2), CKM_Vus2, std::pow(elec_mass,2), fK2);
          }
          return decay_rates;
        }

        // calculates branching ratio of kaon to all possible modes of HNL, so that
        // weights can be applied to final event for given model parameters
        double kaon_pm_lgb_branching_ratio(const model_parameters& params, const GELATO::core::config& conf,
           const std::vector<production_modes>& production_modes_to_use) {
          auto const& decrates = get_kaon_pm_decay_widths(derived_params{params, conf}, production_modes_to_use);
          const double kaon_lgb_decay_rate = sum_decay_widths(decrates);
          const double total_kaon_decay_rate = kaon_lgb_decay_rate +
            (conf.physical_params().hbar / conf.physical_params().find_particle("kaon_pm").lifetime);
          return kaon_lgb_decay_rate / total_kaon_decay_rate;
        }

        GELATO::core::particle_definition make_kaon_definition(const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {
          
          const double boson_mass = dparams.raw_params.boson_mass;

          const int kaon_pm_pdg  = dparams.config.physical_params().find_particle("kaon_pm").pdgcode;
          const int elec_pdg    = dparams.config.physical_params().find_particle("elec").pdgcode;
          const int muon_pdg    = dparams.config.physical_params().find_particle("muon").pdgcode;
          const int nu_pdg    = dparams.config.physical_params().find_particle("nu_e").pdgcode;

          const double kaon_pm_mass = dparams.kaon_pm_mass;
          const double elec_mass    = dparams.elec_mass;
          const double muon_mass    = dparams.muon_mass;
          const double kaon_pm_lt   = dparams.config.physical_params().find_particle("kaon_pm").lifetime;
          
          GELATO::core::particle_definition charged_kaon{kaon_pm_pdg,kaon_pm_mass,kaon_pm_lt};
          
          auto kaon_pm_decay_widths = get_kaon_pm_decay_widths(dparams, production_modes_to_use);
          const double total_kaon_decay_rate = sum_decay_widths(kaon_pm_decay_widths);


          if(kaon_pm_mass > boson_mass + muon_mass 
              && dparams.raw_params.gauge != leptophilic_gauges::Le_minus_Ltau
              && production_mode_enabled(production_modes::k_mu,production_modes_to_use)) {
            auto rw_func = make_diff_rate_meson_to_lep_nu_vec(dparams, std::pow(kaon_pm_mass,2), std::pow(muon_mass,2));
            charged_kaon.add_decay(core::decay_mode{
                kaon_pm_decay_widths[production_modes::k_mu]/total_kaon_decay_rate,
                {{-muon_pdg,final_state},{nu_pdg,final_state},{gauge_boson_pdg,NOT_final_state}}
                }.set_threebody_dalitz_reweighter(core::threebody_dalitz_function{rw_func}));
          }
          if(kaon_pm_mass > boson_mass + elec_mass
              && dparams.raw_params.gauge != leptophilic_gauges::Lmu_minus_Ltau
              && production_mode_enabled(production_modes::k_e,production_modes_to_use)) {
            auto rw_func = make_diff_rate_meson_to_lep_nu_vec(dparams, std::pow(kaon_pm_mass,2), std::pow(elec_mass,2));
            charged_kaon.add_decay(core::decay_mode{
                kaon_pm_decay_widths[production_modes::k_e]/total_kaon_decay_rate,
                {{-elec_pdg,final_state},{nu_pdg,final_state},{gauge_boson_pdg,NOT_final_state}}
                }.set_threebody_dalitz_reweighter(core::threebody_dalitz_function{rw_func}));
          }
          charged_kaon.finalise_decay_table();
          return charged_kaon;
        }




        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////                      PIONS                            /////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////




        decay_width_map_t get_pion_decay_widths(
            const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {

          decay_width_map_t decay_rates;

          const double CKM_Vud2 = std::pow(dparams.config.physical_params().CKM_Vud,2);
          const double fPi2 = std::pow(dparams.config.physical_params().pion_decay_constant,2);

          const double boson_mass = dparams.raw_params.boson_mass;
          const double muon_mass = dparams.muon_mass;
          const double elec_mass = dparams.elec_mass;
          const double pion_mass = dparams.pion_pm_mass;


          if(pion_mass > boson_mass + muon_mass
              && dparams.raw_params.gauge != leptophilic_gauges::Le_minus_Ltau
              && production_mode_enabled(production_modes::pi_mu,production_modes_to_use)) {
            decay_rates[production_modes::pi_mu] = partial_width_meson_to_lep_nu_vec(dparams,
                std::pow(pion_mass,2), CKM_Vud2, std::pow(muon_mass,2), fPi2);
          }
          if(pion_mass > boson_mass + elec_mass 
              && dparams.raw_params.gauge != leptophilic_gauges::Lmu_minus_Ltau
              && production_mode_enabled(production_modes::pi_e,production_modes_to_use)) {
            decay_rates[production_modes::pi_e] = partial_width_meson_to_lep_nu_vec(dparams,
                std::pow(pion_mass,2), CKM_Vud2, std::pow(elec_mass,2), fPi2);
          }
          return decay_rates;
        }

        // calculates branching ratio of pion to all possible modes of HNL, so that
        // weights can be applied to final event for given model parameters
        double pion_lgb_branching_ratio(const model_parameters& params, const GELATO::core::config& conf,
           const std::vector<production_modes>& production_modes_to_use) {
          auto const& decrates = get_pion_decay_widths(derived_params{params, conf}, production_modes_to_use);
          const double pion_lgb_decay_rate = sum_decay_widths(decrates);
          const double total_pion_decay_rate = pion_lgb_decay_rate +
            (conf.physical_params().hbar / conf.physical_params().find_particle("pion_pm").lifetime);
          return pion_lgb_decay_rate / total_pion_decay_rate;
        }

        GELATO::core::particle_definition make_pion_definition(const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {
          
          const double boson_mass = dparams.raw_params.boson_mass;

          const int pion_pdg  = dparams.config.physical_params().find_particle("pion_pm").pdgcode;
          const int elec_pdg    = dparams.config.physical_params().find_particle("elec").pdgcode;
          const int muon_pdg    = dparams.config.physical_params().find_particle("muon").pdgcode;
          const int nu_pdg    = dparams.config.physical_params().find_particle("nu_e").pdgcode;

          const double pion_mass = dparams.pion_pm_mass;
          const double elec_mass    = dparams.elec_mass;
          const double muon_mass    = dparams.muon_mass;
          const double pion_lt   = dparams.config.physical_params().find_particle("pion_pm").lifetime;
          
          GELATO::core::particle_definition charged_pion{pion_pdg,pion_mass,pion_lt};
          
          auto pion_decay_widths = get_pion_decay_widths(dparams, production_modes_to_use);
          const double total_pion_decay_rate = sum_decay_widths(pion_decay_widths);


          if(pion_mass > boson_mass + muon_mass 
              && production_mode_enabled(production_modes::pi_mu,production_modes_to_use)) {
            auto rw_func = make_diff_rate_meson_to_lep_nu_vec(dparams, std::pow(pion_mass,2), std::pow(muon_mass,2));
            charged_pion.add_decay(core::decay_mode{
                pion_decay_widths[production_modes::pi_mu]/total_pion_decay_rate,
                {{-muon_pdg,final_state},{nu_pdg,final_state},{gauge_boson_pdg,NOT_final_state}}
                }.set_threebody_dalitz_reweighter(core::threebody_dalitz_function{rw_func}));
          }
          if(pion_mass > boson_mass + elec_mass
              && production_mode_enabled(production_modes::pi_e,production_modes_to_use)) {
            auto rw_func = make_diff_rate_meson_to_lep_nu_vec(dparams, std::pow(pion_mass,2), std::pow(elec_mass,2));
            charged_pion.add_decay(core::decay_mode{
                pion_decay_widths[production_modes::pi_e]/total_pion_decay_rate,
                {{-elec_pdg,final_state},{nu_pdg,final_state},{gauge_boson_pdg,NOT_final_state}}
                }.set_threebody_dalitz_reweighter(core::threebody_dalitz_function{rw_func}));
          }
          charged_pion.finalise_decay_table();
          return charged_pion;
        }


        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////                      PI0S                             /////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        double pi0_lgb_branching_ratio(const model_parameters& params, const GELATO::core::config& conf,
           const std::vector<production_modes>& production_modes_to_use) {
          if(!production_mode_enabled(production_modes::pi0,production_modes_to_use)) return 0.;

          const double pi0_gamma_gamma_BR = 0.98823;
          const derived_params dparams{params,conf};
          const double pi0_mass2 = std::pow(conf.physical_params().find_particle("pion_0").mass,2);
          const double boson_mass2 = dparams.mB2;

          return neutral_meson_lgb_branching_ratio(dparams, pi0_gamma_gamma_BR, pi0_mass2, boson_mass2);
        }

        GELATO::core::particle_definition make_pi0_definition(const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {

          const double boson_mass = dparams.raw_params.boson_mass;

          auto const& pi0 = dparams.config.physical_params().find_particle("pion_0");
          const int pi0_pdg = pi0.pdgcode;
          const double pi0_mass = pi0.mass;
          const double pi0_lt = pi0.lifetime;

          const int gamma_pdg = dparams.config.physical_params().find_particle("gamma").pdgcode;


          GELATO::core::particle_definition pi0_defn{pi0_pdg,pi0_mass,pi0_lt};

          if(pi0_mass > boson_mass
              && production_mode_enabled(production_modes::pi0,production_modes_to_use)) {
            pi0_defn.add_decay({1.,{{gamma_pdg,final_state},{gauge_boson_pdg,NOT_final_state}}});
          }
          pi0_defn.finalise_decay_table();

          return pi0_defn;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////                      ETAS                             /////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////


        double eta_lgb_branching_ratio(const model_parameters& params, const GELATO::core::config& conf,
           const std::vector<production_modes>& production_modes_to_use) {
          if(!production_mode_enabled(production_modes::eta,production_modes_to_use)) return 0.;
          
          const double eta_gamma_gamma_BR = 0.3941;
          const derived_params dparams{params,conf};
          const double eta_mass2 = std::pow(conf.physical_params().find_particle("eta").mass,2);
          const double boson_mass2 = dparams.mB2;

          return neutral_meson_lgb_branching_ratio(dparams, eta_gamma_gamma_BR, eta_mass2, boson_mass2);
        }

        GELATO::core::particle_definition make_eta_definition(const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {

          const double boson_mass = dparams.raw_params.boson_mass;

          auto const& eta = dparams.config.physical_params().find_particle("eta");
          const int eta_pdg = eta.pdgcode;
          const double eta_mass = eta.mass;
          const double eta_lt = eta.lifetime;

          const int gamma_pdg = dparams.config.physical_params().find_particle("gamma").pdgcode;


          GELATO::core::particle_definition eta_defn{eta_pdg,eta_mass,eta_lt};

          if(eta_mass > boson_mass
              && production_mode_enabled(production_modes::eta,production_modes_to_use)) {
            eta_defn.add_decay({1.,{{gamma_pdg,final_state},{gauge_boson_pdg,NOT_final_state}}});
          }
          eta_defn.finalise_decay_table();

          return eta_defn;
        }


      }
    }
  }
}

#undef DEBUG

#endif // __GELATO_physics_detail_lgb_production_hpp__
