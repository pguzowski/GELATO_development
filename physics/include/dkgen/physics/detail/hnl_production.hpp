#ifndef __dkgen_physics_detail_hnl_production_hpp__
#define __dkgen_physics_detail_hnl_production_hpp__

#include "dkgen/physics/detail/integration.hpp"
#include "dkgen/physics/detail/utils.hpp"
#include "dkgen/physics/detail/hnl_types.hpp"

#include "dkgen/core/particle.hpp"
#include <cmath>
#include <numeric>

#undef DEBUG
#ifdef DEBUG
#include <iostream>
#endif

namespace dkgen {
  namespace physics {
    namespace detail {
      namespace hnl {

        template<typename T, typename U>
        inline bool production_mode_enabled(T pm, const U& production_modes_to_use) {
          return std::find(production_modes_to_use.begin(), production_modes_to_use.end(), pm) != production_modes_to_use.end();
        };

        double I_h1(double x, double y, double z, double f0, double lam_plus, double lam_0, int plus_minus) {
          namespace detint = detail::integration;
          auto s_integrand = [plus_minus,f0,lam_plus,lam_0,x,y,z](double s) {
            auto t_integrand = [plus_minus,s,f0,lam_plus,lam_0,x,y,z](double t) {
              const double u = 1. + x + y + z - s - t;
              const double sqrt_kl_uyz = (std::sqrt(u) < std::sqrt(y) + std::sqrt(z)) ? 0. : plus_minus * detail::utils::sqrtkl(u,y,z);
              const double A = 0.5 * (1. + y - t) * (1. + z - s - plus_minus * detail::utils::sqrtkl(1,z,s))
                - 0.5 * (u - y - z - sqrt_kl_uyz);
              const double B = 0.5 * (y + z) * (u - y - z) + 2 * x * y - 0.5 * (y - z) * sqrt_kl_uyz;
              const double C = z * (1 + y - t) + (y + 0.5 * sqrt_kl_uyz) * (1 + z - s);
              const double F = f0 * (1. + lam_plus * u / x);
              const double G = f0 * (1. + lam_plus * u / x - (lam_plus - lam_0) * (1. + 1. / x));
              /*
              if(std::isnan(A) || std::isnan(B) || std::isnan(C) || std::isnan(F) || std::isnan(G)) {
#ifdef DEBUG
                std::cout<<"s " << s<<" t "<<t<<" u " <<u <<" x "<<x<<" y "<<y<<" z "<<z
                  <<" f0 "<<f0<<" l+ "<<lam_plus<<" l0 "<<lam_0
                  <<" kl_uyz " <<sqrt_kl_uyz<< " kl_1zs "<<detail::utils::sqrtkl(1,z,s)<<" "<<std::endl;
                std::cout << std::sqrt(u) <<" "<< std::sqrt(y) + std::sqrt(z)<<std::endl;
                std::cout << (1. + s*s + t*t + 2.*s* (-1. + t - x) + 2.*x + x*x - 2.*t*(1. + x) - 4*y*z) << std::endl;
#endif
                throw std::runtime_error("nan in integration");
              }
              */
              return F*F*A + G*G*B - F*G*C;
            };
            const double midpt = x + z + (1. - s - z) * (s - y + x) / 2. / s;
            const double delta = detail::utils::sqrtkl(s,x,y) * detail::utils::sqrtkl(1,s,z) / 2. / s;
            /*
            if(std::isnan(midpt) || std::isnan(delta)) {
#ifdef DEBUG
              std::cout << "s "<<s<<" midpt "<<midpt<<" delta "<<delta<< " sqrtkl(s,y,z) "<<utils::sqrtkl(s,y,z)<<" sqrtkl(1,s,z) "<<utils::sqrtkl(1,s,z)<<std::endl;
#endif
              throw std::runtime_error("nan in integration limits");
            }
            */
            return detint::integrate(t_integrand, midpt - delta, midpt + delta);
          };
          return std::abs(detint::integrate(s_integrand, std::pow(std::sqrt(x)+std::sqrt(y),2), std::pow(1.-std::sqrt(z),2)));
        }

        double P_decay_scale_factor_to_leptons(const derived_params& dparams,
            flavour fl, double pseudoscalar_mass, double lep_mass, int plus_minus) {
          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double xi_N = std::pow(HNL_mass/pseudoscalar_mass,2);
          const double xi_l = std::pow(lep_mass/pseudoscalar_mass,2);
          const double sqrt_kl = detail::utils::sqrtkl(1.,xi_N,xi_l);
          const double sub = xi_N - xi_l;
          return dparams.U2(fl) * sqrt_kl * (xi_l + xi_N - std::pow(sub,2) + plus_minus * sub * sqrt_kl)
            / (2.* xi_l * std::pow(1. - xi_l,2));
        };

        double P_decay_rate_to_semilepton_Pprime(const derived_params& dparams, flavour fl,
            double pseudoscalar_mass, double lep_mass, double pprime_mass, double CKM_V2,
            double f0_P_Pprime, double lam_plus_P_Pprime, double lam_0_P_Pprime,
            int plus_minus) {
          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double xi_h = std::pow(pprime_mass/pseudoscalar_mass,2);
          const double xi_N = std::pow(HNL_mass/pseudoscalar_mass,2);
          const double xi_l = std::pow(lep_mass/pseudoscalar_mass,2);
          const double prefac = dparams.gFermi2 * dparams.m5 * dparams.pi3_inv / 128. * dparams.U2(fl) * CKM_V2;
          return prefac * I_h1(xi_h, xi_l, xi_N, f0_P_Pprime, lam_plus_P_Pprime, lam_0_P_Pprime, plus_minus);
        };

        using decay_width_map_t = std::map<hnl_helicities,std::map<production_modes, double>>;

        template<typename T>
          double sum_decay_widths(const T& map1) {
            return std::accumulate(map1.begin(), map1.end(), 0.,
                [](double total1, auto& map2) {
                return std::accumulate(map2.second.begin(), map2.second.end(), total1,
                    [](double total2, auto& pair){
#ifdef DEBUG
                    std::cout << "sum_decay_widths: "<<pair.second << std::endl;
#endif
                      return total2 + pair.second;
                    });
                }
                );
          }

        constexpr double f0_k_pi = 0.9706; // https://arxiv.org/abs/1902.08191
        constexpr double lam_plus_k_pi = 0.0309; // PDG
        constexpr double lam_0_k_pi = 0.0173; // PDG



        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////                      KAONS                            /////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////



        decay_width_map_t get_kaon_pm_decay_widths(
            const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {
          /*
          auto production_mode_enabled = [production_modes_to_use](auto pm) -> bool {
            return std::find(production_modes_to_use.begin(), production_modes_to_use.end(), pm)
              != production_modes_to_use.end();
          };
          */
          decay_width_map_t decay_rates;
          auto& K_decay_rates_poshel = decay_rates[hnl_helicities::pos];
          auto& K_decay_rates_neghel = decay_rates[hnl_helicities::neg];

          const double total_sm_decay_rate = dparams.config.physical_params().hbar
            / dparams.config.physical_params().find_particle("kaon_pm").lifetime;

          const double kaon_pm_m2_BR = 0.6356;
          const double kaon_pm_e2_BR = 0.0016;

          const double CKM_Vus2 = std::pow(dparams.config.physical_params().CKM_Vus,2);

          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double muon_mass = dparams.muon_mass;
          const double elec_mass = dparams.elec_mass;
          const double kaon_pm_mass = dparams.kaon_pm_mass;
          const double pion_0_mass = dparams.pion_0_mass;

          if(kaon_pm_mass > HNL_mass + muon_mass && dparams.U2(flavour::m) > 0.
              && production_mode_enabled(production_modes::k_mu2,production_modes_to_use)) {
            const double sm_decay_rate = kaon_pm_m2_BR * total_sm_decay_rate;
            const double sf_pos = P_decay_scale_factor_to_leptons(dparams, 
                flavour::m, kaon_pm_mass, muon_mass, helicity_plus);
            K_decay_rates_poshel[production_modes::k_mu2] = sm_decay_rate * sf_pos;
#ifdef DEBUG
          std::cout << "Add K->mu_N+ r="<<K_decay_rates_poshel[production_modes::k_mu2]<<std::endl;
#endif
            const double sf_neg = P_decay_scale_factor_to_leptons(dparams, 
                flavour::m, kaon_pm_mass, muon_mass, helicity_minus);
            K_decay_rates_neghel[production_modes::k_mu2] = sm_decay_rate * sf_neg;
#ifdef DEBUG
          std::cout << "Add K->mu_N- r="<<K_decay_rates_neghel[production_modes::k_mu2]<<std::endl;
#endif
          }
          if(kaon_pm_mass > HNL_mass + elec_mass && dparams.U2(flavour::e) > 0.
              && production_mode_enabled(production_modes::k_e2,production_modes_to_use)) {
            const double sm_decay_rate = kaon_pm_e2_BR * total_sm_decay_rate;
            const double sf_pos = P_decay_scale_factor_to_leptons(dparams, 
                flavour::e, kaon_pm_mass, elec_mass, helicity_plus);
            K_decay_rates_poshel[production_modes::k_e2] = sm_decay_rate * sf_pos;
#ifdef DEBUG
          std::cout << "Add K->e_N+ r="<<K_decay_rates_poshel[production_modes::k_e2]<<std::endl;
#endif
            const double sf_neg = P_decay_scale_factor_to_leptons(dparams, 
                flavour::e, kaon_pm_mass, elec_mass, helicity_minus);
            K_decay_rates_neghel[production_modes::k_e2] = sm_decay_rate * sf_neg;
#ifdef DEBUG
          std::cout << "Add K->e_N- r="<<K_decay_rates_neghel[production_modes::k_e2]<<std::endl;
#endif
          }
          if(kaon_pm_mass > HNL_mass + muon_mass + pion_0_mass && dparams.U2(flavour::m) > 0.
              && production_mode_enabled(production_modes::k_mu3,production_modes_to_use)) {
            const double decay_rate_pos = P_decay_rate_to_semilepton_Pprime(dparams,
                flavour::m, kaon_pm_mass, muon_mass, pion_0_mass,
                CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_plus);
            K_decay_rates_poshel[production_modes::k_mu3] = decay_rate_pos;
#ifdef DEBUG
          std::cout << "Add K->pi_mu_N+ r="<<K_decay_rates_poshel[production_modes::k_mu3]<<std::endl;
#endif
            const double decay_rate_neg = P_decay_rate_to_semilepton_Pprime(dparams,
                flavour::m, kaon_pm_mass, muon_mass, pion_0_mass,
                CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_minus);
            K_decay_rates_neghel[production_modes::k_mu3] = decay_rate_neg;
#ifdef DEBUG
          std::cout << "Add K->pi_mu_N- r="<<K_decay_rates_neghel[production_modes::k_mu3]<<std::endl;
#endif
          }
          if(kaon_pm_mass > HNL_mass + elec_mass + pion_0_mass && dparams.U2(flavour::e) > 0.
              && production_mode_enabled(production_modes::k_e3,production_modes_to_use)) {
            const double decay_rate_pos = P_decay_rate_to_semilepton_Pprime(dparams,
                flavour::e, kaon_pm_mass, elec_mass, pion_0_mass,
                CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_plus);
            K_decay_rates_poshel[production_modes::k_e3] = decay_rate_pos;
#ifdef DEBUG
          std::cout << "Add K->pi_e_N+ r="<<K_decay_rates_poshel[production_modes::k_e3]<<std::endl;
#endif
            const double decay_rate_neg = P_decay_rate_to_semilepton_Pprime(dparams,
                flavour::e, kaon_pm_mass, elec_mass, pion_0_mass,
                CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_minus);
            K_decay_rates_neghel[production_modes::k_e3] = decay_rate_neg;
#ifdef DEBUG
          std::cout << "Add K->pi_e_N- r="<<K_decay_rates_neghel[production_modes::k_e3]<<std::endl;
#endif
          }        
          return decay_rates;
        };

        // calculates branching ratio of kaon to all possible modes of HNL, so that
        // weights can be applied to final event for given model parameters
        double kaon_pm_hnl_branching_ratio(const model_parameters& params, const dkgen::core::config& conf,
           const std::vector<production_modes>& production_modes_to_use) {
          auto const& decrates = get_kaon_pm_decay_widths(derived_params{params, conf}, production_modes_to_use);
          const double kaon_hnl_decay_rate = sum_decay_widths(decrates);
          const double total_kaon_decay_rate = kaon_hnl_decay_rate +
            (conf.physical_params().hbar / conf.physical_params().find_particle("kaon_pm").lifetime);
          return kaon_hnl_decay_rate / total_kaon_decay_rate;
        }

        dkgen::core::particle_definition make_kaon_definition(const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {
          
          const double HNL_mass = dparams.raw_params.HNL_mass;

          const int kaon_pm_pdg  = dparams.config.physical_params().find_particle("kaon_pm").pdgcode;
          const int pion_0_pdg  = dparams.config.physical_params().find_particle("pion_0").pdgcode;
          const int elec_pdg    = dparams.config.physical_params().find_particle("elec").pdgcode;
          const int muon_pdg    = dparams.config.physical_params().find_particle("muon").pdgcode;

          const double kaon_pm_mass = dparams.kaon_pm_mass;
          const double pion_0_mass  = dparams.pion_0_mass;
          const double elec_mass    = dparams.elec_mass;
          const double muon_mass    = dparams.muon_mass;
          const double kaon_pm_lt   = dparams.config.physical_params().find_particle("kaon_pm").lifetime;
          
          dkgen::core::particle_definition charged_kaon{kaon_pm_pdg,kaon_pm_mass,kaon_pm_lt};
          
          auto kaon_pm_decay_widths = get_kaon_pm_decay_widths(dparams, production_modes_to_use);
          const double total_kaon_decay_rate = sum_decay_widths(kaon_pm_decay_widths);
          auto& K_decay_rates_poshel = kaon_pm_decay_widths[hnl_helicities::pos];
          auto& K_decay_rates_neghel = kaon_pm_decay_widths[hnl_helicities::neg];


          if(kaon_pm_mass > HNL_mass + muon_mass && dparams.U2(flavour::m) > 0.
              && production_mode_enabled(production_modes::k_mu2,production_modes_to_use)) {
#ifdef DEBUG
            std::cerr << "adding: "<< K_decay_rates_poshel[production_modes::k_mu2]<<"/"<<total_kaon_decay_rate<<std::endl;
#endif
            charged_kaon.add_decay({
                K_decay_rates_poshel[production_modes::k_mu2]/total_kaon_decay_rate,
                {{HNL_pdg_poshel,NOT_final_state},{-muon_pdg,final_state}}
                });
#ifdef DEBUG
            std::cerr << "adding: "<< K_decay_rates_neghel[production_modes::k_mu2]<<"/"<<total_kaon_decay_rate<<std::endl;
#endif
            charged_kaon.add_decay({
                K_decay_rates_neghel[production_modes::k_mu2]/total_kaon_decay_rate,
                {{HNL_pdg_neghel,NOT_final_state},{-muon_pdg,final_state}}
                });
          }
          if(kaon_pm_mass > HNL_mass + elec_mass && dparams.U2(flavour::e) > 0.
              && production_mode_enabled(production_modes::k_e2,production_modes_to_use)) {
#ifdef DEBUG
            std::cerr << "adding: "<< K_decay_rates_poshel[production_modes::k_e2]<<"/"<<total_kaon_decay_rate<<std::endl;
#endif
            charged_kaon.add_decay({
                K_decay_rates_poshel[production_modes::k_e2]/total_kaon_decay_rate,
                {{HNL_pdg_poshel,NOT_final_state},{-elec_pdg,final_state}}
                });
#ifdef DEBUG
            std::cerr << "adding: "<< K_decay_rates_neghel[production_modes::k_e2]<<"/"<<total_kaon_decay_rate<<std::endl;
#endif
            charged_kaon.add_decay({
                K_decay_rates_neghel[production_modes::k_e2]/total_kaon_decay_rate,
                {{HNL_pdg_neghel,NOT_final_state},{-elec_pdg,final_state}}
                });
          }
          if(kaon_pm_mass > HNL_mass + muon_mass + pion_0_mass && dparams.U2(flavour::m) > 0.
              && production_mode_enabled(production_modes::k_mu3,production_modes_to_use)) {
#ifdef DEBUG
            std::cerr << "adding: "<< K_decay_rates_poshel[production_modes::k_mu3]<<"/"<<total_kaon_decay_rate<<std::endl;
#endif
            if(K_decay_rates_poshel[production_modes::k_mu3] < 0 && std::abs(K_decay_rates_poshel[production_modes::k_mu3]/total_kaon_decay_rate) < 1e-5) {
              K_decay_rates_poshel[production_modes::k_mu3] = -K_decay_rates_poshel[production_modes::k_mu3];
            }
            charged_kaon.add_decay({
                K_decay_rates_poshel[production_modes::k_mu3]/total_kaon_decay_rate,
                {{HNL_pdg_poshel,NOT_final_state},{-muon_pdg,final_state},{pion_0_pdg,final_state}}
                });
#ifdef DEBUG
            std::cerr << "adding: "<< K_decay_rates_neghel[production_modes::k_mu3]<<"/"<<total_kaon_decay_rate<<std::endl;
#endif
            if(K_decay_rates_neghel[production_modes::k_mu3] < 0 && std::abs(K_decay_rates_neghel[production_modes::k_mu3]/total_kaon_decay_rate) < 1e-5) {
              K_decay_rates_neghel[production_modes::k_mu3] = -K_decay_rates_neghel[production_modes::k_mu3];
            }

            charged_kaon.add_decay({
                K_decay_rates_neghel[production_modes::k_mu3]/total_kaon_decay_rate,
                {{HNL_pdg_neghel,NOT_final_state},{-muon_pdg,final_state},{pion_0_pdg,final_state}}
                });
          }
          if(kaon_pm_mass > HNL_mass + elec_mass + pion_0_mass && dparams.U2(flavour::e)
              && production_mode_enabled(production_modes::k_e3,production_modes_to_use)) {
#ifdef DEBUG
            std::cerr << "adding: "<< K_decay_rates_poshel[production_modes::k_e3]<<"/"<<total_kaon_decay_rate<<std::endl;
#endif
            if(K_decay_rates_poshel[production_modes::k_e3] < 0 && std::abs(K_decay_rates_poshel[production_modes::k_e3]/total_kaon_decay_rate) < 1e-5) {
              K_decay_rates_poshel[production_modes::k_e3] = -K_decay_rates_poshel[production_modes::k_e3];
            }
            charged_kaon.add_decay({
                K_decay_rates_poshel[production_modes::k_e3]/total_kaon_decay_rate,
                {{HNL_pdg_poshel,NOT_final_state},{-elec_pdg,final_state},{pion_0_pdg,final_state}}
                });
#ifdef DEBUG
            std::cerr << "adding: "<< K_decay_rates_neghel[production_modes::k_e3]<<"/"<<total_kaon_decay_rate<<std::endl;
#endif
            if(K_decay_rates_neghel[production_modes::k_e3] < 0 && std::abs(K_decay_rates_neghel[production_modes::k_e3]/total_kaon_decay_rate) < 1e-5) {
              K_decay_rates_neghel[production_modes::k_e3] = -K_decay_rates_neghel[production_modes::k_e3];
            }
            charged_kaon.add_decay({
                K_decay_rates_neghel[production_modes::k_e3]/total_kaon_decay_rate,
                {{HNL_pdg_neghel,NOT_final_state},{-elec_pdg,final_state},{pion_0_pdg,final_state}}
                });
          }
          charged_kaon.finalise_decay_table();
          return charged_kaon;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////                      K^0LS                            /////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////


        decay_width_map_t get_kaon_0L_decay_widths(
            const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {

          /*
          auto production_mode_enabled = [production_modes_to_use](auto pm) -> bool {
            return std::find(production_modes_to_use.begin(), production_modes_to_use.end(), pm)
              != production_modes_to_use.end();
          };
          */

          decay_width_map_t decay_rates;
          auto& K0_decay_rates_poshel = decay_rates[hnl_helicities::pos];
          auto& K0_decay_rates_neghel = decay_rates[hnl_helicities::neg];

          const double CKM_Vus2 = std::pow(dparams.config.physical_params().CKM_Vus,2);

          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double muon_mass = dparams.muon_mass;
          const double elec_mass = dparams.elec_mass;
          const double pion_pm_mass = dparams.pion_pm_mass;
          const double kaon_0L_mass = dparams.kaon_0L_mass;

          if(kaon_0L_mass > HNL_mass + muon_mass + pion_pm_mass && dparams.U2(flavour::m) > 0.
              && production_mode_enabled(production_modes::k0_mu,production_modes_to_use)) {
            const double decay_rate_pos = P_decay_rate_to_semilepton_Pprime(dparams,
                flavour::m, kaon_0L_mass, muon_mass, pion_pm_mass,
                CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_plus);
            K0_decay_rates_poshel[production_modes::k0_mu] = decay_rate_pos;
#ifdef DEBUG
          std::cout << "Add K0->pi_mu_N+ r="<<K0_decay_rates_neghel[production_modes::k0_mu]<<std::endl;
#endif
            const double decay_rate_neg = P_decay_rate_to_semilepton_Pprime(dparams,
                flavour::m, kaon_0L_mass, muon_mass, pion_pm_mass,
                CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_minus);
            K0_decay_rates_neghel[production_modes::k0_mu] = decay_rate_neg;
#ifdef DEBUG
          std::cout << "Add K0->pi_mu_N- r="<<K0_decay_rates_neghel[production_modes::k0_mu]<<std::endl;
#endif
          }
          if(kaon_0L_mass > HNL_mass + elec_mass + pion_pm_mass && dparams.U2(flavour::e) > 0.
              && production_mode_enabled(production_modes::k0_e,production_modes_to_use)) {
            const double decay_rate_pos = P_decay_rate_to_semilepton_Pprime(dparams,
                flavour::e, kaon_0L_mass, elec_mass, pion_pm_mass,
                CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_plus);
            K0_decay_rates_poshel[production_modes::k0_e] = decay_rate_pos;
#ifdef DEBUG
          std::cout << "Add K0->pi_e_N+ r="<<K0_decay_rates_neghel[production_modes::k0_e]<<std::endl;
#endif
            const double decay_rate_neg = P_decay_rate_to_semilepton_Pprime(dparams,
                flavour::e, kaon_0L_mass, elec_mass, pion_pm_mass,
                CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_minus);
            K0_decay_rates_neghel[production_modes::k0_e] = decay_rate_neg;
#ifdef DEBUG
          std::cout << "Add K0->pi_e_N- r="<<K0_decay_rates_neghel[production_modes::k0_e]<<std::endl;
#endif
          }
          return decay_rates;
        }

        // calculates branching ratio of kaon_0L to all possible modes of HNL, so that
        // weights can be applied to final event for given model parameters
        double kaon_0L_hnl_branching_ratio(const model_parameters& params, const dkgen::core::config& conf,
           const std::vector<production_modes>& production_modes_to_use) {
          auto const& decrates = get_kaon_0L_decay_widths(derived_params{params, conf}, production_modes_to_use);
          const double kaon_hnl_decay_rate = sum_decay_widths(decrates);
          const double total_kaon_decay_rate = kaon_hnl_decay_rate +
            (conf.physical_params().hbar / conf.physical_params().find_particle("kaon_0L").lifetime);
          return kaon_hnl_decay_rate / total_kaon_decay_rate;
        }

        dkgen::core::particle_definition make_K0L_definition(const derived_params& dparams, 
            const std::vector<production_modes>& production_modes_to_use) {

          const double HNL_mass = dparams.raw_params.HNL_mass;

          const int kaon_0L_pdg  = dparams.config.physical_params().find_particle("kaon_0L").pdgcode;
          const int pion_pm_pdg = dparams.config.physical_params().find_particle("pion_pm").pdgcode;
          const int elec_pdg    = dparams.config.physical_params().find_particle("elec").pdgcode;
          const int muon_pdg    = dparams.config.physical_params().find_particle("muon").pdgcode;

          const double kaon_0L_mass = dparams.kaon_0L_mass;
          const double pion_pm_mass = dparams.pion_pm_mass;
          const double elec_mass    = dparams.elec_mass;
          const double muon_mass    = dparams.muon_mass;
          const double kaon_0L_lt   = dparams.config.physical_params().find_particle("kaon_0L").lifetime;
          
          dkgen::core::particle_definition neutral_kaon{kaon_0L_pdg,kaon_0L_mass,kaon_0L_lt};
          
          auto kaon_0L_decay_widths = get_kaon_0L_decay_widths(dparams, production_modes_to_use);
          const double total_k0_decay_rate = sum_decay_widths(kaon_0L_decay_widths);

          auto& K0_decay_rates_poshel = kaon_0L_decay_widths[hnl_helicities::pos];
          auto& K0_decay_rates_neghel = kaon_0L_decay_widths[hnl_helicities::neg];


          if(kaon_0L_mass > HNL_mass + muon_mass + pion_pm_mass && dparams.U2(flavour::m) > 0.
              && production_mode_enabled(production_modes::k0_mu,production_modes_to_use)) {
#ifdef DEBUG
         std::cerr << "adding: "<< K0_decay_rates_poshel[production_modes::k0_mu]<<"/"<<total_k0_decay_rate<<std::endl;
#endif

            neutral_kaon.add_decay({
                .5 * K0_decay_rates_poshel[production_modes::k0_mu]/total_k0_decay_rate,
                {{HNL_pdg_poshel,NOT_final_state},{-muon_pdg,final_state},{-pion_pm_pdg,final_state}}
                });
            neutral_kaon.add_decay({
                .5 * K0_decay_rates_poshel[production_modes::k0_mu]/total_k0_decay_rate,
                {{-HNL_pdg_poshel,NOT_final_state},{muon_pdg,final_state},{pion_pm_pdg,final_state}}
                });
#ifdef DEBUG
         std::cerr << "adding: "<< K0_decay_rates_neghel[production_modes::k0_mu]<<"/"<<total_k0_decay_rate<<std::endl;
#endif
            neutral_kaon.add_decay({
                .5 * K0_decay_rates_neghel[production_modes::k0_mu]/total_k0_decay_rate,
                {{HNL_pdg_neghel,NOT_final_state},{-muon_pdg,final_state},{-pion_pm_pdg,final_state}}
                });
            neutral_kaon.add_decay({
                .5 * K0_decay_rates_neghel[production_modes::k0_mu]/total_k0_decay_rate,
                {{-HNL_pdg_neghel,NOT_final_state},{muon_pdg,final_state},{pion_pm_pdg,final_state}}
                });
          }
          if(kaon_0L_mass > HNL_mass + elec_mass + pion_pm_mass && dparams.U2(flavour::e) > 0.
              && production_mode_enabled(production_modes::k0_e,production_modes_to_use)) {
            neutral_kaon.add_decay({
                .5 * K0_decay_rates_poshel[production_modes::k0_e]/total_k0_decay_rate,
                {{HNL_pdg_poshel,NOT_final_state},{-elec_pdg,final_state},{-pion_pm_pdg,final_state}}
                });
            neutral_kaon.add_decay({
                .5 * K0_decay_rates_poshel[production_modes::k0_e]/total_k0_decay_rate,
                {{-HNL_pdg_poshel,NOT_final_state},{elec_pdg,final_state},{pion_pm_pdg,final_state}}
                });
            neutral_kaon.add_decay({
                .5 * K0_decay_rates_neghel[production_modes::k0_e]/total_k0_decay_rate,
                {{HNL_pdg_neghel,NOT_final_state},{-elec_pdg,final_state},{-pion_pm_pdg,final_state}}
                });
            neutral_kaon.add_decay({
                .5 * K0_decay_rates_neghel[production_modes::k0_e]/total_k0_decay_rate,
                {{-HNL_pdg_neghel,NOT_final_state},{elec_pdg,final_state},{pion_pm_pdg,final_state}}
                });
          }
          neutral_kaon.finalise_decay_table();
          return neutral_kaon;
        }



        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////                      PIONS                            /////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////




        decay_width_map_t get_pion_decay_widths(
            const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {

          /*
          auto production_mode_enabled = [production_modes_to_use](auto pm) -> bool {
            return std::find(production_modes_to_use.begin(), production_modes_to_use.end(), pm)
              != production_modes_to_use.end();
          };
          */

          decay_width_map_t decay_rates;
          auto& pi_decay_rates_poshel = decay_rates[hnl_helicities::pos];
          auto& pi_decay_rates_neghel = decay_rates[hnl_helicities::neg];

          const double total_sm_decay_rate = dparams.config.physical_params().hbar
            / dparams.config.physical_params().find_particle("pion_pm").lifetime;

          const double pion_pm_m2_BR = 0.9998;
          const double pion_pm_e2_BR = 0.0001;

          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double muon_mass = dparams.muon_mass;
          const double elec_mass = dparams.elec_mass;
          const double pion_pm_mass = dparams.pion_pm_mass;

          if(pion_pm_mass > HNL_mass + muon_mass && dparams.U2(flavour::m) > 0.
              && production_mode_enabled(production_modes::pi_mu,production_modes_to_use)) {
            const double sm_decay_rate = pion_pm_m2_BR * total_sm_decay_rate;
            const double sf_pos = P_decay_scale_factor_to_leptons(dparams,
                flavour::m, pion_pm_mass, muon_mass, helicity_plus);
            pi_decay_rates_poshel[production_modes::pi_mu] = sf_pos * sm_decay_rate;
#ifdef DEBUG
          std::cout << "Add pi->mu_N+ r="<<pi_decay_rates_neghel[production_modes::pi_mu]<<std::endl;
#endif
            const double sf_neg = P_decay_scale_factor_to_leptons(dparams,
                flavour::m, pion_pm_mass, muon_mass, helicity_minus);
            pi_decay_rates_neghel[production_modes::pi_mu] = sf_neg * sm_decay_rate;
#ifdef DEBUG
          std::cout << "Add pi->mu_N- r="<<pi_decay_rates_neghel[production_modes::pi_mu]<<std::endl;
#endif
          }
          if(pion_pm_mass > HNL_mass + elec_mass && dparams.U2(flavour::e) > 0.
              && production_mode_enabled(production_modes::pi_e,production_modes_to_use)) {
            const double sm_decay_rate = pion_pm_e2_BR * total_sm_decay_rate;
            const double sf_pos = P_decay_scale_factor_to_leptons(dparams,
                flavour::e, pion_pm_mass, elec_mass, helicity_plus);
            pi_decay_rates_poshel[production_modes::pi_e] = sf_pos * sm_decay_rate;
#ifdef DEBUG
          std::cout << "Add pi->e_N+ r="<<pi_decay_rates_neghel[production_modes::pi_e]<<std::endl;
#endif
            const double sf_neg = P_decay_scale_factor_to_leptons(dparams,
                flavour::e, pion_pm_mass, elec_mass, helicity_minus);
            pi_decay_rates_neghel[production_modes::pi_e] = sf_neg * sm_decay_rate;
#ifdef DEBUG
          std::cout << "Add pi->e_N- r="<<pi_decay_rates_neghel[production_modes::pi_e]<<std::endl;
#endif
          }
          return decay_rates;
        }

        // calculates branching ratio of pion to all possible modes of HNL, so that
        // weights can be applied to final event for given model parameters
        double pion_hnl_branching_ratio(const model_parameters& params, const dkgen::core::config& conf,
           const std::vector<production_modes>& production_modes_to_use) {
          auto const& decrates = get_pion_decay_widths(derived_params{params, conf}, production_modes_to_use);
          const double pion_hnl_decay_rate = sum_decay_widths(decrates);
          const double total_pion_decay_rate = pion_hnl_decay_rate +
            (conf.physical_params().hbar / conf.physical_params().find_particle("pion_pm").lifetime);
          return pion_hnl_decay_rate / total_pion_decay_rate;
        }

        dkgen::core::particle_definition make_pion_definition(const derived_params& dparams,
          const std::vector<production_modes>& production_modes_to_use) {

          const double HNL_mass = dparams.raw_params.HNL_mass;

          const int pion_pm_pdg = dparams.config.physical_params().find_particle("pion_pm").pdgcode;
          const int elec_pdg    = dparams.config.physical_params().find_particle("elec").pdgcode;
          const int muon_pdg    = dparams.config.physical_params().find_particle("muon").pdgcode;

          const double pion_pm_mass = dparams.pion_pm_mass;
          const double elec_mass    = dparams.elec_mass;
          const double muon_mass    = dparams.muon_mass;
          const double pion_pm_lt      = dparams.config.physical_params().find_particle("pion_pm").lifetime;
          
          dkgen::core::particle_definition charged_pion{pion_pm_pdg,pion_pm_mass,pion_pm_lt};
          
          auto pion_pm_decay_widths = get_pion_decay_widths(dparams, production_modes_to_use);
          const double total_pion_decay_rate = sum_decay_widths(pion_pm_decay_widths);

          auto& pi_decay_rates_poshel = pion_pm_decay_widths[hnl_helicities::pos];
          auto& pi_decay_rates_neghel = pion_pm_decay_widths[hnl_helicities::neg];


          if(pion_pm_mass > HNL_mass + muon_mass && dparams.U2(flavour::m) > 0.
              && production_mode_enabled(production_modes::pi_mu,production_modes_to_use)) {
            charged_pion.add_decay({
                pi_decay_rates_poshel[production_modes::pi_mu]/total_pion_decay_rate,
                {{HNL_pdg_poshel,NOT_final_state},{-muon_pdg,final_state}}
                });
            charged_pion.add_decay({
                pi_decay_rates_neghel[production_modes::pi_mu]/total_pion_decay_rate,
                {{HNL_pdg_neghel,NOT_final_state},{-muon_pdg,final_state}}
                });
          }
          if(pion_pm_mass > HNL_mass + elec_mass && dparams.U2(flavour::e) > 0.
              && production_mode_enabled(production_modes::pi_e,production_modes_to_use)) {
            charged_pion.add_decay({
                pi_decay_rates_poshel[production_modes::pi_e]/total_pion_decay_rate,
                {{HNL_pdg_poshel,NOT_final_state},{-elec_pdg,final_state}}
                });
            charged_pion.add_decay({
                pi_decay_rates_neghel[production_modes::pi_e]/total_pion_decay_rate,
                {{HNL_pdg_neghel,NOT_final_state},{-elec_pdg,final_state}}
                });
          }
          charged_pion.finalise_decay_table();
          return charged_pion;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////                      MUONS                            /////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        double I_l(double x, double y, double z, int plus_minus, bool is_anti)  {
          std::function<double(double)> integrand;
          if(is_anti) {
            integrand = [x,y,z,plus_minus](double s) -> double {
              const double sqrt_kl_1sz = detail::utils::sqrtkl(1,s,z);
              return (1.+z-s-plus_minus*sqrt_kl_1sz)*(s-x-y)*detail::utils::sqrtkl(s,x,y)*sqrt_kl_1sz/s;
            };
          } else {
            integrand = [x,y,z,plus_minus](double s) -> double  {
              const double sqrt_kl_sxy = detail::utils::sqrtkl(s,x,y);
              return (1.+z-s)*(s-x-y-plus_minus*sqrt_kl_sxy)*sqrt_kl_sxy*detail::utils::sqrtkl(1,s,z)/s;
            };
          }
          return 12. * detail::integration::integrate(integrand, std::pow(std::sqrt(x)+std::sqrt(y),2), std::pow(1.-std::sqrt(z),2));
        };

        // for mu+ -> e+ + ... (so for antiparticles)
        // split this into UalphaN and UbetaN for later use
        double triple_lepton_decay_rate(const derived_params& dparams, double U2,
            double parent_mass, double child_mass, int plus_minus, bool anti) {
          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double xi_N = std::pow(HNL_mass/parent_mass,2);
          const double xi_l = std::pow(child_mass/parent_mass,2);
          const double prefac = dparams.gFermi2 * dparams.m5 * dparams.pi3_inv / 192. * U2;
          return prefac * (anti ? I_l(0.,xi_l,xi_N,plus_minus,anti) : I_l(xi_N,xi_l,0.,plus_minus,!anti));
        };

        using muon_decay_width_map_t = std::map<hnl_helicities,std::map<muon_production_modes, double>>;
        muon_decay_width_map_t get_muon_decay_widths(const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {

          /*
          auto production_mode_enabled = [production_modes_to_use](auto pm) -> bool {
            return std::find(production_modes_to_use.begin(), production_modes_to_use.end(), pm)
              != production_modes_to_use.end();
          };
          */

          // https://arxiv.org/pdf/1905.00284.pdf p16
          // The decay of a charged lepton (antilepton) of flavour α to a charged lepton (antilepton)
          // of flavour β can be proportional to either |UαN| or |UβN|, producing a heavy Dirac
          // neutrino (antineutrino) in the first case or an antineutrino (neutrino) in the second case
          //
          // Formula 4.4 is for antileptons, hence mu+ -> (antineutrino), e+ -> (neutrino)

          muon_decay_width_map_t muon_decay_rates;
          auto& muon_decay_rates_poshel = muon_decay_rates[hnl_helicities::pos];
          auto& muon_decay_rates_neghel = muon_decay_rates[hnl_helicities::neg];
          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double muon_mass = dparams.muon_mass;
          const double elec_mass = dparams.elec_mass;


          if(muon_mass > HNL_mass + elec_mass && dparams.U2(flavour::m) + dparams.U2(flavour::e) > 0.
              && production_mode_enabled(production_modes::mu_e,production_modes_to_use)) {
            const double decay_rate_pos_Nbar = triple_lepton_decay_rate(dparams, dparams.U2(flavour::m),
                muon_mass, elec_mass, helicity_plus, false);
            muon_decay_rates_poshel[muon_production_modes::mu_e_N] = decay_rate_pos_Nbar;
#ifdef DEBUG
            std::cout << "Add decay_rate_pos_Nbar r="<<decay_rate_pos_Nbar<<std::endl;
#endif
            const double decay_rate_pos_N = triple_lepton_decay_rate(dparams, dparams.U2(flavour::e),
                muon_mass, elec_mass, helicity_plus, true);
            muon_decay_rates_poshel[muon_production_modes::mu_e_Nbar] = decay_rate_pos_N;
#ifdef DEBUG
            std::cout << "Add decay_rate_pos_N r="<<decay_rate_pos_N<<std::endl;
#endif
            const double decay_rate_neg_Nbar = triple_lepton_decay_rate(dparams, dparams.U2(flavour::m),
                muon_mass, elec_mass, helicity_minus, false);
            muon_decay_rates_neghel[muon_production_modes::mu_e_Nbar] = decay_rate_neg_Nbar;
#ifdef DEBUG
            std::cout << "Add decay_rate_neg_Nbar r="<<decay_rate_neg_Nbar<<std::endl;
#endif
            const double decay_rate_neg_N = triple_lepton_decay_rate(dparams, dparams.U2(flavour::e),
                muon_mass, elec_mass, helicity_minus, true);
            muon_decay_rates_neghel[muon_production_modes::mu_e_N] = decay_rate_neg_N;
#ifdef DEBUG
            std::cout << "Add decay_rate_neg_N r="<<decay_rate_neg_N<<std::endl;
#endif
          }
          return muon_decay_rates;
        }


        // calculates branching ratio of pion to all possible modes of HNL, so that
        // weights can be applied to final event for given model parameters
        double muon_hnl_branching_ratio(const model_parameters& params, const dkgen::core::config& conf,
           const std::vector<production_modes>& production_modes_to_use) {
          auto const& decrates = get_muon_decay_widths(derived_params{params, conf}, production_modes_to_use);
          const double muon_hnl_decay_rate = sum_decay_widths(decrates);
          const double total_muon_decay_rate = muon_hnl_decay_rate +
            (conf.physical_params().hbar / conf.physical_params().find_particle("muon").lifetime);
          return muon_hnl_decay_rate / total_muon_decay_rate;
        }


        dkgen::core::particle_definition make_muon_definition(const derived_params& dparams,
            const std::vector<production_modes>& production_modes_to_use) {

          const double HNL_mass = dparams.raw_params.HNL_mass;

          const int elec_pdg    = dparams.config.physical_params().find_particle("elec").pdgcode;
          const int muon_pdg    = dparams.config.physical_params().find_particle("muon").pdgcode;
          const int nu_pdg      = dparams.config.physical_params().find_particle("nu_e").pdgcode;

          const double elec_mass    = dparams.elec_mass;
          const double muon_mass    = dparams.muon_mass;
          const double muon_lt      = dparams.config.physical_params().find_particle("muon").lifetime;
          
          dkgen::core::particle_definition muon_def{muon_pdg,muon_mass,muon_lt};

          auto muon_decay_widths = get_muon_decay_widths(dparams, production_modes_to_use);
          const double total_muon_decay_rate = sum_decay_widths(muon_decay_widths);

          auto& muon_decay_rates_poshel = muon_decay_widths[hnl_helicities::pos];
          auto& muon_decay_rates_neghel = muon_decay_widths[hnl_helicities::neg];


          // need to flip HNL helicities as above decay rates were for antimuons
          if(muon_mass > HNL_mass + elec_mass && dparams.U2(flavour::m) + dparams.U2(flavour::e) > 0.
              && production_mode_enabled(production_modes::mu_e,production_modes_to_use)) {
            // deal with -ve rates due to floating point rounding errors
            if(muon_decay_rates_neghel[muon_production_modes::mu_e_N] < 0
                && std::abs(muon_decay_rates_neghel[muon_production_modes::mu_e_N])
                / (total_muon_decay_rate - muon_decay_rates_neghel[muon_production_modes::mu_e_N]) < 1e-8) {
              muon_decay_rates_neghel[muon_production_modes::mu_e_N] = 0.;
            }
            if(muon_decay_rates_neghel[muon_production_modes::mu_e_N] > 0.) {
              muon_def.add_decay({
                  muon_decay_rates_neghel[muon_production_modes::mu_e_N]/total_muon_decay_rate,
                  {{HNL_pdg_poshel,NOT_final_state},{elec_pdg,final_state},{-nu_pdg,final_state}}
                  });
            }
            if(muon_decay_rates_neghel[muon_production_modes::mu_e_Nbar] < 0
                && std::abs(muon_decay_rates_neghel[muon_production_modes::mu_e_Nbar])
                / (total_muon_decay_rate - muon_decay_rates_neghel[muon_production_modes::mu_e_Nbar]) < 1e-8) {
              muon_decay_rates_neghel[muon_production_modes::mu_e_N] = 0.;
            }
            if(muon_decay_rates_neghel[muon_production_modes::mu_e_Nbar] > 0.) {
              muon_def.add_decay({
                  muon_decay_rates_neghel[muon_production_modes::mu_e_Nbar]/total_muon_decay_rate,
                  {{-HNL_pdg_poshel,NOT_final_state},{elec_pdg,final_state},{nu_pdg,final_state}}
                  });
            }
            if(muon_decay_rates_poshel[muon_production_modes::mu_e_N] < 0
                && std::abs(muon_decay_rates_poshel[muon_production_modes::mu_e_N])
                / (total_muon_decay_rate - muon_decay_rates_poshel[muon_production_modes::mu_e_N]) < 1e-8) {
              muon_decay_rates_neghel[muon_production_modes::mu_e_N] = 0.;
            }
            if(muon_decay_rates_poshel[muon_production_modes::mu_e_N] > 0.) {
              muon_def.add_decay({
                  muon_decay_rates_poshel[muon_production_modes::mu_e_N]/total_muon_decay_rate,
                  {{HNL_pdg_neghel,NOT_final_state},{elec_pdg,final_state},{-nu_pdg,final_state}}
                  });
            }
            if(muon_decay_rates_poshel[muon_production_modes::mu_e_Nbar] < 0
                && std::abs(muon_decay_rates_poshel[muon_production_modes::mu_e_Nbar])
                / (total_muon_decay_rate - muon_decay_rates_poshel[muon_production_modes::mu_e_Nbar]) < 1e-8) {
              muon_decay_rates_neghel[muon_production_modes::mu_e_N] = 0.;
            }
            if(muon_decay_rates_poshel[muon_production_modes::mu_e_Nbar] > 0.) {
              muon_def.add_decay({
                  muon_decay_rates_poshel[muon_production_modes::mu_e_Nbar]/total_muon_decay_rate,
                  {{-HNL_pdg_neghel,NOT_final_state},{elec_pdg,final_state},{nu_pdg,final_state}}
                  });
            }
          }
          muon_def.finalise_decay_table();
          return muon_def;
        }
      }
    }
  }
}

#undef DEBUG

#endif // __dkgen_physics_detail_hnl_production_hpp__
