#ifndef __dkgen_physics_heavy_neutral_leptons_hpp__
#define __dkgen_physics_heavy_neutral_leptons_hpp__

#include "dkgen/core/driver.hpp"
#include "dkgen/core/particle.hpp"

#include <cmath>

#include "dkgen/physics/detail/integration.hpp"
#include "dkgen/physics/detail/hnl.hpp"

//#define DEBUG
#undef DEBUG
#ifdef DEBUG
#include <iostream>
#endif

namespace dkgen {
  namespace physics {
    namespace heavy_neutral_leptons {

      enum class decay_modes { nu_nu_nu, e_e_nu, e_mu_nu, mu_e_nu, mu_mu_nu, e_pi, mu_pi, pi0_nu };
      enum class production_modes { k_mu2, k_mu3, k_e2, k_e3, k0_mu, k0_e, pi_mu, pi_e, mu_e };
      const std::vector<decay_modes> all_decay_modes{
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
      
      struct model_parameters {
        double HNL_mass;
        double U_e4, U_m4, U_t4;
        bool is_majorana;
      };
      
      dkgen::core::driver::particle_map create_particle_content(
          const model_parameters& params,
          const dkgen::core::config& conf,
          // all "visible" final states
          const std::vector<decay_modes>& decay_modes_to_use = all_decay_modes,
          const std::vector<production_modes>& production_modes_to_use = all_production_modes,
          dkgen::core::driver::particle_map input = {}) {

        if(decay_modes_to_use.empty() || production_modes_to_use.empty()) {
          throw std::runtime_error("No HNL production or decay modes requested!");
        }

        auto ret = input;

        auto& pion_pm = conf.physical_params().find_particle("pion_pm");
        auto& pion_0  = conf.physical_params().find_particle("pion_0");
        auto& kaon_pm = conf.physical_params().find_particle("kaon_pm");
        auto& kaon_0L = conf.physical_params().find_particle("kaon_0L");
        auto& elec    = conf.physical_params().find_particle("elec");
        auto& muon    = conf.physical_params().find_particle("muon");
        auto& nu    = conf.physical_params().find_particle("nu_e"); // all neutrinos will be saved as nu_e
        //auto& nu_mu   = conf.physical_params().find_particle("nu_mu");
        //auto& nu_tau  = conf.physical_params().find_particle("nu_tau");

        const int HNL_pdg_poshel  = 91; // free pdg code for 4th generation neutrino = 90 + helicity
        const int HNL_pdg_neghel  = 89; // free pdg code for 4th generation neutrino = 90 + helicity
                                        // antipartciles have opposite helicity
                                        // -91: neghel; -89: poshel (still works with -90 + helicity scheme)
        const int pion_pm_pdg = pion_pm.pdgcode;
        const int pion_0_pdg  = pion_0.pdgcode;
        const int kaon_pm_pdg = kaon_pm.pdgcode;
        const int kaon_0L_pdg = kaon_0L.pdgcode;
        const int elec_pdg    = elec.pdgcode;
        const int muon_pdg    = muon.pdgcode;
        const int nu_pdg   = nu.pdgcode;
        //const int nu_mu_pdg   = nu_mu.pdgcode;
        //const int nu_tau_pdg   = nu_tau.pdgcode;

        // GeV system of units
        const double HNL_mass = params.HNL_mass;
        const double pion_pm_mass = pion_pm.mass;
        const double pion_0_mass  = pion_0.mass;
        const double kaon_pm_mass = kaon_pm.mass;
        const double kaon_0L_mass = kaon_0L.mass;
        const double elec_mass    = elec.mass;
        const double muon_mass    = muon.mass;

        const double kaon_pm_lt = kaon_pm.lifetime;
        const double kaon_0L_lt = kaon_0L.lifetime;
        const double pion_pm_lt = pion_pm.lifetime;
        const double pion_0_lt  = pion_0.lifetime;
        const double elec_lt    = elec.lifetime;
        const double muon_lt    = muon.lifetime;

        const double hbar = conf.physical_params().hbar;
        const double gfermi2 = std::pow(conf.physical_params().gFermi,2);

        const double f2_pion = conf.physical_params().pion_decay_constant;
        const double CKM_Vud2 = std::pow(conf.physical_params().CKM_Vud,2);
        const double CKM_Vus2 = std::pow(conf.physical_params().CKM_Vus,2);
        
        
        const bool self_conjugate = true;
        const bool final_state = true;
        const bool NOT_final_state = !final_state;

        const int helicity_plus = +1;
        const int helicity_minus = -1;
        
        
        // only implementing total decay rates for now;
        constexpr double pi3_inv = 1./M_PI/M_PI/M_PI; // 1/pi^3

        enum class flavour { e, m, t };
        struct derived_params {
          double sum_U2;
          double m5;
          double m3;
          const model_parameters& raw_params;
          derived_params(const model_parameters& p) : raw_params(p) {
            sum_U2 = std::pow(p.U_e4,2) + std::pow(p.U_m4,2) + std::pow(p.U_t4,2);
            m5 = std::pow(p.HNL_mass,5);
            m3 = std::pow(p.HNL_mass,3);
          }
          double U2(flavour fl) const {
            switch(fl) {
              case flavour::e:
                return std::pow(raw_params.U_e4,2);
              case flavour::m:
                return std::pow(raw_params.U_m4,2);
              case flavour::t:
                return std::pow(raw_params.U_t4,2);
            }
          }
        };
        const derived_params dparams(params);

        auto sqrt_kallen_lambda = [](double a, double b, double c)->double {
          return std::sqrt(a*a + b*b + c*c -2*a*b -2*a*c -2*b*c);
        };
        auto sqrtkl = sqrt_kallen_lambda; // alias with shorter name
        
        namespace detint = detail::integration;
        auto I_h1 = [sqrtkl](double x, double y, double z, double f0, double lam_plus, double lam_0, int plus_minus) {
          auto s_integrand = [sqrtkl,plus_minus,f0,lam_plus,lam_0,x,y,z](double s) {
            auto t_integrand = [sqrtkl,plus_minus,s,f0,lam_plus,lam_0,x,y,z](double t) {
              const double u = 1. + x + y + z - s - t;
              if(std::sqrt(u) < std::sqrt(y) + std::sqrt(z)) {
                throw std::runtime_error("will get nan in integration");
              }
              const double sqrt_kl_uyz = (std::sqrt(u) < std::sqrt(y) + std::sqrt(z)) ? 0. : plus_minus * sqrtkl(u,y,z);
              const double A = 0.5 * (1. + y - t) * (1. + z - s - plus_minus * sqrtkl(1,z,s))
                - 0.5 * (u - y - z - sqrt_kl_uyz);
              const double B = 0.5 * (y + z) * (u - y - z) + 2 * x * y - 0.5 * (y - z) * sqrt_kl_uyz;
              const double C = z * (1 + y - t) + (y + 0.5 * sqrt_kl_uyz) * (1 + z - s);
              const double F = f0 * (1. + lam_plus * u / x);
              const double G = f0 * (1. + lam_plus * u / x - (lam_plus - lam_0) * (1. + 1. / x));
              if(std::isnan(A) || std::isnan(B) || std::isnan(C) || std::isnan(F) || std::isnan(G)) {
#ifdef DEBUG
                std::cout<<"s " << s<<" t "<<t<<" u " <<u <<" x "<<x<<" y "<<y<<" z "<<z
                  <<" f0 "<<f0<<" l+ "<<lam_plus<<" l0 "<<lam_0
                  <<" kl_uyz " <<sqrt_kl_uyz<< " kl_1zs "<<sqrtkl(1,z,s)<<" "<<std::endl;
                std::cout << std::sqrt(u) <<" "<< std::sqrt(y) + std::sqrt(z)<<std::endl;
                std::cout << (1. + s*s + t*t + 2.*s* (-1. + t - x) + 2.*x + x*x - 2.*t*(1. + x) - 4*y*z) << std::endl;
#endif
                throw std::runtime_error("nan in integration");
              }
              return F*F*A + G*G*B - F*G*C;
            };
            const double midpt = x + z + (1. - s - z) * (s - y + x) / 2. / s;
            const double delta = sqrtkl(s,x,y) * sqrtkl(1,s,z) / 2. / s;
            if(std::isnan(midpt) || std::isnan(delta)) {
#ifdef DEBUG
              std::cout << "s "<<s<<" midpt "<<midpt<<" delta "<<delta<< " sqrtkl(s,y,z) "<<sqrtkl(s,y,z)<<" sqrtkl(1,s,z) "<<sqrtkl(1,s,z)<<std::endl;
#endif
                throw std::runtime_error("nan in integration limits");
            }
            return detint::integrate(t_integrand, midpt - delta, midpt + delta);
          };
          return detint::integrate(s_integrand, std::pow(std::sqrt(x)+std::sqrt(y),2), std::pow(1.-std::sqrt(z),2));
        };

        auto P_decay_rate_to_leptons = [&dparams,sqrtkl,HNL_mass](flavour fl, double sm_nu_br, double pseudoscalar_mass, double lep_mass, int plus_minus) {
          const double xi_N = std::pow(HNL_mass/pseudoscalar_mass,2);
          const double xi_l = std::pow(lep_mass/pseudoscalar_mass,2);
          const double sqrt_kl = sqrtkl(1.,xi_N,xi_l);
          const double sub = xi_N - xi_l;
          return dparams.U2(fl) * sqrt_kl * (xi_l + xi_N - std::pow(sub,2) + plus_minus * sub * sqrt_kl)
            / (2.* xi_l * std::pow(1. - xi_l,2)) * sm_nu_br;
        };
        auto P_decay_rate_to_semilepton_Pprime = [&dparams,HNL_mass,I_h1](flavour fl, double sm_nu_br,
            double pseudoscalar_mass, double lep_mass, double pprime_mass, double CKM_V2,
            double f0_P_Pprime, double lam_plus_P_Pprime, double lam_0_P_Pprime,
            int plus_minus) {
          const double xi_h = std::pow(pprime_mass/pseudoscalar_mass,2);
          const double xi_N = std::pow(HNL_mass/pseudoscalar_mass,2);
          const double xi_l = std::pow(lep_mass/pseudoscalar_mass,2);
          return dparams.U2(fl) * CKM_V2 * sm_nu_br
            * I_h1(xi_h, xi_l, xi_N, f0_P_Pprime, lam_plus_P_Pprime, lam_0_P_Pprime, plus_minus)
            / I_h1(xi_h, xi_l, 0,    f0_P_Pprime, lam_plus_P_Pprime, lam_0_P_Pprime, plus_minus);
        };

        const double kaon_pm_m2_BR = 0.6356;
        const double kaon_pm_e2_BR = 0.0016;
        const double kaon_pm_m3_BR = 0.0335;
        const double kaon_pm_e3_BR = 0.0507;
        const double kaon_0L_m3_BR = 0.2704;
        const double kaon_0L_e3_BR = 0.4055;
        const double pion_pm_m2_BR = 0.9998;
        const double pion_pm_e2_BR = 0.0001;
        
        const double f0_k_pi = 0.9706; // https://arxiv.org/abs/1902.08191
        const double lam_plus_k_pi = 0.0309; // PDG
        const double lam_0_k_pi = 0.0173; // PDG

        auto production_mode_enabled = [production_modes_to_use](auto pm) -> bool {
          return std::find(production_modes_to_use.begin(), production_modes_to_use.end(), pm) != production_modes_to_use.end();
        };
        
        
        
        
        
        //////////////// KAON DECAYS ///////////////////////////////////////////////////////////////////////////////
        
        double total_kaon_decay_rate = 0.;
        std::map<production_modes,double> K_decay_rates_poshel, K_decay_rates_neghel;
        if(kaon_pm_mass > HNL_mass + muon_mass && params.U_m4 > 0.
            && production_mode_enabled(production_modes::k_mu2)) {
          const double decay_rate_pos = P_decay_rate_to_leptons(
              flavour::m, kaon_pm_m2_BR, kaon_pm_mass, muon_mass, helicity_plus);
          K_decay_rates_poshel[production_modes::k_mu2] = decay_rate_pos;
          const double decay_rate_neg = P_decay_rate_to_leptons(
              flavour::m, kaon_pm_m2_BR, kaon_pm_mass, muon_mass, helicity_minus);
          K_decay_rates_neghel[production_modes::k_mu2] = decay_rate_neg;
          total_kaon_decay_rate += decay_rate_pos + decay_rate_neg;
        }
        if(kaon_pm_mass > HNL_mass + elec_mass && params.U_e4 > 0.
            && production_mode_enabled(production_modes::k_e2)) {
          const double decay_rate_pos = P_decay_rate_to_leptons(
              flavour::e, kaon_pm_e2_BR, kaon_pm_mass, elec_mass, helicity_plus);
          K_decay_rates_poshel[production_modes::k_e2] = decay_rate_pos;
          const double decay_rate_neg = P_decay_rate_to_leptons(
              flavour::e, kaon_pm_e2_BR, kaon_pm_mass, elec_mass, helicity_minus);
          K_decay_rates_neghel[production_modes::k_e2] = decay_rate_neg;
          total_kaon_decay_rate += decay_rate_pos + decay_rate_neg;
        }
        if(kaon_pm_mass > HNL_mass + muon_mass + pion_0_mass && params.U_m4 > 0.
            && production_mode_enabled(production_modes::k_mu3)) {
          const double decay_rate_pos = P_decay_rate_to_semilepton_Pprime(
              flavour::m, kaon_pm_m3_BR, kaon_pm_mass, muon_mass, pion_0_mass,
              CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_plus);
          K_decay_rates_poshel[production_modes::k_mu3] = decay_rate_pos;
          const double decay_rate_neg = P_decay_rate_to_semilepton_Pprime(
              flavour::m, kaon_pm_m3_BR, kaon_pm_mass, muon_mass, pion_0_mass,
              CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_minus);
          K_decay_rates_neghel[production_modes::k_mu3] = decay_rate_neg;
          total_kaon_decay_rate += decay_rate_pos + decay_rate_neg;
        }
        if(kaon_pm_mass > HNL_mass + elec_mass + pion_0_mass && params.U_e4 > 0.
            && production_mode_enabled(production_modes::k_e3)) {
          const double decay_rate_pos = P_decay_rate_to_semilepton_Pprime(
              flavour::e, kaon_pm_e3_BR, kaon_pm_mass, elec_mass, pion_0_mass,
              CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_plus);
          K_decay_rates_poshel[production_modes::k_e3] = decay_rate_pos;
          const double decay_rate_neg = P_decay_rate_to_semilepton_Pprime(
              flavour::e, kaon_pm_e3_BR, kaon_pm_mass, elec_mass, pion_0_mass,
              CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_minus);
          K_decay_rates_neghel[production_modes::k_e3] = decay_rate_neg;
          total_kaon_decay_rate += decay_rate_pos + decay_rate_neg;
        }        
        
        
        ret.push_back(dkgen::core::particle_definition{kaon_pm_pdg,kaon_pm_mass,kaon_pm_lt});
        auto& charged_kaon = ret.back();
        
        if(kaon_pm_mass > HNL_mass + muon_mass && params.U_m4 > 0.
            && production_mode_enabled(production_modes::k_mu2)) {
          charged_kaon.add_decay({
              K_decay_rates_poshel[production_modes::k_mu2]/total_kaon_decay_rate,
              {{HNL_pdg_poshel,NOT_final_state},{-muon_pdg,final_state}}
          });
          charged_kaon.add_decay({
              K_decay_rates_neghel[production_modes::k_mu2]/total_kaon_decay_rate,
              {{HNL_pdg_neghel,NOT_final_state},{-muon_pdg,final_state}}
          });
        }
        if(kaon_pm_mass > HNL_mass + elec_mass && params.U_e4 > 0.
            && production_mode_enabled(production_modes::k_e2)) {
          charged_kaon.add_decay({
              K_decay_rates_poshel[production_modes::k_e2]/total_kaon_decay_rate,
              {{HNL_pdg_poshel,NOT_final_state},{-elec_pdg,final_state}}
          });
          charged_kaon.add_decay({
              K_decay_rates_neghel[production_modes::k_e2]/total_kaon_decay_rate,
              {{HNL_pdg_neghel,NOT_final_state},{-elec_pdg,final_state}}
          });
        }
        if(kaon_pm_mass > HNL_mass + muon_mass + pion_0_mass && params.U_m4 > 0.
            && production_mode_enabled(production_modes::k_mu3)) {
          charged_kaon.add_decay({
              K_decay_rates_poshel[production_modes::k_mu3]/total_kaon_decay_rate,
              {{HNL_pdg_poshel,NOT_final_state},{-muon_pdg,final_state},{pion_0_pdg,final_state}}
          });
          charged_kaon.add_decay({
              K_decay_rates_neghel[production_modes::k_mu3]/total_kaon_decay_rate,
              {{HNL_pdg_neghel,NOT_final_state},{-muon_pdg,final_state},{pion_0_pdg,final_state}}
          });
        }
        if(kaon_pm_mass > HNL_mass + elec_mass + pion_0_mass && params.U_e4 > 0.
            && production_mode_enabled(production_modes::k_e3)) {
          charged_kaon.add_decay({
              K_decay_rates_poshel[production_modes::k_e3]/total_kaon_decay_rate,
              {{HNL_pdg_poshel,NOT_final_state},{-elec_pdg,final_state},{pion_0_pdg,final_state}}
          });
          charged_kaon.add_decay({
              K_decay_rates_neghel[production_modes::k_e3]/total_kaon_decay_rate,
              {{HNL_pdg_neghel,NOT_final_state},{-elec_pdg,final_state},{pion_0_pdg,final_state}}
          });
        }
        charged_kaon.finalise_decay_table();

        
        //////////////// K0L DECAYS ///////////////////////////////////////////////////////////////////////////////
        
        
        
        double total_k0_decay_rate = 0.;
        std::map<production_modes,double> K0_decay_rates_poshel, K0_decay_rates_neghel;
        if(kaon_0L_mass > HNL_mass + muon_mass + pion_0_mass && params.U_m4 > 0.
            && production_mode_enabled(production_modes::k0_mu)) {
          const double decay_rate_pos = P_decay_rate_to_semilepton_Pprime(
              flavour::m, kaon_0L_m3_BR, kaon_0L_mass, muon_mass, pion_pm_mass,
              CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_plus);
          K0_decay_rates_poshel[production_modes::k0_mu] = decay_rate_pos;
          const double decay_rate_neg = P_decay_rate_to_semilepton_Pprime(
              flavour::m, kaon_0L_m3_BR, kaon_pm_mass, muon_mass, pion_0_mass,
              CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_minus);
          K0_decay_rates_neghel[production_modes::k0_mu] = decay_rate_neg;
          total_k0_decay_rate += decay_rate_pos + decay_rate_neg;
        }
        if(kaon_0L_mass > HNL_mass + elec_mass + pion_0_mass && params.U_e4 > 0.
            && production_mode_enabled(production_modes::k0_e)) {
          const double decay_rate_pos = P_decay_rate_to_semilepton_Pprime(
              flavour::e, kaon_0L_e3_BR, kaon_pm_mass, elec_mass, pion_0_mass,
              CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_plus);
          K0_decay_rates_poshel[production_modes::k0_e] = decay_rate_pos;
          const double decay_rate_neg = P_decay_rate_to_semilepton_Pprime(
              flavour::e, kaon_0L_e3_BR, kaon_pm_mass, elec_mass, pion_0_mass,
              CKM_Vus2, f0_k_pi, lam_plus_k_pi, lam_0_k_pi, helicity_minus);
          K0_decay_rates_neghel[production_modes::k0_e] = decay_rate_neg;
          total_k0_decay_rate += decay_rate_pos + decay_rate_neg;
        }        
        
        
        ret.push_back(dkgen::core::particle_definition{kaon_0L_pdg,kaon_0L_mass,kaon_0L_lt});
        auto& neutral_kaon = ret.back();
        
        if(kaon_pm_mass > HNL_mass + muon_mass + pion_0_mass && params.U_m4 > 0.
            && production_mode_enabled(production_modes::k_mu3)) {
          neutral_kaon.add_decay({
              .5 * K0_decay_rates_poshel[production_modes::k0_mu]/total_k0_decay_rate,
              {{HNL_pdg_poshel,NOT_final_state},{-muon_pdg,final_state},{pion_0_pdg,final_state}}
          });
          neutral_kaon.add_decay({
              .5 * K0_decay_rates_poshel[production_modes::k0_mu]/total_k0_decay_rate,
              {{-HNL_pdg_poshel,NOT_final_state},{muon_pdg,final_state},{pion_0_pdg,final_state}}
          });
          neutral_kaon.add_decay({
              .5 * K0_decay_rates_neghel[production_modes::k0_mu]/total_k0_decay_rate,
              {{HNL_pdg_neghel,NOT_final_state},{-muon_pdg,final_state},{pion_0_pdg,final_state}}
          });
          neutral_kaon.add_decay({
              .5 * K0_decay_rates_neghel[production_modes::k0_mu]/total_k0_decay_rate,
              {{-HNL_pdg_neghel,NOT_final_state},{muon_pdg,final_state},{pion_0_pdg,final_state}}
          });
        }
        if(kaon_pm_mass > HNL_mass + elec_mass + pion_0_mass && params.U_e4 > 0.
            && production_mode_enabled(production_modes::k_e3)) {
          neutral_kaon.add_decay({
              .5 * K0_decay_rates_poshel[production_modes::k0_e]/total_k0_decay_rate,
              {{HNL_pdg_poshel,NOT_final_state},{-elec_pdg,final_state},{pion_0_pdg,final_state}}
          });
          neutral_kaon.add_decay({
              .5 * K0_decay_rates_poshel[production_modes::k0_e]/total_k0_decay_rate,
              {{-HNL_pdg_poshel,NOT_final_state},{elec_pdg,final_state},{pion_0_pdg,final_state}}
          });
          neutral_kaon.add_decay({
              .5 * K0_decay_rates_neghel[production_modes::k0_e]/total_k0_decay_rate,
              {{HNL_pdg_neghel,NOT_final_state},{-elec_pdg,final_state},{pion_0_pdg,final_state}}
          });
          neutral_kaon.add_decay({
              .5 * K0_decay_rates_neghel[production_modes::k0_e]/total_k0_decay_rate,
              {{-HNL_pdg_neghel,NOT_final_state},{elec_pdg,final_state},{pion_0_pdg,final_state}}
          });
        }
        neutral_kaon.finalise_decay_table();

        
        
        
        
        //////////////// PION DECAYS ///////////////////////////////////////////////////////////////////////////////
        





        double total_pion_decay_rate = 0.;
        std::map<production_modes,double> pi_decay_rates_poshel, pi_decay_rates_neghel;
        if(pion_pm_mass > HNL_mass + muon_mass && params.U_m4 > 0.
            && production_mode_enabled(production_modes::pi_mu)) {
          const double decay_rate_pos = P_decay_rate_to_leptons(
              flavour::m, pion_pm_m2_BR, pion_pm_mass, muon_mass, helicity_plus);
          pi_decay_rates_poshel[production_modes::pi_mu] = decay_rate_pos;
          const double decay_rate_neg = P_decay_rate_to_leptons(
              flavour::m, pion_pm_m2_BR, pion_pm_mass, muon_mass, helicity_minus);
          pi_decay_rates_neghel[production_modes::pi_mu] = decay_rate_neg;
          total_pion_decay_rate += decay_rate_pos + decay_rate_neg;
        }
        if(pion_pm_mass > HNL_mass + elec_mass && params.U_e4 > 0.
            && production_mode_enabled(production_modes::pi_e)) {
          const double decay_rate_pos = P_decay_rate_to_leptons(
              flavour::e, pion_pm_e2_BR, pion_pm_mass, elec_mass, helicity_plus);
          pi_decay_rates_poshel[production_modes::pi_e] = decay_rate_pos;
          const double decay_rate_neg = P_decay_rate_to_leptons(
              flavour::e, pion_pm_e2_BR, pion_pm_mass, elec_mass, helicity_minus);
          pi_decay_rates_neghel[production_modes::pi_e] = decay_rate_neg;
          total_pion_decay_rate += decay_rate_pos + decay_rate_neg;
        }
        
        ret.push_back(dkgen::core::particle_definition{pion_pm_pdg,pion_pm_mass,pion_pm_lt});
        auto& charged_pion = ret.back();
        
        if(pion_pm_mass > HNL_mass + muon_mass && params.U_m4 > 0.
            && production_mode_enabled(production_modes::pi_mu)) {
          charged_pion.add_decay({
              pi_decay_rates_poshel[production_modes::pi_mu]/total_pion_decay_rate,
              {{HNL_pdg_poshel,NOT_final_state},{-muon_pdg,final_state}}
          });
          charged_pion.add_decay({
              pi_decay_rates_neghel[production_modes::pi_mu]/total_pion_decay_rate,
              {{HNL_pdg_neghel,NOT_final_state},{-muon_pdg,final_state}}
          });
        }
        if(pion_pm_mass > HNL_mass + elec_mass && params.U_e4 > 0.
            && production_mode_enabled(production_modes::pi_e)) {
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

        

        //////////////// MUON DECAYS ///////////////////////////////////////////////////////////////////////////////


        auto I_l = [sqrtkl] (double x, double y, double z, int plus_minus, bool is_anti) -> double  {
          std::function<double(double)> integrand;
          if(is_anti) {
            integrand = [sqrtkl,x,y,z,plus_minus](double s) -> double {
              const double sqrt_kl_1sz = sqrtkl(1,s,z);
              return (1.+z-s-plus_minus*sqrt_kl_1sz)*(s-x-y)*sqrtkl(s,x,y)*sqrt_kl_1sz/s;
            };
          } else {
            integrand = [sqrtkl,x,y,z,plus_minus](double s) -> double  {
              const double sqrt_kl_sxy = sqrtkl(s,x,y);
              return (1.+z-s)*(s-x-y-plus_minus*sqrt_kl_sxy)*sqrt_kl_sxy*sqrtkl(1,s,z)/s;
            };
          }
          return 12. * detint::integrate(integrand, std::pow(std::sqrt(x)+std::sqrt(y),2), std::pow(1.-std::sqrt(z),2));
        };
        
        // for mu+ -> e+ + ... (so for antiparticles)
        // split this into UalphaN and UbetaN for later use
        auto triple_lepton_decay_rate = [HNL_mass,I_l] (double U2,
            double parent_mass, double child_mass, int plus_minus, bool anti) {
          const double xi_N = std::pow(HNL_mass/parent_mass,2);
          const double xi_l = std::pow(child_mass/parent_mass,2);
          // only care about relative branching ratios, gFerm * m5/192 pi^3 not needed
          return U2 * (anti ? I_l(0.,xi_l,xi_N,plus_minus,anti) : I_l(xi_N,xi_l,0.,plus_minus,!anti));
        };

        double total_muon_decay_rate = 0.;
        enum class mudecay { mu_e_N, mu_e_Nbar };

        // https://arxiv.org/pdf/1905.00284.pdf p16
        // The decay of a charged lepton (antilepton) of flavour α to a charged lepton (antilepton)
        // of flavour β can be proportional to either |UαN| or |UβN|, producing a heavy Dirac
        // neutrino (antineutrino) in the first case or an antineutrino (neutrino) in the second case
        //
        // Formula 4.4 is for antileptons, hence mu+ -> (antineutrino), e+ -> (neutrino)

        std::map<mudecay,double> muon_decay_rates_poshel, muon_decay_rates_neghel;
        if(muon_mass > HNL_mass + elec_mass && params.U_m4 + params.U_e4 > 0.
            && production_mode_enabled(production_modes::mu_e)) {
          const double decay_rate_pos_Nbar = triple_lepton_decay_rate(dparams.U2(flavour::m),
              muon_mass, elec_mass, helicity_plus, false);
          muon_decay_rates_poshel[mudecay::mu_e_N] = decay_rate_pos_Nbar;
#ifdef DEBUG
          std::cout << "Add decay_rate_pos_Nbar r="<<decay_rate_pos_Nbar<<std::endl;
#endif
          const double decay_rate_pos_N = triple_lepton_decay_rate(dparams.U2(flavour::e),
              muon_mass, elec_mass, helicity_plus, true);
          muon_decay_rates_poshel[mudecay::mu_e_Nbar] = decay_rate_pos_N;
#ifdef DEBUG
          std::cout << "Add decay_rate_pos_N r="<<decay_rate_pos_N<<std::endl;
#endif
          const double decay_rate_neg_Nbar = triple_lepton_decay_rate(dparams.U2(flavour::m),
              muon_mass, elec_mass, helicity_minus, false);
          muon_decay_rates_neghel[mudecay::mu_e_Nbar] = decay_rate_neg_Nbar;
#ifdef DEBUG
          std::cout << "Add decay_rate_neg_Nbar r="<<decay_rate_neg_Nbar<<std::endl;
#endif
          const double decay_rate_neg_N = triple_lepton_decay_rate(dparams.U2(flavour::e),
              muon_mass, elec_mass, helicity_minus, true);
          muon_decay_rates_neghel[mudecay::mu_e_N] = decay_rate_neg_N;
#ifdef DEBUG
          std::cout << "Add decay_rate_neg_N r="<<decay_rate_neg_N<<std::endl;
#endif
          total_muon_decay_rate += decay_rate_pos_N + decay_rate_neg_N + decay_rate_pos_Nbar + decay_rate_neg_Nbar;
        }
        
        ret.push_back(dkgen::core::particle_definition{muon_pdg,muon_mass,muon_lt});
        auto& muon_def = ret.back();
        
        // need to flip HNL helicities as above decay rates were for antimuons
        if(muon_mass > HNL_mass + elec_mass && params.U_m4 + params.U_e4 > 0.
            && production_mode_enabled(production_modes::mu_e)) {
          // deal with -ve rates due to floating point rounding errors
          if(muon_decay_rates_neghel[mudecay::mu_e_N] < 0
              && std::abs(muon_decay_rates_neghel[mudecay::mu_e_N])
              / (total_muon_decay_rate - muon_decay_rates_neghel[mudecay::mu_e_N]) < 1e-8) {
            muon_decay_rates_neghel[mudecay::mu_e_N] = 0.;
          }
          if(muon_decay_rates_neghel[mudecay::mu_e_N] > 0.) {
            muon_def.add_decay({
                muon_decay_rates_neghel[mudecay::mu_e_N]/total_muon_decay_rate,
                {{HNL_pdg_poshel,NOT_final_state},{elec_pdg,final_state},{-nu_pdg,final_state}}
                });
          }
          if(muon_decay_rates_neghel[mudecay::mu_e_Nbar] < 0
              && std::abs(muon_decay_rates_neghel[mudecay::mu_e_Nbar])
              / (total_muon_decay_rate - muon_decay_rates_neghel[mudecay::mu_e_Nbar]) < 1e-8) {
            muon_decay_rates_neghel[mudecay::mu_e_N] = 0.;
          }
          if(muon_decay_rates_neghel[mudecay::mu_e_Nbar] > 0.) {
            muon_def.add_decay({
                muon_decay_rates_neghel[mudecay::mu_e_Nbar]/total_muon_decay_rate,
                {{-HNL_pdg_poshel,NOT_final_state},{elec_pdg,final_state},{nu_pdg,final_state}}
                });
          }
          if(muon_decay_rates_poshel[mudecay::mu_e_N] < 0
              && std::abs(muon_decay_rates_poshel[mudecay::mu_e_N])
              / (total_muon_decay_rate - muon_decay_rates_poshel[mudecay::mu_e_N]) < 1e-8) {
            muon_decay_rates_neghel[mudecay::mu_e_N] = 0.;
          }
          if(muon_decay_rates_poshel[mudecay::mu_e_N] > 0.) {
            muon_def.add_decay({
                muon_decay_rates_poshel[mudecay::mu_e_N]/total_muon_decay_rate,
                {{HNL_pdg_neghel,NOT_final_state},{elec_pdg,final_state},{-nu_pdg,final_state}}
                });
          }
          if(muon_decay_rates_poshel[mudecay::mu_e_Nbar] < 0
              && std::abs(muon_decay_rates_poshel[mudecay::mu_e_Nbar])
              / (total_muon_decay_rate - muon_decay_rates_poshel[mudecay::mu_e_Nbar]) < 1e-8) {
            muon_decay_rates_neghel[mudecay::mu_e_N] = 0.;
          }
          if(muon_decay_rates_poshel[mudecay::mu_e_Nbar] > 0.) {
            muon_def.add_decay({
                muon_decay_rates_poshel[mudecay::mu_e_Nbar]/total_muon_decay_rate,
                {{-HNL_pdg_neghel,NOT_final_state},{elec_pdg,final_state},{nu_pdg,final_state}}
                });
          }
        }
        muon_def.finalise_decay_table();








        

        if(total_kaon_decay_rate + total_k0_decay_rate + total_pion_decay_rate + total_muon_decay_rate == 0.) {
          throw std::runtime_error("Unable to produce HNLs for these parameters and options");
        }





        //////////////// HNL DECAYS ///////////////////////////////////////////////////////////////////////////////






        

        auto I1_2 = [sqrtkl](double x, double y) {
          return sqrtkl(1.,x,y) * (std::pow(1.-x,2) - y*(1.+x));
        };
        auto I1_2_1 = [sqrtkl](double x, double y, double cos_theta, int plus_minus) {
          const double sqrt_kl = sqrtkl(1.,x,y);
          return sqrt_kl / 4./M_PI * (std::pow(1.-x,2) - y*(1.+x) + plus_minus * (x-1.) * sqrt_kl * cos_theta);
        };
        auto I1_3 = [sqrtkl](double x, double y, double z) {
          if(x<0.) x = 0.;
          if(y<0.) y = 0.;
          if(z<0.) z = 0.;
          if(x==0. && y==0. && z==0.) return 1.;
          if((x==0. && y==z) || (y == 0. && z==x) || (z==0. && x==y)) {
            const double w = (x==0.) ? y : ((y==0.) ? z : x);
            // assume w < 0.25 (kinematically allowed) otherwise sqrt will be imaginary
            const double w2 = w*w;
            return std::sqrt(1.-4.*w) * (1.-2.*w * (7. + w + 6.*w2))
              + 24. * w2 * (w2-1.) * (-2. * std::log(1. + std::sqrt(1.-4*w)) + std::log(4*w));
          }
          auto integrand = [sqrtkl,x,y,z](double s) {
            return (s-x-y)*(1+z-s)*sqrtkl(s,x,y)*sqrtkl(1,s,z)/s;
          };
          
          return 12. * detint::integrate(integrand, std::pow(std::sqrt(x)+std::sqrt(y),2), std::pow(1-std::sqrt(z),2));
        };
        auto I2_3 = [sqrtkl](double x, double y, double z) {
          if(x<0.) x = 0.;
          if(y<0.) y = 0.;
          if(z<0.) z = 0.;
          if(y==0. || z==0.) return 0.;
          if(x==0. && y==z) {
            const double w = y;
            return -8. * w * (
                std::sqrt(1.-4.*w)*(w-1.)*(1.+6.*w)
                + 3.*w*(
                  (4.+8.*(w-1.)*w)*std::log(1.+std::sqrt(1.-4.*w))
                  - (2.+4.*(w-1.)*w)*std::log(4.*w)
                  )
                );
          }
          auto integrand = [sqrtkl,x,y,z](double s) {
            return (1+x-s)*sqrtkl(s,y,z)*sqrtkl(1,s,x)/s;
          };
          return 24. * std::sqrt(y*z)
            * detint::integrate(integrand, std::pow(std::sqrt(y)+std::sqrt(z),2), std::pow(1-std::sqrt(x),2));
        };

        
        const double decay_rate_to_3nu = dparams.sum_U2 * gfermi2 * dparams.m5 * pi3_inv / 96.;

        auto decay_rate_to_nu_pseudoscalar = [&dparams, gfermi2, HNL_mass, I1_2](double pseudoscalar_mass, double f2) {
          const double xi_p = std::pow(pseudoscalar_mass / HNL_mass, 2);
          return dparams.sum_U2 * gfermi2 * f2 * dparams.m3 * I1_2(0.,xi_p) / 16./ M_PI;
        };
        
        auto decay_rate_to_lepton_pseudoscalar =
          [&dparams, gfermi2, I1_2](double lep_mass, double pseudoscalar_mass, double f2,
              double CKM_V2, double U2) {
          const double xi_l = std::pow(lep_mass / dparams.raw_params.HNL_mass,2);
          const double xi_p = std::pow(pseudoscalar_mass / dparams.raw_params.HNL_mass,2);
          return U2 * CKM_V2 * gfermi2 * f2 * dparams.m3 * I1_2(xi_l,xi_p) / 16./ M_PI;
        };
        auto make_diff_decay_rate_function_to_lepton_pseudoscalar = [I1_2,I1_2_1,HNL_mass=dparams.raw_params.HNL_mass](
          double lep_mass, double pseudoscalar_mass, int plus_minus /* helicity */) {
          const double xi_l = std::pow(lep_mass / HNL_mass,2);
          const double xi_p = std::pow(pseudoscalar_mass / HNL_mass,2);
          auto diff_decay_rate = [I1_2_1,i1_2=I1_2(xi_l, xi_p),xi_l,xi_p,plus_minus](double cos_theta_l) {
            return I1_2_1(xi_l, xi_p, cos_theta_l, plus_minus)/i1_2;
          };
          return diff_decay_rate;
        };

        auto decay_rate_to_nu_2lep = [&dparams, gfermi2, HNL_mass, I1_3, I2_3](
            double lep_minus_mass, double lep_plus_mass, double c1, double c2, double c3) {
          const double xi_m = std::pow(lep_minus_mass / HNL_mass,2);
          const double xi_p = std::pow(lep_plus_mass / HNL_mass,2);
          if(lep_minus_mass == lep_plus_mass) {
            return gfermi2 * dparams.m5 * pi3_inv / 192.
            * ((c1+c2) * I1_3(0,xi_m,xi_p) + c3 * I2_3(0,xi_m,xi_p));
          }
          return gfermi2 * dparams.m5 * pi3_inv / 192.
            * (c1 * I1_3(0,xi_m,xi_p) + c2 * I1_3(0,xi_p,xi_m) + c3 * I2_3(0,xi_m,xi_p));
        };
        auto make_diff_decay_rate_function_to_nu_2lep = [sqrtkl,HNL_mass,I1_3,I2_3](double lep_minus_mass, double lep_plus_mass,
            double c1, double c2, double c3, double c4, double c5, double c6, int plus_minus /* helicity */) {
          const double xi_m = std::pow(lep_minus_mass / HNL_mass,2);
          const double xi_p = std::pow(lep_plus_mass / HNL_mass,2);
          const double total_rate = (lep_minus_mass == lep_plus_mass)
            ?
            ((c1+c2) * I1_3(0,xi_m,xi_p) + c3 * I2_3(0,xi_m,xi_p)) 
            :
            (c1 * I1_3(0,xi_m,xi_p) + c2 * I1_3(0,xi_p,xi_m) + c3 * I2_3(0,xi_m,xi_p));
          auto diff_decay_rate = [sqrtkl,xi_m, xi_p, c1, c2, c3, c4, c5, c6, plus_minus, total_rate]
            (double s1, double s2, double cos_theta_m,
             double /*dummy phi_m*/, double cos_theta_p, double /*dummy phi_mp*/) -> double {
            const double A02 = c1 * (s2 - xi_m)*(1+xi_p-s2) + c2 * (s1 - xi_p)*(1+xi_m-s1) + 2.*c3*xi_m*xi_p*(s1+s2-xi_m-xi_p);
            const double A12 = (c4 *(s2-xi_m) -2.*c6*xi_m*xi_p)*sqrtkl(1.,s2,xi_p)*cos_theta_p
              + (c5*(s1-xi_p)-2.*c6*xi_m*xi_p)*sqrtkl(1.,s1,xi_m)*cos_theta_m;
            return (A02 + plus_minus * A12)/total_rate/12.;
          };
          return diff_decay_rate;
        };

        const double gL = conf.physical_params().sin2thW  - 0.5;
        const double gR = conf.physical_params().sin2thW;
        auto c1_nu_dirac = [&params,gL](flavour l1, flavour l2) { 
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U = (fl==flavour::e ? params.U_e4 : (fl==flavour::m? params.U_m4 : params.U_t4));
            //                                        vv  ----- l1 (cf l2 in c2_nubar_dirac)
            tot += U*U *((l1==l2 ? gL*gL : 0.) + (fl==l1?1:0)*(1.+(l1==l2?gL:0.)));
          }
          return tot;
        };
        auto c1_nubar_dirac = [&dparams,gR](flavour l1, flavour l2) {
          return l1==l2 ? dparams.sum_U2 * gR*gR : 0.;
        };
        auto c2_nu_dirac = [&dparams,gR](flavour l1, flavour l2) {
          return l1==l2 ? dparams.sum_U2 * gR*gR : 0.;
        };
        auto c2_nubar_dirac = [&params,gL](flavour l1, flavour l2) { 
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U = (fl==flavour::e ? params.U_e4 : (fl==flavour::m? params.U_m4 : params.U_t4));
            //                                        vv  ----- l2 (cf l1 in c1_nu_dirac)
            tot += U*U *((l1==l2 ? gL*gL : 0.) + (fl==l2?1:0)*(1.+(l1==l2?gL:0.)));
          }
          return tot;
        };
        auto c3_nu_dirac = [&params,gL,gR](flavour l1, flavour l2) { 
          if(l1 != l2) return 0.;
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U = (fl==flavour::e ? params.U_e4 : (fl==flavour::m? params.U_m4 : params.U_t4));
            //                 vv  ----- l2 (cf l1 in c3_nubar_dirac)
            tot += U*U * ((fl==l2? 1. : 0.) + gL);
          }
          return tot*gR;
        };
        auto c3_nubar_dirac = [&params,gL,gR](flavour l1, flavour l2) { 
          if(l1 != l2) return 0.;
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U = (fl==flavour::e ? params.U_e4 : (fl==flavour::m? params.U_m4 : params.U_t4));
            //                 vv  ----- l1 (cf l2 in c3_nu_dirac)
            tot += U*U * ((fl==l1? 1. : 0.) + gL);
          }
          return tot*gR;
        };


        auto c4_nu_dirac = c1_nu_dirac;
        auto c5_nu_dirac = c2_nu_dirac;
        auto c6_nu_dirac = c3_nu_dirac;
        auto c4_nubar_dirac = [c1_nubar_dirac](flavour l1, flavour l2) { return -c1_nubar_dirac(l1,l2); };
        auto c5_nubar_dirac = [c2_nubar_dirac](flavour l1, flavour l2) { return -c2_nubar_dirac(l1,l2); };
        auto c6_nubar_dirac = [c3_nubar_dirac](flavour l1, flavour l2) { return -c3_nubar_dirac(l1,l2); };
        
        auto c1_nu_major = [c1_nu_dirac,c1_nubar_dirac](flavour l1, flavour l2) {
          return c1_nu_dirac(l1,l2) + c1_nubar_dirac(l1,l2); };
        auto c2_nu_major = [c2_nu_dirac,c2_nubar_dirac](flavour l1, flavour l2) {
          return c2_nu_dirac(l1,l2) + c2_nubar_dirac(l1,l2); };
        auto c3_nu_major = [c3_nu_dirac,c3_nubar_dirac](flavour l1, flavour l2) {
          return c3_nu_dirac(l1,l2) + c3_nubar_dirac(l1,l2); };
        auto c4_nu_major = [c4_nu_dirac,c4_nubar_dirac](flavour l1, flavour l2) {
          return c4_nu_dirac(l1,l2) - c4_nubar_dirac(l1,l2); };
        auto c5_nu_major = [c5_nu_dirac,c5_nubar_dirac](flavour l1, flavour l2) {
          return c5_nu_dirac(l1,l2) - c5_nubar_dirac(l1,l2); };
        auto c6_nu_major = [c6_nu_dirac,c6_nubar_dirac](flavour l1, flavour l2) {
          return c6_nu_dirac(l1,l2) - c6_nubar_dirac(l1,l2); };

        auto c1_final = [c1_nu_major,c1_nu_dirac](bool majorana, flavour l1, flavour l2) {
          return majorana ? c1_nu_major(l1,l2) : c1_nu_dirac(l1,l2);
        };
        auto c2_final = [c2_nu_major,c2_nu_dirac](bool majorana, flavour l1, flavour l2) {
          return majorana ? c2_nu_major(l1,l2) : c2_nu_dirac(l1,l2);
        };
        auto c3_final = [c3_nu_major,c3_nu_dirac](bool majorana, flavour l1, flavour l2) {
          return majorana ? c3_nu_major(l1,l2) : c3_nu_dirac(l1,l2);
        };
        auto c4_final = [c4_nu_major,c4_nu_dirac](bool majorana, flavour l1, flavour l2) {
          return majorana ? c4_nu_major(l1,l2) : c4_nu_dirac(l1,l2);
        };
        auto c5_final = [c5_nu_major,c5_nu_dirac](bool majorana, flavour l1, flavour l2) {
          return majorana ? c5_nu_major(l1,l2) : c5_nu_dirac(l1,l2);
        };
        auto c6_final = [c6_nu_major,c6_nu_dirac](bool majorana, flavour l1, flavour l2) {
          return majorana ? c6_nu_major(l1,l2) : c6_nu_dirac(l1,l2);
        };


        // total decay rates are helicity-independent (HNL_lifetime will apply to both helicities)
        std::map<decay_modes, double> decay_rates;
        // always have 3nu mode
        double total_decay_rate = (params.is_majorana ? 1. : 0.5) * decay_rate_to_3nu;
        decay_rates[decay_modes::nu_nu_nu] = total_decay_rate;
#ifdef DEBUG
          std::cout << "Add hnl->nu_nu_nu r="<<total_decay_rate<<std::endl;
#endif
        if(HNL_mass > 2*elec_mass) {
          const double decay_rate = params.is_majorana
            ?
            decay_rate_to_nu_2lep(elec_mass, elec_mass,
                c1_nu_major(flavour::e,flavour::e),
                c2_nu_major(flavour::e,flavour::e),
                c3_nu_major(flavour::e,flavour::e))
            :
            decay_rate_to_nu_2lep(elec_mass, elec_mass,
                c1_nu_dirac(flavour::e,flavour::e),
                c2_nu_dirac(flavour::e,flavour::e),
                c3_nu_dirac(flavour::e,flavour::e));
          total_decay_rate += decay_rate;
          decay_rates[decay_modes::e_e_nu] = decay_rate;
#ifdef DEBUG
          std::cout << "Add hnl->e_e_nu r="<<decay_rate<<std::endl;
#endif
        }
        if(HNL_mass > elec_mass + muon_mass) {
          const double decay_rate1 = params.is_majorana
            ? decay_rate_to_nu_2lep(elec_mass, muon_mass,
                c1_nu_major(flavour::e,flavour::m),
                c2_nu_major(flavour::e,flavour::m),
                c3_nu_major(flavour::e,flavour::m))
            : decay_rate_to_nu_2lep(elec_mass, muon_mass,
                c1_nu_dirac(flavour::e,flavour::m),
                c2_nu_dirac(flavour::e,flavour::m),
                c3_nu_dirac(flavour::e,flavour::m));
          const double decay_rate2 = params.is_majorana
            ? decay_rate_to_nu_2lep(muon_mass, elec_mass,
                c1_nu_major(flavour::m,flavour::e),
                c2_nu_major(flavour::m,flavour::e),
                c3_nu_major(flavour::m,flavour::e))
            : decay_rate_to_nu_2lep(muon_mass, elec_mass,
                c1_nu_dirac(flavour::m,flavour::e),
                c2_nu_dirac(flavour::m,flavour::e),
                c3_nu_dirac(flavour::m,flavour::e));
          total_decay_rate += decay_rate1 + decay_rate2;
          decay_rates[decay_modes::e_mu_nu] = decay_rate1;
          decay_rates[decay_modes::mu_e_nu] = decay_rate2;
#ifdef DEBUG
          std::cout << "Add hnl->e_mu_nu r="<<decay_rate1<<std::endl;
          std::cout << "Add hnl->mu_e_nu r="<<decay_rate2<<std::endl;
#endif
        }
        if(HNL_mass > pion_0_mass) {
          const double decay_rate = (params.is_majorana ? 1. : 0.5) * decay_rate_to_nu_pseudoscalar(pion_0_mass, f2_pion);
          total_decay_rate += decay_rate;
          decay_rates[decay_modes::pi0_nu] = decay_rate;
#ifdef DEBUG
          std::cout << "Add hnl->pi0_nu r="<<decay_rate<<std::endl;
#endif
        }
        if(HNL_mass > elec_mass+pion_pm_mass) {
          const double decay_rate = (params.is_majorana ? 2. : 1.) * decay_rate_to_lepton_pseudoscalar(
              elec_mass, pion_pm_mass, f2_pion, CKM_Vud2, params.U_e4*params.U_e4);
          total_decay_rate += decay_rate;
          decay_rates[decay_modes::e_pi] = decay_rate;
#ifdef DEBUG
          std::cout << "Add hnl->e_pi r="<<decay_rate<<std::endl;
#endif
        }
        if(HNL_mass > 2*muon_mass) {
          const double decay_rate = params.is_majorana
            ?
            decay_rate_to_nu_2lep(muon_mass, muon_mass,
                c1_nu_major(flavour::m,flavour::m),
                c2_nu_major(flavour::m,flavour::m),
                c3_nu_major(flavour::m,flavour::m))
            :
            decay_rate_to_nu_2lep(muon_mass, muon_mass,
                c1_nu_dirac(flavour::m,flavour::m),
                c2_nu_dirac(flavour::m,flavour::m),
                c3_nu_dirac(flavour::m,flavour::m));
          total_decay_rate += decay_rate;
          decay_rates[decay_modes::mu_mu_nu] = decay_rate;
#ifdef DEBUG
          std::cout << "Add hnl->mu_mu_nu r="<<decay_rate<<std::endl;
#endif
        }
        if(HNL_mass > muon_mass+pion_pm_mass) {
          const double decay_rate = (params.is_majorana ? 2. : 1.) * decay_rate_to_lepton_pseudoscalar(
              muon_mass, pion_pm_mass, f2_pion, CKM_Vud2, params.U_m4*params.U_m4);
          total_decay_rate += decay_rate;
          decay_rates[decay_modes::mu_pi] = decay_rate;
#ifdef DEBUG
          std::cout << "Add hnl->mu_pi r="<<decay_rate<<std::endl;
#endif
        }
        // total decay rate is helicity-independent according to https://arxiv.org/abs/1905.00284
        const double HNL_lifetime = hbar / total_decay_rate;


        // have to define separate helicity states, so that they can have different angular distributions
        ret.push_back(dkgen::core::particle_definition{
            HNL_pdg_poshel, HNL_mass, HNL_lifetime, params.is_majorana? self_conjugate : !self_conjugate
            });
        const size_t HNL_poshel_loc = ret.size()-1;
        ret.push_back(dkgen::core::particle_definition{
            HNL_pdg_neghel, HNL_mass, HNL_lifetime, params.is_majorana? self_conjugate : !self_conjugate
            });
        const size_t HNL_neghel_loc = ret.size()-1;

        // cannot directly use references to ret.back(), because push_back() can invalidate iterators
        auto& HNL_poshel_info = *std::next(ret.begin(), HNL_poshel_loc);
        auto& HNL_neghel_info = *std::next(ret.begin(), HNL_neghel_loc);

        auto decay_mode_enabled = [decay_modes_to_use](auto dm) -> bool {
          return std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), dm) != decay_modes_to_use.end();
        };

        // don't simulate triple-neutrino decay mode
        /*
        if(HNL_mass > 0.) {
          HNL_poshel_info.add_decay(
              { decay_rates[decay_modes::nu_nu_nu]/total_decay_rate,
                { {nu_pdg,final_state},{nu_pdg,final_state},{nu_pdg,final_state} } // isotropic kinematics
              });
        }
        */

        if(HNL_mass > 2*elec_mass && decay_mode_enabled(decay_modes::e_e_nu)) {
          if(decay_rates[decay_modes::e_e_nu] > 0.) {
            const double c1 = c1_final(params.is_majorana, flavour::e, flavour::e);
            const double c2 = c2_final(params.is_majorana, flavour::e, flavour::e);
            const double c3 = c3_final(params.is_majorana, flavour::e, flavour::e);
            const double c4 = c4_final(params.is_majorana, flavour::e, flavour::e);
            const double c5 = c5_final(params.is_majorana, flavour::e, flavour::e);
            const double c6 = c6_final(params.is_majorana, flavour::e, flavour::e);
            // double lep_minus_mass, double lep_plus_mass, double c1, double c2, double c3, double c4, double c5, double c6
            auto rw_fn_pos = make_diff_decay_rate_function_to_nu_2lep(elec_mass, elec_mass, c1, c2, c3, c4, c5, c6, helicity_plus);
            auto rw_fn_neg = make_diff_decay_rate_function_to_nu_2lep(elec_mass, elec_mass, c1, c2, c3, c4, c5, c6, helicity_minus);
            HNL_poshel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::e_e_nu]/total_decay_rate,
                { {nu_pdg,final_state},{elec_pdg,final_state},{-elec_pdg,final_state} }
                }.set_angular_threebody_dalitz_reweighter(core::angular_threebody_dalitz_function{rw_fn_pos}));
            HNL_neghel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::e_e_nu]/total_decay_rate,
                { {-nu_pdg,final_state},{elec_pdg,final_state},{-elec_pdg,final_state} }
                }.set_angular_threebody_dalitz_reweighter(core::angular_threebody_dalitz_function{rw_fn_neg}));
          }
        }
        
        if(HNL_mass > elec_mass + muon_mass
            && (decay_mode_enabled(decay_modes::e_mu_nu) || decay_mode_enabled(decay_modes::mu_e_nu))) {
          if(decay_rates[decay_modes::e_mu_nu] > 0.) {
            const double c1 = c1_final(params.is_majorana, flavour::e, flavour::m);
            const double c2 = c2_final(params.is_majorana, flavour::e, flavour::m);
            const double c3 = c3_final(params.is_majorana, flavour::e, flavour::m);
            const double c4 = c4_final(params.is_majorana, flavour::e, flavour::m);
            const double c5 = c5_final(params.is_majorana, flavour::e, flavour::m);
            const double c6 = c6_final(params.is_majorana, flavour::e, flavour::m);
            auto rw_fn_pos = make_diff_decay_rate_function_to_nu_2lep(elec_mass, muon_mass, c1, c2, c3, c4, c5, c6, helicity_plus);
            auto rw_fn_neg = make_diff_decay_rate_function_to_nu_2lep(elec_mass, muon_mass, c1, c2, c3, c4, c5, c6, helicity_minus);
            HNL_poshel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::e_mu_nu]/total_decay_rate,
                {{nu_pdg,final_state},{elec_pdg,final_state},{-muon_pdg,final_state}}
                }.set_angular_threebody_dalitz_reweighter(core::angular_threebody_dalitz_function{rw_fn_pos}));
            HNL_neghel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::e_mu_nu]/total_decay_rate,
                {{-nu_pdg,final_state},{elec_pdg,final_state},{-muon_pdg,final_state}}
                }.set_angular_threebody_dalitz_reweighter(core::angular_threebody_dalitz_function{rw_fn_neg}));
          }
          if(decay_rates[decay_modes::mu_e_nu] > 0.) {
            const double c1 = c1_final(params.is_majorana, flavour::m, flavour::e);
            const double c2 = c2_final(params.is_majorana, flavour::m, flavour::e);
            const double c3 = c3_final(params.is_majorana, flavour::m, flavour::e);
            const double c4 = c4_final(params.is_majorana, flavour::m, flavour::e);
            const double c5 = c5_final(params.is_majorana, flavour::m, flavour::e);
            const double c6 = c6_final(params.is_majorana, flavour::m, flavour::e);
            auto rw_fn_pos = make_diff_decay_rate_function_to_nu_2lep(muon_mass, elec_mass, c1, c2, c3, c4, c5, c6, helicity_plus);
            auto rw_fn_neg = make_diff_decay_rate_function_to_nu_2lep(muon_mass, elec_mass, c1, c2, c3, c4, c5, c6, helicity_minus);
            HNL_poshel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::mu_e_nu]/total_decay_rate,
                {{nu_pdg,final_state},{muon_pdg,final_state},{-elec_pdg,final_state}}
                }.set_angular_threebody_dalitz_reweighter(core::angular_threebody_dalitz_function{rw_fn_pos}));
            HNL_neghel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::mu_e_nu]/total_decay_rate,
                {{-nu_pdg,final_state},{muon_pdg,final_state},{-elec_pdg,final_state}}
                }.set_angular_threebody_dalitz_reweighter(core::angular_threebody_dalitz_function{rw_fn_neg}));
          }
        }
        
        if(HNL_mass > pion_0_mass && decay_mode_enabled(decay_modes::pi0_nu)) {
          if(decay_rates[decay_modes::pi0_nu] > 0.) {
            if(params.is_majorana) {
              HNL_poshel_info.add_decay({
                  decay_rates[decay_modes::pi0_nu]/total_decay_rate,
                  {{nu_pdg,final_state},{pion_0_pdg,final_state}}
                  });
              HNL_neghel_info.add_decay({
                  decay_rates[decay_modes::pi0_nu]/total_decay_rate,
                  {{-nu_pdg,final_state},{pion_0_pdg,final_state}}
                  });
            } else {
              auto rw_fn_pos = make_diff_decay_rate_function_to_lepton_pseudoscalar(0., pion_0_mass, helicity_plus);
              auto rw_fn_neg = make_diff_decay_rate_function_to_lepton_pseudoscalar(0., pion_0_mass, helicity_minus);
              HNL_poshel_info.add_decay(
                  core::decay_mode{ decay_rates[decay_modes::pi0_nu]/total_decay_rate,
                  {{nu_pdg,final_state},{pion_0_pdg,final_state}}
                  }.set_twobody_dalitz_reweighter(core::twobody_dalitz_function{rw_fn_pos}));
              HNL_neghel_info.add_decay(
                  core::decay_mode{ decay_rates[decay_modes::pi0_nu]/total_decay_rate,
                  {{-nu_pdg,final_state},{pion_0_pdg,final_state}}
                  }.set_twobody_dalitz_reweighter(core::twobody_dalitz_function{rw_fn_neg}));
            }
          }
        }
        
        if(HNL_mass > elec_mass+pion_pm_mass && decay_mode_enabled(decay_modes::e_pi)) {
          if(decay_rates[decay_modes::e_pi] > 0.) {
            auto rw_fn_pos = make_diff_decay_rate_function_to_lepton_pseudoscalar(elec_mass, pion_pm_mass, helicity_plus);
            auto rw_fn_neg = make_diff_decay_rate_function_to_lepton_pseudoscalar(elec_mass, pion_pm_mass, helicity_minus);
            HNL_poshel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::e_pi]/total_decay_rate,
                {{elec_pdg,final_state},{pion_pm_pdg,final_state}}
                }.set_twobody_dalitz_reweighter(core::twobody_dalitz_function{rw_fn_pos}));
            HNL_neghel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::e_pi]/total_decay_rate,
                {{elec_pdg,final_state},{pion_pm_pdg,final_state}}
                }.set_twobody_dalitz_reweighter(core::twobody_dalitz_function{rw_fn_neg}));
          }
        }

        if(HNL_mass > 2*muon_mass && decay_mode_enabled(decay_modes::mu_mu_nu)) {
          if(decay_rates[decay_modes::mu_mu_nu] > 0.) {
            const double c1 = c1_final(params.is_majorana, flavour::m, flavour::m);
            const double c2 = c2_final(params.is_majorana, flavour::m, flavour::m);
            const double c3 = c3_final(params.is_majorana, flavour::m, flavour::m);
            const double c4 = c4_final(params.is_majorana, flavour::m, flavour::m);
            const double c5 = c5_final(params.is_majorana, flavour::m, flavour::m);
            const double c6 = c6_final(params.is_majorana, flavour::m, flavour::m);
            auto rw_fn_pos = make_diff_decay_rate_function_to_nu_2lep(
                muon_mass, muon_mass, c1, c2, c3, c4, c5, c6, helicity_plus);
            auto rw_fn_neg = make_diff_decay_rate_function_to_nu_2lep(
                muon_mass, muon_mass, c1, c2, c3, c4, c5, c6, helicity_minus);
            HNL_poshel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::mu_mu_nu]/total_decay_rate,
                {{nu_pdg,final_state},{muon_pdg,final_state},{-muon_pdg,final_state}}
                }.set_angular_threebody_dalitz_reweighter(core::angular_threebody_dalitz_function{rw_fn_pos}));
            HNL_neghel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::mu_mu_nu]/total_decay_rate,
                {{-nu_pdg,final_state},{muon_pdg,final_state},{-muon_pdg,final_state}}
                }.set_angular_threebody_dalitz_reweighter(core::angular_threebody_dalitz_function{rw_fn_neg}));
          }
        }

        if(HNL_mass > muon_mass+pion_pm_mass && decay_mode_enabled(decay_modes::mu_pi)) {
          if(decay_rates[decay_modes::mu_pi] > 0.) {
            auto rw_fn_pos = make_diff_decay_rate_function_to_lepton_pseudoscalar(muon_mass, pion_pm_mass, helicity_plus);
            auto rw_fn_neg = make_diff_decay_rate_function_to_lepton_pseudoscalar(muon_mass, pion_pm_mass, helicity_minus);
            HNL_poshel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::mu_pi]/total_decay_rate,
                {{muon_pdg,final_state},{pion_pm_pdg,final_state}}
                }.set_twobody_dalitz_reweighter(core::twobody_dalitz_function{rw_fn_pos}));
            HNL_neghel_info.add_decay(
                core::decay_mode{ decay_rates[decay_modes::mu_pi]/total_decay_rate,
                {{muon_pdg,final_state},{pion_pm_pdg,final_state}}
                }.set_twobody_dalitz_reweighter(core::twobody_dalitz_function{rw_fn_neg}));
          }

        }
        if(HNL_poshel_info.get_decay_table().empty() && HNL_neghel_info.get_decay_table().empty()) {
          throw std::runtime_error("Unable to have decaying HNLs for these parameters and options");
        }

        HNL_poshel_info.finalise_decay_table();
        HNL_neghel_info.finalise_decay_table();
        
        
        
        
        
        //////////////// REMAINING PARTICLES ////////////////////////////////////////////////////////////////////
        
        
        
        
        
        
        ret.push_back(dkgen::core::particle_definition{pion_0_pdg,pion_0_mass,pion_0_lt,self_conjugate});
        ret.push_back(dkgen::core::particle_definition{elec_pdg,elec_mass,elec_lt});
        ret.push_back(dkgen::core::particle_definition{nu_pdg,0.,-1.});

        return ret;
      }

      /*
      double kaon_branching_ratio(const model_parameters& params, int kaon_pdg,
          const dkgen::core::config& conf) {

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
      */
    }
  }
}

#endif
