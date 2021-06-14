#ifndef __dkgen_physics_detail_hnl_decays_hpp__
#define __dkgen_physics_detail_hnl_decays_hpp__


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

              
        struct c_struct {
          double c1, c2, c3, c4, c5, c6;
        };

        double I1_2(double x, double y) {
          return utils::sqrtkl(1.,x,y) * (std::pow(1.-x,2) - y*(1.+x));
        }
        double I1_2_1(double x, double y, double cos_theta, int plus_minus) {
          const double sqrt_kl = utils::sqrtkl(1.,x,y);
          return sqrt_kl / 4./M_PI * (std::pow(1.-x,2) - y*(1.+x) + plus_minus * (x-1.) * sqrt_kl * cos_theta);
        }
        double I1_3(double x, double y, double z) {
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
          auto integrand = [x,y,z](double s) {
            return (s-x-y)*(1+z-s)*utils::sqrtkl(s,x,y)*utils::sqrtkl(1,s,z)/s;
          };

          return 12. * detail::integration::integrate(integrand, std::pow(std::sqrt(x)+std::sqrt(y),2), std::pow(1-std::sqrt(z),2));
        }
        double I2_3(double x, double y, double z) {
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
          auto integrand = [x,y,z](double s) {
            return (1+x-s)*utils::sqrtkl(s,y,z)*utils::sqrtkl(1,s,x)/s;
          };
          return 24. * std::sqrt(y*z)
            * detail::integration::integrate(integrand, std::pow(std::sqrt(y)+std::sqrt(z),2), std::pow(1-std::sqrt(x),2));
        }


        double decay_rate_to_3nu(const derived_params& dparams) {
          return dparams.sum_U2 * dparams.gFermi2 * dparams.m5 * dparams.pi3_inv / 96.;
        }

        double decay_rate_to_nu_pseudoscalar(const derived_params& dparams,
            double pseudoscalar_mass, double f2) {
          const double xi_p = std::pow(pseudoscalar_mass / dparams.raw_params.HNL_mass, 2);
          return dparams.sum_U2 * dparams.gFermi2 * f2 * dparams.m3 * std::pow(1-xi_p,2) / 64./ M_PI;
        }

        double decay_rate_to_lepton_pseudoscalar(const derived_params& dparams,
            double lep_mass, double pseudoscalar_mass, double f2,
            double CKM_V2, double U2) {
          const double xi_l = std::pow(lep_mass / dparams.raw_params.HNL_mass,2);
          const double xi_p = std::pow(pseudoscalar_mass / dparams.raw_params.HNL_mass,2);
          return U2 * CKM_V2 * dparams.gFermi2 * f2 * dparams.m3 * I1_2(xi_l,xi_p) / 16./ M_PI;
        }
        auto make_diff_decay_rate_function_to_lepton_pseudoscalar(const derived_params& dparams,
            double lep_mass, double pseudoscalar_mass, int plus_minus /* helicity */) {
          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double xi_l = std::pow(lep_mass / HNL_mass,2);
          const double xi_p = std::pow(pseudoscalar_mass / HNL_mass,2);
          auto diff_decay_rate = [i1_2=I1_2(xi_l, xi_p),xi_l,xi_p,plus_minus](double cos_theta_l) {
            return I1_2_1(xi_l, xi_p, cos_theta_l, plus_minus)/i1_2;
          };
          return diff_decay_rate;
        }

        double decay_rate_to_nu_2lep(const derived_params& dparams, double lep_minus_mass,
            double lep_plus_mass, const c_struct& cs) {
          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double xi_m = std::pow(lep_minus_mass / HNL_mass,2);
          const double xi_p = std::pow(lep_plus_mass / HNL_mass,2);
          if(lep_minus_mass == lep_plus_mass) {
            return dparams.gFermi2 * dparams.m5 * dparams.pi3_inv / 192.
              * ((cs.c1 + cs.c2) * I1_3(0,xi_m,xi_p) + cs.c3 * I2_3(0,xi_m,xi_p));
          }
          return dparams.gFermi2 * dparams.m5 * dparams.pi3_inv / 192.
            * (cs.c1 * I1_3(0,xi_m,xi_p) + cs.c2 * I1_3(0,xi_p,xi_m) + cs.c3 * I2_3(0,xi_m,xi_p));
        }
        auto make_diff_decay_rate_function_to_nu_2lep(const derived_params& dparams, double lep_minus_mass, double lep_plus_mass,
            const c_struct& cs, int plus_minus /* helicity */) {
          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double xi_m = std::pow(lep_minus_mass / HNL_mass,2);
          const double xi_p = std::pow(lep_plus_mass / HNL_mass,2);
          const double total_rate = (lep_minus_mass == lep_plus_mass)
            ?
            ((cs.c1+cs.c2) * I1_3(0,xi_m,xi_p) + cs.c3 * I2_3(0,xi_m,xi_p)) 
            :
            (cs.c1 * I1_3(0,xi_m,xi_p) + cs.c2 * I1_3(0,xi_p,xi_m) + cs.c3 * I2_3(0,xi_m,xi_p));
          auto diff_decay_rate =
            [xi_m, xi_p, c1=cs.c1, c2=cs.c2, c3=cs.c3, c4=cs.c4, c5=cs.c5, c6=cs.c6, plus_minus, total_rate]
            (double s1, double s2, double cos_theta_m,
             double /*dummy phi_m*/, double cos_theta_p, double /*dummy phi_mp*/) -> double {
              const double A02 = c1 * (s2 - xi_m)*(1+xi_p-s2) + c2 * (s1 - xi_p)*(1+xi_m-s1) + 2.*c3*xi_m*xi_p*(s1+s2-xi_m-xi_p);
              const double A12 = (c4 *(s2-xi_m) -2.*c6*xi_m*xi_p)*utils::sqrtkl(1.,s2,xi_p)*cos_theta_p
                + (c5*(s1-xi_p)-2.*c6*xi_m*xi_p)*utils::sqrtkl(1.,s1,xi_m)*cos_theta_m;
              return (A02 + plus_minus * A12)/total_rate/12.;
            };
          return diff_decay_rate;
        }

        double c1_nu_dirac(const derived_params& dparams, flavour l1, flavour l2) { 
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U2 = dparams.U2(fl);
            //                                       vv  ----- l1 (cf l2 in c2_nubar_dirac)
            tot += U2 *((l1==l2 ? dparams.gL*dparams.gL : 0.) + (fl==l1?1:0)*(1.+(l1==l2?2.*dparams.gL:0.)));
            // (1 + 2 gL) in all the other formulas, unlike 1905.00284.pdf
          }
          return tot;
        }
        double c1_nubar_dirac(const derived_params& dparams, flavour l1, flavour l2) {
          return l1==l2 ? dparams.sum_U2 * dparams.gR*dparams.gR : 0.;
        }
        double c2_nu_dirac(const derived_params& dparams, flavour l1, flavour l2) {
          return l1==l2 ? dparams.sum_U2 * dparams.gR*dparams.gR : 0.;
        }
        double c2_nubar_dirac(const derived_params& dparams, flavour l1, flavour l2) { 
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U2 = dparams.U2(fl);
            //                                       vv  ----- l2 (cf l1 in c1_nu_dirac)
            tot += U2 *((l1==l2 ? dparams.gL*dparams.gL : 0.) + (fl==l2?1:0)*(1.+(l1==l2?2.*dparams.gL:0.)));
            // (1 + 2 gL) in all the other formulas, unlike 1905.00284.pdf
          }
          return tot;
        }
        double c3_nu_dirac(const derived_params& dparams, flavour l1, flavour l2) { 
          if(l1 != l2) return 0.;
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U2 = dparams.U2(fl);
            //                vv  ----- l2 (cf l1 in c3_nubar_dirac)
            tot += U2 * ((fl==l2? 1. : 0.) + dparams.gL);
          }
          return tot*dparams.gR;
        }
        double c3_nubar_dirac(const derived_params& dparams, flavour l1, flavour l2) { 
          if(l1 != l2) return 0.;
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U2 = dparams.U2(fl);
            //                vv  ----- l1 (cf l2 in c3_nu_dirac)
            tot += U2 * ((fl==l1? 1. : 0.) + dparams.gL);
          }
          return tot*dparams.gR;
        }

        auto c4_nu_dirac = c1_nu_dirac;
        auto c5_nu_dirac = c2_nu_dirac;
        auto c6_nu_dirac = c3_nu_dirac;
        double c4_nubar_dirac(const derived_params& dparams, flavour l1, flavour l2) { return -c1_nubar_dirac(dparams,l1,l2); }
        double c5_nubar_dirac(const derived_params& dparams, flavour l1, flavour l2) { return -c2_nubar_dirac(dparams,l1,l2); }
        double c6_nubar_dirac(const derived_params& dparams, flavour l1, flavour l2) { return -c3_nubar_dirac(dparams,l1,l2); }

        double c1_nu_major(const derived_params& dparams, flavour l1, flavour l2) {
          return c1_nu_dirac(dparams,l1,l2) + c1_nubar_dirac(dparams,l1,l2); }
        double c2_nu_major(const derived_params& dparams, flavour l1, flavour l2) {
          return c2_nu_dirac(dparams,l1,l2) + c2_nubar_dirac(dparams,l1,l2); }
        double c3_nu_major(const derived_params& dparams, flavour l1, flavour l2) {
          return c3_nu_dirac(dparams,l1,l2) + c3_nubar_dirac(dparams,l1,l2); }
        double c4_nu_major(const derived_params& dparams, flavour l1, flavour l2) {
          return c4_nu_dirac(dparams,l1,l2) - c4_nubar_dirac(dparams,l1,l2); }
        double c5_nu_major(const derived_params& dparams, flavour l1, flavour l2) {
          return c5_nu_dirac(dparams,l1,l2) - c5_nubar_dirac(dparams,l1,l2); }
        double c6_nu_major(const derived_params& dparams, flavour l1, flavour l2) {
          return c6_nu_dirac(dparams,l1,l2) - c6_nubar_dirac(dparams,l1,l2); }

        double c1_final(const derived_params& dparams, bool majorana, flavour l1, flavour l2) {
          return majorana ? c1_nu_major(dparams,l1,l2) : c1_nu_dirac(dparams,l1,l2);
        }
        double c2_final(const derived_params& dparams, bool majorana, flavour l1, flavour l2) {
          return majorana ? c2_nu_major(dparams,l1,l2) : c2_nu_dirac(dparams,l1,l2);
        }
        double c3_final(const derived_params& dparams, bool majorana, flavour l1, flavour l2) {
          return majorana ? c3_nu_major(dparams,l1,l2) : c3_nu_dirac(dparams,l1,l2);
        }
        double c4_final(const derived_params& dparams, bool majorana, flavour l1, flavour l2) {
          return majorana ? c4_nu_major(dparams,l1,l2) : c4_nu_dirac(dparams,l1,l2);
        }
        double c5_final(const derived_params& dparams, bool majorana, flavour l1, flavour l2) {
          return majorana ? c5_nu_major(dparams,l1,l2) : c5_nu_dirac(dparams,l1,l2);
        }
        double c6_final(const derived_params& dparams, bool majorana, flavour l1, flavour l2) {
          return majorana ? c6_nu_major(dparams,l1,l2) : c6_nu_dirac(dparams,l1,l2);
        }

        c_struct c_all_final(const derived_params& dparams, 
            bool majorana, flavour l1, flavour l2) {
          return {
            c1_final(dparams,majorana, l1, l2), c2_final(dparams,majorana, l1, l2),
            c3_final(dparams,majorana, l1, l2), c4_final(dparams,majorana, l1, l2),
            c5_final(dparams,majorana, l1, l2), c6_final(dparams,majorana, l1, l2)
          };
        }

        using hnl_decay_width_map_t = std::map<decay_modes,double>;
        
        double sum_decay_widths(const hnl_decay_width_map_t& map) {
          return std::accumulate(map.begin(), map.end(), 0.,
              [](double total, auto& pair) {
                return total + pair.second;
              });
        }
        
        hnl_decay_width_map_t get_hnl_decay_widths(const derived_params& dparams) 
        {
          const bool is_majorana = dparams.raw_params.is_majorana;
          const double HNL_mass = dparams.raw_params.HNL_mass;
          const double muon_mass = dparams.muon_mass;
          const double elec_mass = dparams.elec_mass;
          const double pion_pm_mass = dparams.pion_pm_mass;
          const double pion_0_mass = dparams.pion_0_mass;
          const double f2_pion = std::pow(dparams.config.physical_params().pion_decay_constant,2);
          const double CKM_Vud2 = std::pow(dparams.config.physical_params().CKM_Vud,2);
          
          hnl_decay_width_map_t decay_rates;
          decay_rates[decay_modes::nu_nu_nu] = (is_majorana ? 1. : 0.5) * decay_rate_to_3nu(dparams);
#ifdef DEBUG
          std::cout << "Add hnl->nu_nu_nu r="<<decay_rates[decay_modes::nu_nu_nu]<<std::endl;
#endif

          if(HNL_mass > 2*elec_mass) {
            const double decay_rate = decay_rate_to_nu_2lep(dparams, elec_mass, elec_mass,
                c_all_final(dparams, is_majorana, flavour::e, flavour::e));
            decay_rates[decay_modes::e_e_nu] = decay_rate;
#ifdef DEBUG
            std::cout << "Add hnl->e_e_nu r="<<decay_rate<<std::endl;
#endif
          }
          if(HNL_mass > elec_mass + muon_mass) {
            {
              const double decay_rate1 = decay_rate_to_nu_2lep(dparams, elec_mass, muon_mass,
                  c_all_final(dparams, is_majorana, flavour::e, flavour::m));
              decay_rates[decay_modes::e_mu_nu] = decay_rate1;
#ifdef DEBUG
              std::cout << "Add hnl->e_mu_nu r="<<decay_rate1<<std::endl;
#endif
            }
            {
              const double decay_rate2 = decay_rate_to_nu_2lep(dparams, elec_mass, muon_mass,
                  c_all_final(dparams, is_majorana, flavour::m, flavour::e));
              decay_rates[decay_modes::mu_e_nu] = decay_rate2;
#ifdef DEBUG
              std::cout << "Add hnl->e_mu_nu r="<<decay_rate2<<std::endl;
#endif
            }
          }
          if(HNL_mass > pion_0_mass) {
            const double decay_rate = (is_majorana ? 1. : 0.5) * decay_rate_to_nu_pseudoscalar(dparams, pion_0_mass, f2_pion);
            decay_rates[decay_modes::pi0_nu] = decay_rate;
#ifdef DEBUG
            std::cout << "Add hnl->pi0_nu r="<<decay_rate<<std::endl;
#endif
          }
          if(HNL_mass > elec_mass+pion_pm_mass) {
            const double decay_rate = (is_majorana ? 2. : 1.) * decay_rate_to_lepton_pseudoscalar(dparams, 
                elec_mass, pion_pm_mass, f2_pion, CKM_Vud2, dparams.U2(flavour::e));
            decay_rates[decay_modes::e_pi] = decay_rate;
#ifdef DEBUG
            std::cout << "Add hnl->e_pi r="<<decay_rate<<std::endl;
#endif
          }
          if(HNL_mass > 2*muon_mass) {
            const double decay_rate = decay_rate_to_nu_2lep(dparams, muon_mass, muon_mass,
                c_all_final(dparams, is_majorana, flavour::m, flavour::m));
            decay_rates[decay_modes::mu_mu_nu] = decay_rate;
#ifdef DEBUG
            std::cout << "Add hnl->mu_mu_nu r="<<decay_rate<<std::endl;
#endif
          }
          if(HNL_mass > muon_mass+pion_pm_mass) {
            const double decay_rate = (is_majorana ? 2. : 1.) * decay_rate_to_lepton_pseudoscalar(dparams, 
                muon_mass, pion_pm_mass, f2_pion, CKM_Vud2, dparams.U2(flavour::m));
            decay_rates[decay_modes::mu_pi] = decay_rate;
#ifdef DEBUG
            std::cout << "Add hnl->mu_pi r="<<decay_rate<<std::endl;
#endif
          }

          return decay_rates;
        }

        using both_hnls_t = std::pair<dkgen::core::particle_definition,dkgen::core::particle_definition>;
        both_hnls_t make_hnl_definitions(const derived_params& dparams,
           const std::vector<decay_modes>& decay_modes_to_use) {
          
          const double HNL_mass = dparams.raw_params.HNL_mass;

          
          const int pion_pm_pdg = dparams.config.physical_params().find_particle("pion_pm").pdgcode;
          const int pion_0_pdg  = dparams.config.physical_params().find_particle("pion_0").pdgcode;
          const int elec_pdg    = dparams.config.physical_params().find_particle("elec").pdgcode;
          const int muon_pdg    = dparams.config.physical_params().find_particle("muon").pdgcode;
          const int nu_pdg      = dparams.config.physical_params().find_particle("nu_e").pdgcode;
          
          
          const double pion_pm_mass = dparams.pion_pm_mass;
          const double pion_0_mass  = dparams.pion_0_mass;
          const double elec_mass    = dparams.elec_mass;
          const double muon_mass    = dparams.muon_mass;
          
          const double hbar = dparams.config.physical_params().hbar;
          const bool is_majorana = dparams.raw_params.is_majorana;

          auto decay_rates = get_hnl_decay_widths(dparams);
          const double total_decay_rate = sum_decay_widths(decay_rates);
          // total decay rate is helicity-independent according to https://arxiv.org/abs/1905.00284
          const double HNL_lifetime = hbar / total_decay_rate;

          // total decay rates are helicity-independent (HNL_lifetime will apply to both helicities)
          
          dkgen::core::particle_definition HNL_poshel_info{
            HNL_pdg_poshel, HNL_mass, HNL_lifetime, is_majorana? self_conjugate : !self_conjugate
          };
          dkgen::core::particle_definition HNL_neghel_info{
            HNL_pdg_neghel, HNL_mass, HNL_lifetime, is_majorana? self_conjugate : !self_conjugate
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

          auto decay_mode_enabled = [decay_modes_to_use](auto dm) -> bool {
            return std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), dm) != decay_modes_to_use.end();
          };


          if(HNL_mass > 2*elec_mass && decay_mode_enabled(decay_modes::e_e_nu)) {
            if(decay_rates[decay_modes::e_e_nu] > 0.) {
              // double lep_minus_mass, double lep_plus_mass, double c1, double c2, double c3, double c4, double c5, double c6
              auto rw_fn_pos = make_diff_decay_rate_function_to_nu_2lep(dparams, elec_mass, elec_mass,
                  c_all_final(dparams, is_majorana, flavour::e, flavour::e), helicity_plus);
              auto rw_fn_neg = make_diff_decay_rate_function_to_nu_2lep(dparams, elec_mass, elec_mass,
                  c_all_final(dparams, is_majorana, flavour::e, flavour::e), helicity_minus);
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
              auto rw_fn_pos = make_diff_decay_rate_function_to_nu_2lep(dparams, elec_mass, muon_mass,
                  c_all_final(dparams, is_majorana, flavour::e, flavour::m), helicity_plus);
              auto rw_fn_neg = make_diff_decay_rate_function_to_nu_2lep(dparams, elec_mass, muon_mass,
                  c_all_final(dparams, is_majorana, flavour::e, flavour::m), helicity_minus);
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
              auto rw_fn_pos = make_diff_decay_rate_function_to_nu_2lep(dparams, muon_mass, elec_mass,
                  c_all_final(dparams, is_majorana, flavour::m, flavour::e), helicity_plus);
              auto rw_fn_neg = make_diff_decay_rate_function_to_nu_2lep(dparams, muon_mass, elec_mass,
                  c_all_final(dparams, is_majorana, flavour::m, flavour::e), helicity_minus);
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
              if(is_majorana) {
                HNL_poshel_info.add_decay({
                    decay_rates[decay_modes::pi0_nu]/total_decay_rate,
                    {{pion_0_pdg,final_state},{nu_pdg,final_state}}
                    });
                HNL_neghel_info.add_decay({
                    decay_rates[decay_modes::pi0_nu]/total_decay_rate,
                    {{pion_0_pdg,final_state},{-nu_pdg,final_state}}
                    });
              } else {
                auto rw_fn_pos = make_diff_decay_rate_function_to_lepton_pseudoscalar(dparams, 0., pion_0_mass, helicity_plus);
                auto rw_fn_neg = make_diff_decay_rate_function_to_lepton_pseudoscalar(dparams, 0., pion_0_mass, helicity_minus);
                // formula 3.22 of arXiv:1905.00284 uses θ_P, so have to have pion first in
                // decay list to match twobody_dalitz_function def'n (cos_theta_1)
                HNL_poshel_info.add_decay(
                    core::decay_mode{ decay_rates[decay_modes::pi0_nu]/total_decay_rate,
                    {{pion_0_pdg,final_state},{nu_pdg,final_state}}
                    }.set_twobody_dalitz_reweighter(core::twobody_dalitz_function{rw_fn_pos}));
                HNL_neghel_info.add_decay(
                    core::decay_mode{ decay_rates[decay_modes::pi0_nu]/total_decay_rate,
                    {{pion_0_pdg,final_state},{-nu_pdg,final_state}}
                    }.set_twobody_dalitz_reweighter(core::twobody_dalitz_function{rw_fn_neg}));
              }
            }
          }

          if(HNL_mass > elec_mass+pion_pm_mass && decay_mode_enabled(decay_modes::e_pi)) {
            if(decay_rates[decay_modes::e_pi] > 0.) {
              auto rw_fn_pos = make_diff_decay_rate_function_to_lepton_pseudoscalar(dparams, elec_mass, pion_pm_mass, helicity_plus);
              auto rw_fn_neg = make_diff_decay_rate_function_to_lepton_pseudoscalar(dparams, elec_mass, pion_pm_mass, helicity_minus);
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
              auto rw_fn_pos = make_diff_decay_rate_function_to_nu_2lep(dparams, muon_mass, muon_mass,
                  c_all_final(dparams, is_majorana, flavour::m, flavour::m), helicity_plus);
              auto rw_fn_neg = make_diff_decay_rate_function_to_nu_2lep(dparams, muon_mass, muon_mass,
                  c_all_final(dparams, is_majorana, flavour::m, flavour::m), helicity_minus);
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
              auto rw_fn_pos = make_diff_decay_rate_function_to_lepton_pseudoscalar(dparams, muon_mass, pion_pm_mass, helicity_plus);
              auto rw_fn_neg = make_diff_decay_rate_function_to_lepton_pseudoscalar(dparams, muon_mass, pion_pm_mass, helicity_minus);
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

          HNL_poshel_info.finalise_decay_table();
          HNL_neghel_info.finalise_decay_table();


          return {HNL_poshel_info,HNL_neghel_info};

        }

      }
    }
  }
}
#undef DEBUG

#endif // __dkgen_physics_detail_hnl_hpp__
