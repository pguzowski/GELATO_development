#ifndef __dkgen_physics_heavy_neutral_leptons_hpp__
#define __dkgen_physics_heavy_neutral_leptons_hpp__

#include "dkgen/core/driver.hpp"
#include "dkgen/core/particle.hpp"

#include <cmath>

#include <gsl/gsl_integration.h>

#define DEBUG
#ifdef DEBUG
#include <iostream>
#endif

template< typename F >  class gsl_function_pp : public gsl_function {
  public:
    gsl_function_pp(const F& func) : _func(func) {
      function = &gsl_function_pp::invoke;
      params=this;
    }
  private:
    const F& _func;
    static double invoke(double x, void *params) {
      return static_cast<gsl_function_pp*>(params)->_func(x);
    }
};
namespace dkgen {
  namespace physics {
    namespace heavy_neutral_leptons {
      enum class decay_modes { nu_nu_nu, e_e_nu, e_mu_nu, mu_e_nu, mu_mu_nu, e_pi, mu_pi, pi0_nu };
      enum class production_modes { k_mu2, k_e2, k_e3, pi_mu, mu_e3 };
      struct model_parameters {
        double HNL_mass;
        double U_e4, U_m4, U_t4;
        bool is_majorana;
      };
      dkgen::core::driver::particle_map create_particle_content(
          const model_parameters& params,
          const dkgen::core::config& conf,
          // all "visible" final states
          std::vector<decay_modes> decay_modes_to_use = { 
            decay_modes::e_e_nu,
            decay_modes::e_mu_nu,
            decay_modes::mu_mu_nu,
            decay_modes::e_pi,
            decay_modes::mu_pi,
            decay_modes::pi0_nu,  
          },
          dkgen::core::driver::particle_map input = {}) {

        auto ret = input;

        auto& pion_pm = conf.physical_params().find_particle("pion_pm");
        auto& pion_0  = conf.physical_params().find_particle("pion_0");
        auto& kaon_pm = conf.physical_params().find_particle("kaon_pm");
        //auto& kaon_0L = conf.physical_params().find_particle("kaon_0L");
        auto& elec    = conf.physical_params().find_particle("elec");
        auto& muon    = conf.physical_params().find_particle("muon");
        auto& nu_e    = conf.physical_params().find_particle("nu_e");
        auto& nu_mu   = conf.physical_params().find_particle("nu_mu");
        auto& nu_tau  = conf.physical_params().find_particle("nu_tau");

        const int HNL_pdg  = 18; // free pdg code for 4th generation neutrino
        const int pion_pm_pdg = pion_pm.pdgcode;
        const int pion_0_pdg  = pion_0.pdgcode;
        const int kaon_pm_pdg = kaon_pm.pdgcode;
        //const int kaon_0L_pdg = kaon_0L.pdgcode;
        const int elec_pdg    = elec.pdgcode;
        const int muon_pdg    = muon.pdgcode;
        const int nu_e_pdg   = nu_e.pdgcode;
        const int nu_mu_pdg   = nu_mu.pdgcode;
        const int nu_tau_pdg   = nu_tau.pdgcode;

        // GeV system of units
        const double HNL_mass = params.HNL_mass;
        const double pion_pm_mass = pion_pm.mass;
        const double pion_0_mass  = pion_0.mass;
        const double kaon_pm_mass = kaon_pm.mass;
        //const double kaon_0L_mass = kaon_0L.mass;
        const double elec_mass    = elec.mass;
        const double muon_mass    = muon.mass;

        const double kaon_pm_lt = kaon_pm.lifetime;
        //const double kaon_0L_lt = kaon_0L.lifetime;
        const double pion_pm_lt = pion_pm.lifetime;
        const double pion_0_lt  = pion_0.lifetime;
        const double elec_lt    = elec.lifetime;
        const double muon_lt    = muon.lifetime;

        const double hbar = conf.physical_params().hbar;
        const double gfermi2 = std::pow(conf.physical_params().gFermi,2);

        // only implementing total decay rates for now;
        constexpr double pi3_inv = 1./M_PI/M_PI/M_PI; // 1/pi^3

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
        };
        const derived_params dparams(params);

        auto sqrt_kallen_lambda = [](double a, double b, double c)->double {
          return std::sqrt(a*a + b*b + c*c -2*a*b -2*a*c -2*b*c);
        };
        auto sqrtkl = sqrt_kallen_lambda;
        auto I1_2 = [sqrtkl](double x, double y) {
          return sqrtkl(1.,x,y)*(std::pow(1.-x,2) - y*(1.+x));
        };
        auto I1_2_1 = [sqrtkl](double x, double y, double cos_theta, int plus_minus) {
          const double sqrt_kl = sqrtkl(1.,x,y);
          return sqrt_kl / 4./M_PI * (std::pow(1.-x,2) - y*(1.+x) + plus_minus * (x-1.) * sqrt_kl * cos_theta);
        };
        auto integrate = [](auto fn, double low, double high) -> double { 
          double result, error;
          gsl_function_pp<decltype(fn)> Fp(fn);
          gsl_function *func = static_cast<gsl_function*>(&Fp);   
          size_t neval;
          gsl_integration_qng(func, low, high, 1e-5, 1e-5, &result, &error, &neval);
          return result;
        };
        auto I1_3 = [sqrtkl,integrate](double x, double y, double z) {
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
          
          return 12. * integrate(integrand, std::pow(std::sqrt(x)+std::sqrt(y),2), std::pow(1-std::sqrt(z),2));
        };
        auto I2_3 = [sqrtkl,integrate](double x, double y, double z) {
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
          return 24. * std::sqrt(y*z) * integrate(integrand, std::pow(std::sqrt(y)+std::sqrt(z),2), std::pow(1-std::sqrt(x),2));
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
        auto diff_decay_rate_to_lepton_pseudoscalar = [I1_2,I1_2_1,HNL_mass=dparams.raw_params.HNL_mass](
          double lep_mass, double pseudoscalar_mass, double cos_theta_l, int plus_minus) {
          const double xi_l = std::pow(lep_mass / HNL_mass,2);
          const double xi_p = std::pow(pseudoscalar_mass / HNL_mass,2);
          return I1_2_1(xi_l, xi_p, cos_theta_l, plus_minus)/I1_2(xi_l, xi_p);
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
        auto diff_decay_rate_to_nu_2lep = [sqrtkl,HNL_mass,I1_3,I2_3](double lep_minus_mass, double lep_plus_mass,
            double c1, double c2, double c3, double c4, double c5, double c6,
            double s1, double s2, double cos_theta_m, double cos_theta_p) {
          const double xi_m = std::pow(lep_minus_mass / HNL_mass,2);
          const double xi_p = std::pow(lep_plus_mass / HNL_mass,2);
          const double A02 = c1 * (s2 - xi_m)*(1+xi_p-s2) + c2 * (s1 - xi_p)*(1+xi_m-s1) + 2.*c3*xi_m*xi_p*(s1+s2-xi_m-xi_p);
          const double A12 = (c4 *(s2-xi_m) -2.*c6*xi_m*xi_p)*sqrtkl(1,s2,xi_p)*cos_theta_p
            + (c5*(s1-xi_p)-2.*c6*xi_m*xi_p)*sqrtkl(1.,s1,xi_m)*cos_theta_m;
          const double total_rate = (lep_minus_mass == lep_plus_mass)
            ?
            ((c1+c2) * I1_3(0,xi_m,xi_p) + c3 * I2_3(0,xi_m,xi_p)) 
            :
            (c1 * I1_3(0,xi_m,xi_p) + c2 * I1_3(0,xi_p,xi_m) + c3 * I2_3(0,xi_m,xi_p));
          return (A02+A12)/total_rate/12.;
        };

        const double gL = conf.physical_params().sin2thW  - 0.5;
        const double gR = conf.physical_params().sin2thW;
        enum class flavour { e, m, t };
        auto c1_nu_dirac = [&params,gL](flavour l1, flavour l2) { 
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U = (fl==flavour::e ? params.U_e4 : (fl==flavour::m? params.U_m4 : params.U_t4));
            //                                        vv  ----- l1 vs l2 in c2_nubar_dirac
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
            //                                        vv  ----- l2 vs l1 in c1_nu_dirac
            tot += U*U *((l1==l2 ? gL*gL : 0.) + (fl==l2?1:0)*(1.+(l1==l2?gL:0.)));
          }
          return tot;
        };
        auto c3_nu_dirac = [&params,gL,gR](flavour l1, flavour l2) { 
          if(l1 != l2) return 0.;
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U = (fl==flavour::e ? params.U_e4 : (fl==flavour::m? params.U_m4 : params.U_t4));
            //                                        vv  ----- l1 vs l2 in c2_nubar_dirac
            tot += U*U * ((fl==l2? 1. : 0.) + gL);
          }
          return tot*gR;
        };
        auto c3_nubar_dirac = [&params,gL,gR](flavour l1, flavour l2) { 
          if(l1 != l2) return 0.;
          double tot = 0.;
          for(auto fl: { flavour::e, flavour::m, flavour::t }) {
            const double U = (fl==flavour::e ? params.U_e4 : (fl==flavour::m? params.U_m4 : params.U_t4));
            //                                        vv  ----- l1 vs l2 in c2_nubar_dirac
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
        
        auto c1_nu_major = [c1_nu_dirac,c1_nubar_dirac](flavour l1, flavour l2){
          return c1_nu_dirac(l1,l2) + c1_nubar_dirac(l1,l2); };
        auto c2_nu_major = [c2_nu_dirac,c2_nubar_dirac](flavour l1, flavour l2){
          return c2_nu_dirac(l1,l2) + c2_nubar_dirac(l1,l2); };
        auto c3_nu_major = [c3_nu_dirac,c3_nubar_dirac](flavour l1, flavour l2){
          return c3_nu_dirac(l1,l2) + c3_nubar_dirac(l1,l2); };
        auto c4_nu_major = [c4_nu_dirac,c4_nubar_dirac](flavour l1, flavour l2){
          return c4_nu_dirac(l1,l2) - c4_nubar_dirac(l1,l2); };
        auto c5_nu_major = [c5_nu_dirac,c5_nubar_dirac](flavour l1, flavour l2){
          return c5_nu_dirac(l1,l2) - c5_nubar_dirac(l1,l2); };
        auto c6_nu_major = [](flavour , flavour ){ return 0.; };

        const double f2_pion = conf.physical_params().pion_decay_constant;
        const double CKM_Vud2 = std::pow(conf.physical_params().CKM_Vud,2);

        std::map<decay_modes, double> decay_rates;
        // always have 3nu mode
        double total_decay_rate = (params.is_majorana ? 1. : 0.5) * decay_rate_to_3nu;
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
        const double HNL_lifetime = hbar / total_decay_rate;

        const bool self_conjugate = true;
        const bool final_state = true;
        const bool NOT_final_state = !final_state;
        ret.push_back(dkgen::core::particle_definition{kaon_pm_pdg,kaon_pm_mass,kaon_pm_lt}
            .add_decay({1.,{{muon_pdg,final_state},{HNL_pdg,NOT_final_state}}})
            .finalise_decay_table());

        ret.push_back(dkgen::core::particle_definition{
            HNL_pdg, HNL_mass, HNL_lifetime, params.is_majorana? self_conjugate : !self_conjugate
            });
        auto& HNL_info = ret.back();

        if(HNL_mass > 2*elec_mass &&
            std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), decay_modes::e_e_nu) != decay_modes_to_use.end()) {
          HNL_info.add_decay(
              { decay_rates[decay_modes::e_e_nu]/total_decay_rate,
                {{nu_e_pdg,final_state},{elec_pdg,final_state},{-elec_pdg,final_state}}
              });
        }
        if(HNL_mass > elec_mass + muon_mass &&
            (std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), decay_modes::e_mu_nu) != decay_modes_to_use.end()
           || std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), decay_modes::mu_e_nu) != decay_modes_to_use.end())) {
          HNL_info.add_decay(
              { decay_rates[decay_modes::e_mu_nu]/total_decay_rate,
                {{nu_e_pdg,final_state},{elec_pdg,final_state},{-muon_pdg,final_state}}
              });
          HNL_info.add_decay(
              { decay_rates[decay_modes::mu_e_nu]/total_decay_rate,
                {{nu_e_pdg,final_state},{muon_pdg,final_state},{-elec_pdg,final_state}}
              });
        }
        if(HNL_mass > pion_0_mass &&
            std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), decay_modes::pi0_nu) != decay_modes_to_use.end()) {
          HNL_info.add_decay(
              { decay_rates[decay_modes::pi0_nu]/total_decay_rate, {{nu_e_pdg,final_state},{pion_0_pdg,final_state}} }
              );
        }
        if(HNL_mass > elec_mass+pion_pm_mass &&
            std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), decay_modes::e_pi) != decay_modes_to_use.end()) {
          const double majorana_factor = (params.is_majorana ? 0.5 : 1.);
          HNL_info.add_decay(
              { majorana_factor * decay_rates[decay_modes::e_pi]/total_decay_rate,
              {{elec_pdg,final_state},{pion_pm_pdg,final_state}}
              });
          if(params.is_majorana) {
            HNL_info.add_decay(
                { majorana_factor * decay_rates[decay_modes::e_pi]/total_decay_rate,
                {{-elec_pdg,final_state},{-pion_pm_pdg,final_state}}
                });
          }
        }
        if(HNL_mass > 2*muon_mass &&
            std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), decay_modes::mu_mu_nu) != decay_modes_to_use.end()) {
          HNL_info.add_decay(
              { decay_rates[decay_modes::mu_mu_nu]/total_decay_rate, {{nu_mu_pdg,final_state},{muon_pdg,final_state},{-muon_pdg,final_state}} }
              );
        }
        if(HNL_mass > muon_mass+pion_pm_mass &&
            std::find(decay_modes_to_use.begin(), decay_modes_to_use.end(), decay_modes::mu_pi) != decay_modes_to_use.end()) {
          const double majorana_factor = (params.is_majorana ? 0.5 : 1.);
          HNL_info.add_decay(
              { majorana_factor * decay_rates[decay_modes::mu_pi]/total_decay_rate,
                {{muon_pdg,final_state},{pion_pm_pdg,final_state}}
              });
          if(params.is_majorana) {
            HNL_info.add_decay(
                { majorana_factor * decay_rates[decay_modes::mu_pi]/total_decay_rate,
                  {{-muon_pdg,final_state},{-pion_pm_pdg,final_state}}
                });
          }
        }

        HNL_info.finalise_decay_table();
        
        ret.push_back(dkgen::core::particle_definition{pion_pm_pdg,pion_pm_mass,pion_pm_lt});
        ret.push_back(dkgen::core::particle_definition{pion_0_pdg,pion_0_mass,pion_0_lt,self_conjugate});
        ret.push_back(dkgen::core::particle_definition{elec_pdg,elec_mass,elec_lt});
        ret.push_back(dkgen::core::particle_definition{muon_pdg,muon_mass,muon_lt});
        ret.push_back(dkgen::core::particle_definition{nu_e_pdg,0.,-1.});
        ret.push_back(dkgen::core::particle_definition{nu_mu_pdg,0.,-1.});
        ret.push_back(dkgen::core::particle_definition{nu_tau_pdg,0.,-1.});

        return ret;

        // evade unused variable warnings in compiler
        (void)diff_decay_rate_to_lepton_pseudoscalar;
        (void)c4_nu_major;
        (void)c5_nu_major;
        (void)c6_nubar_dirac;
        (void)diff_decay_rate_to_nu_2lep;
        (void)c6_nu_dirac;
        (void)c6_nu_major;
        //if() {}
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
