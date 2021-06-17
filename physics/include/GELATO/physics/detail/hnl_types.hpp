#ifndef __GELATO_physics_detail_hnl_types_hpp__
#define __GELATO_physics_detail_hnl_types_hpp__

#include <cmath>

namespace GELATO {
  namespace physics {
    namespace detail {
      namespace hnl {

        // have to define separate helicity states, so that they can have different angular distributions
        constexpr int HNL_pdg_poshel  = 91; // free pdg code for 4th generation neutrino = 90 + helicity
        constexpr int HNL_pdg_neghel  = 89; // free pdg code for 4th generation neutrino = 90 + helicity
                                        // antipartciles have opposite helicity
                                        // -91: neghel; -89: poshel (still works with -90 + helicity scheme)
          
        constexpr bool self_conjugate = true;
        constexpr bool final_state = true;
        constexpr bool NOT_final_state = !final_state;

        constexpr int helicity_plus  = +1;
        constexpr int helicity_minus = -1;

        enum class flavour { e, m, t };
        enum class hnl_helicities { neg = -1, pos = 1 };
        enum class muon_production_modes { mu_e_N, mu_e_Nbar };

        enum class decay_modes { nu_nu_nu, e_e_nu, e_mu_nu, mu_e_nu, mu_mu_nu, e_pi, mu_pi, pi0_nu };
        enum class production_modes { k_mu2, k_mu3, k_e2, k_e3, k0_mu, k0_e, pi_mu, pi_e, mu_e };


        struct model_parameters {
          double HNL_mass;
          double U_e4, U_m4, U_t4;
          bool is_majorana;

          // have to specify these because of the other constructor
          model_parameters() = default;
          ~model_parameters() = default;
          model_parameters(const model_parameters&) = default;
          model_parameters(model_parameters&&) = default;
          model_parameters& operator=(const model_parameters&) = default;
          model_parameters& operator=(model_parameters&&) = default;

          // copy from GELATO::physics::heavy_neutral_leptons::model_parameters
          // without defining that struct in this header file or includes
          template<typename T> model_parameters(const T& p)
            : HNL_mass{p.HNL_mass}, U_e4{p.U_e4}, U_m4{p.U_m4}, U_t4{p.U_t4}, is_majorana{p.is_majorana} { }
        };
        // so we don't have to constantly recalculate these
        struct derived_params {
          const double U2_e;
          const double U2_m;
          const double U2_t;
          const double sum_U2;
          const double m5;
          const double m3;
          const double gFermi2;
          const double elec_mass;
          const double muon_mass;
          const double pion_pm_mass;
          const double pion_0_mass;
          const double kaon_pm_mass;
          const double kaon_0L_mass;
          const double gL;
          const double gR;
          const double pi3_inv;
          const model_parameters& raw_params;
          const GELATO::core::config& config;
          derived_params(const model_parameters& p,
              const GELATO::core::config& c) :
            U2_e{std::pow(p.U_e4,2)}, U2_m{std::pow(p.U_m4,2)},
            U2_t{std::pow(p.U_t4,2)}, sum_U2{std::pow(p.U_e4,2) + std::pow(p.U_m4,2) + std::pow(p.U_t4,2)},
            m5{std::pow(p.HNL_mass,5)}, m3{std::pow(p.HNL_mass,3)},
            gFermi2{std::pow(c.physical_params().gFermi,2)},
            elec_mass{c.physical_params().find_particle("elec").mass},
            muon_mass{c.physical_params().find_particle("muon").mass},
            pion_pm_mass{c.physical_params().find_particle("pion_pm").mass},
            pion_0_mass{c.physical_params().find_particle("pion_0").mass},
            kaon_pm_mass{c.physical_params().find_particle("kaon_pm").mass},
            kaon_0L_mass{c.physical_params().find_particle("kaon_0L").mass},
            gL{c.physical_params().sin2thW  - 0.5},
            gR{c.physical_params().sin2thW},
            pi3_inv{1./std::pow(M_PI,3)},
            raw_params(p), config(c) {      
            }
          double U2(flavour fl) const {
            switch(fl) {
              case flavour::e:
                return U2_e;
              case flavour::m:
                return U2_m;
              case flavour::t:
                return U2_t;
            }
            // should never reach here, but gcc warns about
            // control reaching end of non-void funciton
            return 0.;
          }
        };
      }
    }
  }
}
#endif
