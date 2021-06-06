#ifndef __dkgen_core_particle_hpp__
#define __dkgen_core_particle_hpp__

#include <vector>
#include <functional>

#include "dkgen/core/random_number_generators.hpp"

namespace dkgen {
  namespace core {

    class dalitz_function {
      public:
        using inv_mass_1_2_squared = double; // invariant mass of 1st and 2nd particle system
        using inv_mass_1_3_squared = double; // invariant mass of 1st and 3rd partcile system
        using functype = std::function<double(inv_mass_1_2_squared,inv_mass_1_3_squared)>;

        dalitz_function() : enabled(false) {}
        dalitz_function(functype fn) : enabled(true), /*maxweight(0.),*/ func(fn) {}
        
        double operator()(inv_mass_1_2_squared m12, inv_mass_1_3_squared m13) const {
          if(!enabled) return 1.;
          double ret = func(m12, m13);
          return ret;
        }
        bool is_enabled() const { return enabled; }
      private:
        bool enabled;
        functype func;
    };

    struct decay_mode {
      using pdg_code = int;
      using is_final_state = bool;
      using daughter_pair_t = std::pair<pdg_code,is_final_state>;
      using daughter_vector_t = std::vector<daughter_pair_t>;

      double branching_ratio;
      daughter_vector_t daughters;
      dalitz_function threebody_dalitz_reweighter;

      decay_mode() : branching_ratio(-1.) {};
      decay_mode(double br, daughter_vector_t dgt, dalitz_function rw = dalitz_function());

      decay_mode& set_daughters(const std::vector<std::pair<pdg_code,is_final_state>>& dgt) {
        daughters = dgt;
        return *this;
      }
      decay_mode& set_daughters(std::vector<std::pair<pdg_code,is_final_state>>&& dgt) {
        daughters = std::move(dgt);
        return *this;
      }

      decay_mode& set_threebody_dalitz_reweighter(const dalitz_function& rw) {
        threebody_dalitz_reweighter = rw;
        return *this;
      }
      decay_mode& set_threebody_dalitz_reweighter(dalitz_function&& rw) {
        threebody_dalitz_reweighter = std::move(rw);
        return *this;
      }
      bool is_pure_final_state() const;
      bool is_null() const { return branching_ratio < 0.; }
    };

    class particle_definition {
      public:
        particle_definition() = default;
        particle_definition(int pdgc, double m, double lt = -1.); // negative lifetime for stable or final state particles

        int pdg() const { return pdg_code; };
        double mass() const { return m; };
        double lifetime() const { return ltime; };

        particle_definition& add_decay(const decay_mode& dm);
        particle_definition& add_decay(decay_mode&& dm);

        const decay_mode& generate_decay_mode(random_uniform_0_1_generator rng) const;
        const decay_mode& generate_weighted_decay_mode(random_uniform_0_1_generator rng, double& weight) const;

        // sums up branching ratios
        particle_definition& finalise_decay_table(); 
      private:
        int pdg_code;
        double m;
        double ltime;
        std::vector<decay_mode> decay_table;
        std::vector<double> sum_branching_ratios;

        static decay_mode null_decay;
    };
  }
}

#endif
