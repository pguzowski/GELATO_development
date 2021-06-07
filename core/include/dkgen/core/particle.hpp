#ifndef __dkgen_core_particle_hpp__
#define __dkgen_core_particle_hpp__

#include <vector>
#include <functional>

#include "dkgen/core/random_number_generators.hpp"

namespace dkgen {
  namespace core {

    // for fermions, can have angular distribution
    class twobody_dalitz_function {
      public:
        using cos_theta = double; // angle w.r.t. parent momentum (or parent's parent momentum if decaying at rest)
        using functype = std::function<double(cos_theta cth)>;

        twobody_dalitz_function() : enabled(false) {}
        twobody_dalitz_function(functype fn) : enabled(true),  func(fn) {}
        
        double operator()(cos_theta cth) const {
          if(!enabled) return 1.;
          return func(cth);
        }
        bool is_enabled() const { return enabled; }
      private:
        bool enabled;
        functype func;
    };

    class threebody_dalitz_function {
      public:
        using reduced_inv_mass_1_2_squared = double; // invariant mass of 1st and 2nd particle system
        using reduced_inv_mass_1_3_squared = double; // invariant mass of 1st and 3rd partcile system
        using functype = std::function<double(reduced_inv_mass_1_2_squared,reduced_inv_mass_1_3_squared)>;

        threebody_dalitz_function() : enabled(false) {}
        threebody_dalitz_function(functype fn) : enabled(true), /*maxweight(0.),*/ func(fn) {}
        
        double operator()(reduced_inv_mass_1_2_squared m12, reduced_inv_mass_1_3_squared m13) const {
          if(!enabled) return 1.;
          double ret = func(m12, m13);
          return ret;
        }
        bool is_enabled() const { return enabled; }
      private:
        bool enabled;
        functype func;
    };

    class angular_threebody_dalitz_function {
      public:
        using cos_theta_2 = double; // angle w.r.t. parent momentum (or parent's parent momentum if decaying at rest)
        using phi_2 = double; // azimuthal angle of 2nd particle w.r.t. parent
        using cos_theta_3 = double; // angle w.r.t. parent momentum (or parent's parent momentum if decaying at rest)
        using phi_23 = double; // azimuthal angle between 2nd and 3rd particle
        using reduced_inv_mass_1_2_squared = double; // invariant mass of 1st and 2nd particle system
        using reduced_inv_mass_1_3_squared = double; // invariant mass of 1st and 3rd partcile system
        using functype = std::function<double(reduced_inv_mass_1_2_squared,reduced_inv_mass_1_3_squared,
            cos_theta_2,phi_2,cos_theta_3,phi_23)>;

        angular_threebody_dalitz_function() : enabled(false) {}
        angular_threebody_dalitz_function(functype fn) : enabled(true), /*maxweight(0.),*/ func(fn) {}
        
        double operator()(reduced_inv_mass_1_2_squared m12, reduced_inv_mass_1_3_squared m13,
            cos_theta_2 cth2, phi_2 ph2, cos_theta_3 cth3, phi_23 ph23) const {
          if(!enabled) return 1.;
          double ret = func(m12, m13, cth2, ph2, cth3, ph23);
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
      twobody_dalitz_function twobody_dalitz_reweighter;
      threebody_dalitz_function threebody_dalitz_reweighter;
      angular_threebody_dalitz_function angular_threebody_dalitz_reweighter;

      enum class reweighter_type {
        none,
        twobody,
        threebody,
        threebody_angular
      };
      reweighter_type reweighting_type;

      decay_mode() : branching_ratio(-1.), reweighting_type(reweighter_type::none) {};
      decay_mode(double br, daughter_vector_t dgt);
      

      decay_mode& set_daughters(const std::vector<std::pair<pdg_code,is_final_state>>& dgt) {
        daughters = dgt;
        return *this;
      }
      decay_mode& set_daughters(std::vector<std::pair<pdg_code,is_final_state>>&& dgt) {
        daughters = std::move(dgt);
        return *this;
      }

      decay_mode& set_twobody_dalitz_reweighter(const twobody_dalitz_function& rw) {
        twobody_dalitz_reweighter = rw;
        reweighting_type = reweighter_type::twobody;
        return *this;
      }
      decay_mode& set_twobody_dalitz_reweighter(twobody_dalitz_function&& rw) {
        twobody_dalitz_reweighter = std::move(rw);
        reweighting_type = reweighter_type::twobody;
        return *this;
      }
      decay_mode& set_threebody_dalitz_reweighter(const threebody_dalitz_function& rw) {
        threebody_dalitz_reweighter = rw;
        reweighting_type = reweighter_type::threebody;
        return *this;
      }
      decay_mode& set_threebody_dalitz_reweighter(threebody_dalitz_function&& rw) {
        threebody_dalitz_reweighter = std::move(rw);
        reweighting_type = reweighter_type::threebody;
        return *this;
      }
      decay_mode& set_angular_dalitz_reweighter(const angular_threebody_dalitz_function& rw) {
        angular_threebody_dalitz_reweighter = rw;
        reweighting_type = reweighter_type::threebody_angular;
        return *this;
      }
      decay_mode& set_threebody_dalitz_reweighter(angular_threebody_dalitz_function&& rw) {
        angular_threebody_dalitz_reweighter = std::move(rw);
        reweighting_type = reweighter_type::threebody_angular;
        return *this;
      }
      bool is_pure_final_state() const;
      bool is_null() const { return branching_ratio < 0.; }
    };

    class particle_definition {
      public:
        particle_definition() = default;
        // negative lifetime for stable or final state particles
        particle_definition(int pdgc, double m, double lt = -1., bool scj = false);

        int pdg() const { return pdg_code; };
        int antipdg() const { return self_conjugate ? pdg_code : -pdg_code; };
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
        bool self_conjugate;

        static decay_mode null_decay;
    };
  }
}

#endif
