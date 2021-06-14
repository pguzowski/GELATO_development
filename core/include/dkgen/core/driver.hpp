#ifndef __dkgen_core_driver_hpp__
#define __dkgen_core_driver_hpp__

#include "dkgen/core/particle_history.hpp"
#include "dkgen/core/random_number_generators.hpp"
#include "dkgen/core/geometry.hpp"
#include "dkgen/core/config.hpp"

#include <map>

namespace dkgen {
  namespace core {
    class particle_info;
    class particle_definition;
    class driver {
      public:
        using particle_map = std::vector<particle_definition>;
        
        driver& add_particle_definition(const particle_definition& p);
        driver& set_particle_content(const particle_map& p);
        driver& set_particle_content(particle_map&& p);
        driver& set_geometry(const geometry& geom);
        driver& set_config(const config& conf);
        
        // main event loop: generate decay of initial particle
        particle_history generate_decays(const particle_info& initial_decay, random_uniform_0_1_generator rng) const;
        
      private:
        using decaying_particle_info_ptr = decaying_particle_info*;
        
        particle_map particle_content;
        geometry geo;
        ::dkgen::core::config config;

        // recursively generate decay positions, forcing some into the detector if
        // possible and "force_decays_inside_detector" config flag is on
        // returns true if any final state particles are produced inside detector
        bool generate_decay_position(decaying_particle_info_ptr parent,
            random_uniform_0_1_generator rng) const;
        /*
        bool generate_decay_position(decaying_particle_info_ptr parent,
            const decaying_particle_info_ptr forced_decay,
            random_uniform_0_1_generator rng) const;
        */

        void sort_particles();
        const particle_definition& find_particle(int pdg) const;

        // recursively generates daughters from decay tables, and momenta
        // doesn't do positions (use generate_decay_position(...) for that)
        void generate_decay_tree(decaying_particle_info_ptr parent,
            std::vector<decaying_particle_info_ptr>& queue,
            std::vector<decaying_particle_info_ptr>& pure_final_states,
            random_uniform_0_1_generator rng) const;
    };
  }
}

#endif
