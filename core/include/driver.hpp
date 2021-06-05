#ifndef __driver_hpp__
#define __driver_hpp__

#include "particle_history.hpp"
#include "random_number_generators.hpp"
#include "geometry.hpp"
#include "config.hpp"

#include <map>

namespace decaygen {
  class decaying_particle_info;
  class particle;
  class driver {
    public:
      using abs_particle_pdg = unsigned int;
      using particle_map = std::map<abs_particle_pdg, particle>;
      driver();
      ~driver();
      particle_history generate_decays(decaying_particle_info&& parent_meson, random_uniform_0_1_generator rng) const;
      driver& add_particle_definition(const particle& p);
      driver& set_particle_content(const particle_map& p);
      driver& set_particle_content(particle_map&& p);
      driver& set_geometry(const geometry& geom);
      driver& set_config(const config& conf);
    private:
      particle_map particle_content;
      geometry geo;
      config config;
      
      // returns true if particle or any of its daughters decays inside detector
      // forced_decay: decaying particle that is forced to decay inside detector (could be nullptr)
      bool generate_decay_position(decaying_particle_info_ptr parent,
          const decaying_particle_info_ptr forced_decay,
          random_uniform_0_1_generator rng) const;

      const particle& find_particle(int pdg) const;

  };
}

#endif
