#ifndef __particle_history_hpp__
#define __particle_history_hpp__

#include <memory>
#include <string>
#include <map>
#include <vector>

#include "vectors.hpp"

namespace decaygen {
  class decaying_particle_info;
  using decaying_particle_info_ptr = decaying_particle_info*;
  struct hepevt_particle {
    int status;
    int pdg;
    int first_mother, last_mother, first_daughter, last_daughter;
    fourvector momentum;
    double mass;
    fourvector position;
  };
  struct hepevt_info {
    unsigned int event_counter;
    double total_weight;
    std::vector<hepevt_particle> particle_info;
    std::string build_text() const;
  };
  class particle_history {
    public:
      particle_history();
      ~particle_history();
      particle_history(const particle_history&) = delete;
      particle_history(particle_history&&) = default;
      particle_history& operator=(const particle_history&) = delete;
      particle_history& operator=(particle_history&&) = default;

      double get_total_weight() const { return total_weight; };

      // transfer parent ownership to this object
      particle_history& build_hierarchy(std::unique_ptr<decaying_particle_info>&& p);

      hepevt_info build_hepevt_output() const;
      
      // allow translating states from beamline coordinates to final lab coordinates
      hepevt_info build_hepevt_output_with_translation_followed_by_rotation(
          const fourvector& translation, const rotation& rotation) const;

      explicit operator bool() const { return parent ? true : false; }

    private:
      std::unique_ptr<decaying_particle_info> parent; // need to keep ownership so that hierarchy remains valid
      std::map<int,decaying_particle_info_ptr> hierarchy;
      std::map<decaying_particle_info_ptr,int> particles_to_particle_ids;
      int last_particle_id;
      double total_weight;
      static unsigned int event_counter;
  };
}

#endif
