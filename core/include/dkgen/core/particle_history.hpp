#ifndef __dkgen_core_particle_history_hpp__
#define __dkgen_core_particle_history_hpp__

#include <memory>
#include <string>
#include <vector>

#include "dkgen/core/vectors.hpp"

namespace dkgen {
  namespace core {
    class decaying_particle_info;
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
      std::string build_text(const char* metadata = "") const;
    };
    class particle_history {
      public:
        particle_history()  : total_weight(1.) {}

        double get_total_weight() const { return total_weight; };

        // transfer parent ownership to this object
        particle_history& build_hierarchy(std::unique_ptr<decaying_particle_info>&& p);

        hepevt_info build_hepevt_output() const;

        // allow translating states from beamline coordinates to final lab coordinates
        hepevt_info build_hepevt_output_with_translation_followed_by_rotation(
            const fourvector& translation, const rotation& rotation) const;

        explicit operator bool() const { return parent ? true : false; }

      private:
        using decaying_particle_info_ptr = decaying_particle_info*;
        
        std::unique_ptr<decaying_particle_info> parent; // need to keep ownership so that hierarchy remains valid
        std::vector<decaying_particle_info_ptr> hierarchy;
        double total_weight;
        static unsigned int event_counter;
    };
  }
}

#endif
