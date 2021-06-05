#ifndef __decaying_particle_info_hpp__
#define __decaying_particle_info_hpp__

#include <memory>
#include <vector>

#include "vectors.hpp"

namespace decaygen {
  class decaying_particle_info {
    public:
      using child_t = std::unique_ptr<decaying_particle_info>;
      using child_vector_t = std::vector<child_t>;
      using decaying_particle_info_ptr = decaying_particle_info*;

      // final states are non-decaying (within this framework)
      // pre final states are states that only decay to final states
      // non final states are the rest
      typedef enum { non_final, pre_final_state, final_state } state_type;
      
      decaying_particle_info();
      // Parent should be nullptr or ptr to parent (never new/delete)
      decaying_particle_info(decaying_particle_info* parent, int pdg, fourvector prod_pos,
          fourvector prod_mom, state_type state);
      // constructor for initial parent particle
      decaying_particle_info(int pdg, fourvector prod_pos,
          fourvector dec_pos, fourvector prod_mom, state_type state);
      ~decaying_particle_info();
      decaying_particle_info(decaying_particle_info&&) = default;
      decaying_particle_info(const decaying_particle_info&) = delete;
      decaying_particle_info& operator=(decaying_particle_info&&) = default;
      decaying_particle_info& operator=(const decaying_particle_info&) = delete;

      int pdg() const { return pdg_code; };

      const fourvector& production_position() const { return prod_pos; };
      const fourvector& production_momentum() const { return momentum; };
      const fourvector& decay_position() const { return dec_pos; };
      const fourvector& decay_momentum() const { return momentum; };

      double decay_weight() const { return decay_wt; };

      decaying_particle_info& set_decay_pos_from_tof(double tof, double speed_of_light);
      decaying_particle_info& set_production_position(const fourvector& pos);
      decaying_particle_info& set_production_position(fourvector&& pos);

      decaying_particle_info& set_decay_weight(double wt) { decay_wt = wt; return *this; };
      decaying_particle_info& multiply_decay_weight(double wt) { decay_wt *= wt; return *this; };
      
      decaying_particle_info* get_parent() const { return parent; }
      const child_vector_t& get_children() const { return children; }

      decaying_particle_info& add_child(child_t&& child) { children.push_back(std::move(child)); return *this; }
      // Release all children. Making this public allows to destroy children
      // backwards through hierarchy to avoid stack overflows in recursion
      // if there are too many layers of children
      void destroy_all_children() { children.clear(); }

      bool is_final_state() const { return state == final_state; }
      // all daughters are final states
      bool is_pre_final_state() const { return state == pre_final_state; }
      decaying_particle_info& set_pre_final_state() { state = pre_final_state; return *this; }
      decaying_particle_info& unset_pre_final_state() { state = non_final; return *this; }

      bool is_decay_set() const;
    private:
      int pdg_code;
      decaying_particle_info *parent; // not owned by this object. deleting this has no effect on parent
      child_vector_t children; // owned by this object. deleting this deletes all children
      fourvector prod_pos;
      fourvector dec_pos;
      fourvector momentum;
      double decay_wt;
      state_type state;

      void reset_decay_position();
  };
}

#endif
