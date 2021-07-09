#ifndef __GELATO_core_decaying_particle_info_hpp__
#define __GELATO_core_decaying_particle_info_hpp__

#include <memory>
#include <vector>

#include "GELATO/core/vectors.hpp"

namespace GELATO {
  namespace core {
    class particle_info {
      public:
        particle_info() : pdg_code{0}, log_weight{0.} { }
        particle_info(int pdg, fourvector prod_pos, fourvector prod_mom, double log_wt = 0.) 
          : pdg_code{pdg}, prod_pos{std::move(prod_pos)},  momentum{std::move(prod_mom)}, log_weight{log_wt} {
          reset_decay_position();
        };
        particle_info(int pdg, fourvector prod_pos, fourvector dec_pos, fourvector prod_mom, double log_wt = 0.)
          : pdg_code{pdg}, prod_pos{std::move(prod_pos)},
          dec_pos{std::move(dec_pos)}, momentum{std::move(prod_mom)}, log_weight{log_wt}   { }

        int pdg() const { return pdg_code; };
        double get_log_weight() const { return log_weight; };

        const fourvector& production_position() const { return prod_pos; };
        const fourvector& production_momentum() const { return momentum; };
        const fourvector& decay_position() const { return dec_pos; };
        const fourvector& decay_momentum() const { return momentum; };

      private:
        int pdg_code;
        fourvector prod_pos;
        fourvector dec_pos;
        fourvector momentum;
        double log_weight;
        void reset_decay_position();
        friend class decaying_particle_info;
    };
    class decaying_particle_info {
      public:
        using child_t = std::unique_ptr<decaying_particle_info>;
        using child_vector_t = std::vector<child_t>;
        using decaying_particle_info_ptr = decaying_particle_info*;

        // final states are non-decaying (within this framework)
        // pre final states are states that only decay to final states
        // non final states are the rest
        enum class state_type { non_final, pre_final_state, final_state };

        // no default constructor
        decaying_particle_info() = delete;
        //decaying_particle_info() : part_info{}, parent{nullptr}, state{state_type::final_state} { }
        // promote particle_info to decaying_particle_info as parent
        decaying_particle_info(particle_info pi, state_type state = state_type::non_final)
          : part_info{std::move(pi)}, parent{nullptr}, state{state} { }
        // Parent should be nullptr or ptr to parent 
        decaying_particle_info(decaying_particle_info_ptr parent, int pdg, const fourvector& prod_pos,
            const fourvector& prod_mom, state_type state)
          : part_info{pdg, prod_pos, prod_mom}, parent{parent}, state{state} { }
        /*
        // constructor for initial parent particle
        decaying_particle_info(int pdg, fourvector prod_pos,
            fourvector dec_pos, fourvector prod_mom, state_type state);
        */

        // need custom destructor to efficiently delete children without stack overflows in recersive functions
        ~decaying_particle_info();
        // and because of that, have to specify these ones:
        decaying_particle_info(decaying_particle_info&&) = default;
        decaying_particle_info(const decaying_particle_info&) = delete;
        decaying_particle_info& operator=(decaying_particle_info&&) = default;
        decaying_particle_info& operator=(const decaying_particle_info&) = delete;

        int pdg() const { return part_info.pdg_code; };

        const fourvector& production_position() const { return part_info.prod_pos; };
        const fourvector& production_momentum() const { return part_info.momentum; };
        const fourvector& decay_position() const { return part_info.dec_pos; };
        const fourvector& decay_momentum() const { return part_info.momentum; };

        double decay_log_weight() const { return part_info.log_weight; };

        decaying_particle_info& set_decay_pos_from_tof(double tof, double speed_of_light);
        decaying_particle_info& set_production_position(const fourvector& pos);
        decaying_particle_info& set_production_position(fourvector&& pos);

        decaying_particle_info& set_decay_log_weight(double log_wt) { part_info.log_weight = log_wt; return *this; };
        decaying_particle_info& multiply_decay_log_weight(double log_wt) { part_info.log_weight += log_wt; return *this; };

        decaying_particle_info* get_parent() const { return parent; }
        const child_vector_t& get_children() const { return children; }

        decaying_particle_info& add_child(child_t&& child) { children.push_back(std::move(child)); return *this; }

        bool is_final_state() const { return state == state_type::final_state; }
        // all daughters are final states
        bool is_pre_final_state() const { return state == state_type::pre_final_state; }
        decaying_particle_info& set_pre_final_state() { state = state_type::pre_final_state; return *this; }
        decaying_particle_info& unset_pre_final_state() { state = state_type::non_final; return *this; }

        bool is_decay_position_set() const { return !(part_info.dec_pos.t() < part_info.prod_pos.t()); }

        size_t get_number_of_particles_in_hierarchy() const;
        
      private:
        particle_info part_info;
        decaying_particle_info_ptr parent; // not owned by this object. deleting "this" has no effect on parent
        child_vector_t children; // owned by this object. deleting "this" deletes all children
        state_type state;
        
    };
  }
}

#endif
