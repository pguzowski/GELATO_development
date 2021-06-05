#include "decaying_particle_info.hpp"

#include <algorithm>

decaygen::decaying_particle_info::decaying_particle_info() 
  : pdg_code{0}, parent{nullptr}, decay_wt{1.}, state{final_state}
{
}

decaygen::decaying_particle_info::decaying_particle_info(decaying_particle_info* parent, int pdg,
    fourvector prodpos, fourvector prod_mom, state_type state) 
  : pdg_code{pdg}, parent{parent}, prod_pos{std::move(prodpos)},
  momentum{std::move(prod_mom)}, decay_wt{1.}, state{state}
{
  reset_decay_position();
}

//Release all children. Work backwards through hierarchy to avoid
// stack overflows in recursion if there are too many layers of children
decaygen::decaying_particle_info::~decaying_particle_info() {
  while(!children.empty()) {
    decaying_particle_info* p = this;

    auto has_nonempty_children = [](auto& child) -> bool {
      return child->get_children().size() > 0;
    };

    // returns true if any children have non-empty children of their own.
    // loop breaks when all children of p are final states, so we can safely clear p's children
    // without calling further destructors on the stack
    while(std::any_of(p->children.begin(), p->children.end(), has_nonempty_children)) {
      // finds the first child that has children of its own
      p = std::find_if(p->children.begin(), p->children.end(), has_nonempty_children)->get();
    }
    p->children.clear();
  }
}


decaygen::decaying_particle_info&
decaygen::decaying_particle_info::set_production_position(const fourvector& prodpos) {
  prod_pos = prodpos;
  reset_decay_position();
  return *this;
}

decaygen::decaying_particle_info&
decaygen::decaying_particle_info::set_production_position(fourvector&& prodpos) {
  prod_pos = std::move(prodpos);
  reset_decay_position();
  return *this;
}

void decaygen::decaying_particle_info::reset_decay_position() {
  dec_pos = prod_pos;
  // set decay time to less than production time, to invalidate
  if(dec_pos.t() >= 0.) {
    dec_pos.set_t(-1.);
  }
  else {
    // deal with some floating point rounding conditions;
    dec_pos.set_t(prod_pos.t() - std::abs(prod_pos.t()) * 1e-10);
    if(!(dec_pos.t() < prod_pos.t())) {
      // invalid decay time couldn't be set, probably due to the
      // production time being close to the largest negative floating point
      throw std::runtime_error("particle production time is too large and negative");
    } 
  }
}


bool decaygen::decaying_particle_info::is_decay_set() const {
  return !(dec_pos.t() < prod_pos.t());
}

decaygen::decaying_particle_info::decaying_particle_info(int pdg,
          fourvector prod_pos, fourvector dec_pos, fourvector prod_mom, state_type state) 
  : pdg_code{pdg}, parent{nullptr}, prod_pos{std::move(prod_pos)}, dec_pos{std::move(dec_pos)},
  momentum{std::move(prod_mom)}, decay_wt{1.}, state{state}
{
}

decaygen::decaying_particle_info& decaygen::decaying_particle_info::set_decay_pos_from_tof(double tof, double speed_of_light) {
  const double dist = momentum.beta() * tof * speed_of_light;
  const vector3& direction = momentum.vect().unit();
  const vector3& decpos = prod_pos.vect() + direction * dist;
  dec_pos.set_vec_time(decpos, prod_pos.t() + tof);
  return *this;
}

