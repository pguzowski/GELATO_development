#include "GELATO/core/decaying_particle_info.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>

void GELATO::core::particle_info::reset_decay_position() {
  dec_pos = prod_pos;
  // set decay time to less than production time, to invalidate
  if(dec_pos.t() >= 0.) {
#ifdef EXPOSE_PHYSICS_VECTORS
    dec_pos.setT(-1.);
#else
    dec_pos.set_t(-1.);
#endif
  }
  else {
    // deal with some floating point rounding conditions;
#ifdef EXPOSE_PHYSICS_VECTORS
    dec_pos.setT(prod_pos.t() - std::abs(prod_pos.t()) * 1e-10);
#else
    dec_pos.set_t(prod_pos.t() - std::abs(prod_pos.t()) * 1e-10);
#endif
    if(!(dec_pos.t() < prod_pos.t())) {
      // invalid decay time couldn't be set, probably due to the
      // production time being close to the largest negative floating point
      throw std::runtime_error("particle production time is too large and negative");
    } 
  }
}



////////// now decaying: ////////////////////////////////////////////////////////////////////////////

//Release all children. Work backwards through hierarchy to avoid
// stack overflows in recursion if there are too many layers of children
GELATO::core::decaying_particle_info::~decaying_particle_info() {
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


GELATO::core::decaying_particle_info&
GELATO::core::decaying_particle_info::set_production_position(const fourvector& prodpos) {
  part_info.prod_pos = prodpos;
  part_info.reset_decay_position();
  return *this;
}

GELATO::core::decaying_particle_info&
GELATO::core::decaying_particle_info::set_production_position(fourvector&& prodpos) {
  part_info.prod_pos = std::move(prodpos);
  part_info.reset_decay_position();
  return *this;
}


GELATO::core::decaying_particle_info& GELATO::core::decaying_particle_info::set_decay_pos_from_tof(double tof, double speed_of_light) {
  const double dist = part_info.momentum.beta() * tof * speed_of_light;
  const vector3& direction = part_info.momentum.vect().mag() > 0. ? part_info.momentum.vect().unit() : vector3{0.,0.,0.};
  const vector3& decpos = part_info.prod_pos.vect() + direction * dist;
#ifdef EXPOSE_PHYSICS_VECTORS
  part_info.dec_pos.set(decpos, part_info.prod_pos.t() + tof);
#else
  part_info.dec_pos.set_vec_time(decpos, part_info.prod_pos.t() + tof);
#endif
  return *this;
}

size_t GELATO::core::decaying_particle_info::get_number_of_particles_in_hierarchy() const {
  const size_t initial = 1; // 1 for this particle
  return std::accumulate(children.begin(), children.end(), initial,
      [](size_t s, auto& c) {
        return s + c->get_number_of_particles_in_hierarchy();
      });
}
