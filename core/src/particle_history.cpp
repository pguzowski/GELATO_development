#include "dkgen/core/particle_history.hpp"

#include <cstdio> // for std::snprintf

#include "dkgen/core/decaying_particle_info.hpp"


dkgen::core::particle_history& dkgen::core::particle_history::build_hierarchy(std::unique_ptr<decaying_particle_info>&& p) {
  // first make sure p is not a descendent of current parent,
  // otherwise it will be deleted when current parent is reassigned and destructed
  decaying_particle_info_ptr curr = p.get();
  if(curr == parent.get()) return *this;
  while(curr != nullptr) {
    if(curr == parent.get()) {
      break;
    }
    curr = curr->get_parent();
  }
  if(curr != nullptr) {
    // not implemented
    throw std::runtime_error("particle_history::build_hierarchy: trying to add daughter of current parent. This shouldn't happen with unique_ptrs");
  }
  parent = std::move(p);
  if(parent) {
    event_counter++;
    hierarchy.reserve(parent->get_number_of_particles_in_hierarchy());
    hierarchy.push_back(parent.get());
    total_weight *= parent->decay_weight();
    std::function<void(decaying_particle_info_ptr)> traverse = [this,&traverse](decaying_particle_info_ptr current_parent) -> void {
      auto const& daughters = current_parent->get_children();

      for(auto const& d : daughters) {
        hierarchy.push_back(d.get());
        total_weight *= d->decay_weight();
      }
      // double-traversal because we have to flatten the hierarchy horizontally for hepevt format
      for(auto const& d : daughters) {
        traverse(d.get());
      }
    };
    traverse(parent.get());
  }
  return *this;
}

dkgen::core::hepevt_info dkgen::core::particle_history::build_hepevt_output() const {
  return build_hepevt_output_with_translation_followed_by_rotation(dkgen::core::fourvector{},dkgen::core::rotation{});
}


dkgen::core::hepevt_info dkgen::core::particle_history::build_hepevt_output_with_translation_followed_by_rotation(
    const dkgen::core::fourvector& translation, const dkgen::core::rotation& rot) const {
  if(hierarchy.empty()) {
    throw std::runtime_error("particle_history: tried to build hepevt output for empty hierarchy");
  }
  dkgen::core::hepevt_info ret;
  ret.event_counter = event_counter;
  ret.total_weight = total_weight;
  for(auto const& p : hierarchy) {

    auto find_pid = [this](auto p) -> int {
      auto pitr = std::find(hierarchy.begin(), hierarchy.end(), p);
      if (pitr == hierarchy.end()) {
        throw std::runtime_error("Undefined pointer found in particle_history::hierarchy");
      }
      return std::distance(hierarchy.begin(),pitr)+1; // need +1 for hepevt format
    };

    const bool has_parent = p->get_parent() != nullptr;
    const int mother = has_parent ? find_pid(p->get_parent()) : 0;
    
    const bool has_daughters = p->get_children().size() > 0;
    const int first_daughter = has_daughters ? find_pid(p->get_children().front().get()) : 0;
    const int last_daughter = has_daughters ? find_pid(p->get_children().back().get()) : 0;
    
    const dkgen::core::vector3& momentum_new = rot * p->production_momentum().vect();
    const dkgen::core::vector3& position_new = rot * (p->production_position().vect() + translation.vect());

    ret.particle_info.push_back({});
    hepevt_particle& part = ret.particle_info.back();
    part.status = has_parent ? (has_daughters ? 2 : 1) : 0;
    part.pdg = p->pdg();
    part.first_mother = mother;
    part.last_mother = mother;
    part.first_daughter = first_daughter;
    part.last_daughter = last_daughter;
    part.mass = p->production_momentum().m();
#ifdef EXPOSE_PHYSICS_VECTORS
    part.momentum.set(momentum_new, p->production_momentum().e());
    part.position.set(position_new, p->production_position().t() + translation.t());
#else
    part.momentum.set_vec_time(momentum_new, p->production_momentum().e());
    part.position.set_vec_time(position_new, p->production_position().t() + translation.t());
#endif

  }
  return ret;



}

std::string dkgen::core::hepevt_info::build_text(const char* metadata) const {
  const size_t SIZEOF_BUFFER = 1000;
  char buffer[SIZEOF_BUFFER+1]; // +1 for extra terminating null character
  std::snprintf(buffer,SIZEOF_BUFFER,"%d %lu weight=%g %s\n",event_counter,particle_info.size(), total_weight, metadata);
  std::string ret = buffer;
  for (auto const& p: particle_info) {
    std::snprintf(buffer,SIZEOF_BUFFER,
        "%d %d %d %d %d %d %f %f %f %f %f %f %f %f %f\n",
        p.status, p.pdg, p.first_mother, p.last_mother, p.first_daughter, p.last_daughter,
        p.momentum.px(), p.momentum.py(), p.momentum.pz(), p.momentum.e(), p.mass,
        p.position.x(), p.position.y(), p.position.z(), p.position.t()
        );
    ret += buffer;
  }
  return ret;
}

// initialize static member
unsigned int dkgen::core::particle_history::event_counter = 0;
