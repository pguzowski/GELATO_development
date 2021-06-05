#include "particle_history.hpp"

#include <cstdio>

#include "decaying_particle_info.hpp"

/*
decaygen::particle_history& decaygen::particle_history::add_final_state(const std::shared_ptr<decaygen::decaying_particle_info>& p) {
  final_states[++last_daughter] = p;
  return *this;
}
decaygen::particle_history& decaygen::particle_history::add_final_state(std::shared_ptr<decaygen::decaying_particle_info>&& p) {
  final_states[++last_daughter] = std::move(p);
  return *this;
}

decaygen::particle_history& decaygen::particle_history::prepend_particle(const std::shared_ptr<decaygen::decaying_particle_info>& p) {
  total_weight *= p->decay_weight();
  parents[--last_parent] = p;
  return *this;
}
decaygen::particle_history& decaygen::particle_history::prepend_particle(std::shared_ptr<decaygen::decaying_particle_info>&& p) {
  total_weight *= p->decay_weight();
  parents[--last_parent] = std::move(p);
  return *this;
}

decaygen::particle_history& decaygen::particle_history::prepend_particle(std::shared_ptr<decaygen::decaying_particle_info>&& p) {
  total_weight *= p->decay_weight();
  parents[--last_parent] = std::move(p);
  return *this;
}
*/

decaygen::particle_history::particle_history() : total_weight(1.) {
}

decaygen::particle_history::~particle_history() {
}

decaygen::particle_history& decaygen::particle_history::build_hierarchy(std::unique_ptr<decaying_particle_info>&& p) {
  parent = std::move(p);
  if(parent) {
    event_counter++;
    last_particle_id=1;
    hierarchy[last_particle_id] = parent.get();
    particles_to_particle_ids[parent.get()] = last_particle_id;
    total_weight *= parent->decay_weight();
    std::function<void(decaying_particle_info_ptr)> traverse = [this,&traverse](decaying_particle_info_ptr current_parent) -> void {
      auto const& daughters = current_parent->get_children();

      for(auto const& d : daughters) {
        last_particle_id++;
        hierarchy[last_particle_id] = d.get();
        particles_to_particle_ids[d.get()] = last_particle_id;
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

decaygen::hepevt_info decaygen::particle_history::build_hepevt_output() const {
  return build_hepevt_output_with_translation_followed_by_rotation(decaygen::fourvector{},decaygen::rotation{});
}


decaygen::hepevt_info decaygen::particle_history::build_hepevt_output_with_translation_followed_by_rotation(
    const decaygen::fourvector& translation, const decaygen::rotation& rot) const {
  if(hierarchy.empty()) {
    throw std::runtime_error("particle_history: tried to build hepevt output for empty hierarchy");
  }
  decaygen::hepevt_info ret;
  ret.event_counter = event_counter;
  ret.total_weight = total_weight;
  for(auto const& h : hierarchy) {
    auto p = h.second;
    

    auto find_pid = [this](auto p) {
      auto pid = particles_to_particle_ids.find(p);
      if (pid == particles_to_particle_ids.end()) {
        throw std::runtime_error("Undefined pointer found in particle_history::particles_to_particle_ids");
      }
      return pid->second;
    };

    const bool has_parent = p->get_parent() != nullptr;
    const int mother = has_parent ? find_pid(p->get_parent()) : 0;
    
    const bool has_daughters = p->get_children().size() > 0;
    const int first_daughter = has_daughters ? find_pid(p->get_children().front().get()) : 0;
    const int last_daughter = has_daughters ? find_pid(p->get_children().back().get()) : 0;
    
    const decaygen::vector3& momentum_new = rot * p->production_momentum().vect();
    const decaygen::vector3& position_new = rot * (p->production_position().vect() + translation.vect());

    ret.particle_info.push_back({});
    hepevt_particle& part = ret.particle_info.back();
    part.status = has_parent ? (has_daughters ? 2 : 1) : 0;
    part.pdg = p->pdg();
    part.first_mother = mother;
    part.last_mother = mother;
    part.first_daughter = first_daughter;
    part.last_daughter = last_daughter;
    part.momentum.set_vec_time(momentum_new, p->production_momentum().e());
    part.mass = p->production_momentum().m();
    part.position.set_vec_time(position_new, p->production_position().t() + translation.t());

  }
  return ret;



}

std::string decaygen::hepevt_info::build_text() const {
  const size_t SIZEOF_BUFFER = 1000;
  char buffer[SIZEOF_BUFFER+1]; // +1 for extra terminating null character
  std::snprintf(buffer,SIZEOF_BUFFER,"%d %lu %f\n",event_counter,particle_info.size(), total_weight);
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

unsigned int decaygen::particle_history::event_counter = 0;
