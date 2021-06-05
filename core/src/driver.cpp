#include "dkgen/core/driver.hpp"

#include <queue>
#include <memory>
#include <cmath>
#include <numeric>

#include "dkgen/core/decaying_particle_info.hpp"
#include "dkgen/core/particle.hpp"


#ifdef DEBUG
#include <iostream>
#endif

dkgen::core::driver::driver() {
}

dkgen::core::driver::~driver() {
}

dkgen::core::particle_history dkgen::core::driver::generate_decays(
    dkgen::core::decaying_particle_info&& parent_meson,
    dkgen::core::random_uniform_0_1_generator rng
    ) const {

  particle_history ret;

  // will keep ownership of all children
  std::unique_ptr<decaying_particle_info> top_parent = std::make_unique<decaying_particle_info>(std::move(parent_meson));

  std::queue<decaying_particle_info_ptr> queue_to_decay;
  std::vector<decaying_particle_info_ptr> pure_final_states;

  // mutli stage process, will have to refactor into multiple functions
  // 
  // stage 1: generate decay chain of particles
  // establish how many truly-final states there are
  // randomly choose one truly-final state (A) to place inside detector
  //
  // stage 2: generate momenta of all particles
  //
  // stage 3: generate positions of all particles
  // including forcing (A) into the detector and applying weights
  //

  // stage 1: produce decay chain, and momenta (not positions)
  queue_to_decay.push(top_parent.get());
  while(!queue_to_decay.empty()) {
    auto decpar = queue_to_decay.front();
    queue_to_decay.pop();

    if(decpar->is_final_state()) continue;

    const particle_definition& parent_particle_info = find_particle(decpar->pdg());
    auto const& dm = parent_particle_info.generate_decay_mode(rng);
    if(dm.is_null()) continue;

    if(config.force_decays_inside_detector() && dm.is_pure_final_state()) {
      decpar->set_pre_final_state();
      pure_final_states.push_back(decpar);
    }

    struct daughter_info {
      int pdg;
      fourvector momentum;
      bool final_state;
    };
    std::vector<daughter_info> daughters;
    daughters.reserve(dm.daughters.size());

    auto twobody_decay_momentum_in_com_frame = [](double parent_mass, double m1, double m2) -> double {
      return std::sqrt((parent_mass+m1+m2)*(parent_mass-m1-m2)*(parent_mass+m1-m2)*(parent_mass-m1+m2))/(2.*parent_mass);
    };

    // need to flip signs of daughters if parent decay was antiparticle
    const int sign_flip = (decpar->pdg() != std::abs(decpar->pdg())) ? -1 : 1;

    if(dm.daughters.size() == 2) {
      auto const& d1 = find_particle(dm.daughters[0].first);
      auto const& d2 = find_particle(dm.daughters[1].first);
      const double m1 = d1.mass();
      const double m2 = d2.mass();
      const double cos_theta_cm = 1.-2.*rng();
      const double phi_cm = 2*M_PI*rng();
      const double decay_momentum = twobody_decay_momentum_in_com_frame(parent_particle_info.mass(), m1, m2);
      fourvector p1, p2;
#ifdef EXPOSE_PHYSICS_VECTORS
      p1.setVectM({
          decay_momentum*std::sin(std::acos(cos_theta_cm))*std::cos(phi_cm),
          decay_momentum*std::sin(std::acos(cos_theta_cm))*std::sin(phi_cm),
          decay_momentum*cos_theta_cm
          },m1);
      p2.setVectM({
          -decay_momentum*std::sin(std::acos(cos_theta_cm))*std::cos(phi_cm),
          -decay_momentum*std::sin(std::acos(cos_theta_cm))*std::sin(phi_cm),
          -decay_momentum*cos_theta_cm
          },m1);
      p1.boost(decpar->decay_momentum().boostVector());
      p2.boost(decpar->decay_momentum().boostVector());
#else
      p1.set_mass_momentum_theta_phi(m1,decay_momentum,std::acos(cos_theta_cm),phi_cm);
      p2.set_mass_momentum_theta_phi(m2,-decay_momentum,std::acos(cos_theta_cm),phi_cm);
      p1.boost(decpar->decay_momentum().get_boost_vector());
      p2.boost(decpar->decay_momentum().get_boost_vector());
#endif

      daughters.push_back({sign_flip * dm.daughters[0].first, p1, dm.daughters[0].second});
      daughters.push_back({sign_flip * dm.daughters[1].first, p2, dm.daughters[1].second});

    }
    else if(dm.daughters.size() == 3) {
      auto const& d1 = find_particle(dm.daughters[0].first);
      auto const& d2 = find_particle(dm.daughters[1].first);
      auto const& d3 = find_particle(dm.daughters[2].first);

      const double m1 = d1.mass();
      const double m2 = d2.mass();
      const double m3 = d3.mass();
      const double M  = parent_particle_info.mass();
      const double tcm = M - m1 - m2 - m3;

      // actually larger than the max weight
      const double decay_max_weight =
        twobody_decay_momentum_in_com_frame(tcm+m1,0,m1)
        * twobody_decay_momentum_in_com_frame(tcm+m1+m2,m1,m2);

      fourvector p1, p2, p3;
      while(true) {
        const double three_body_split = rng();
        const double m23 = m2 + m3 + three_body_split * tcm;

        const double decay_1_momentum = twobody_decay_momentum_in_com_frame(M, m23, m1);
        const double decay_23_momentum = twobody_decay_momentum_in_com_frame(m23, m2, m3);

        const double decay_weight = decay_1_momentum*decay_23_momentum;

        if(decay_weight < rng() * decay_max_weight) {
          continue; // while(true)
        }

        const double cos_theta_cm = 1.-2.*rng();
        const double phi_cm = 2*M_PI*rng();
        const double cos_theta_23 = 1.-2.*rng();
        const double phi_23 = 2*M_PI*rng();


#ifdef EXPOSE_PHYSICS_VECTORS
        p1.setVectM({
            decay_1_momentum*std::sin(std::acos(cos_theta_cm))*std::cos(phi_cm),
            decay_1_momentum*std::sin(std::acos(cos_theta_cm))*std::sin(phi_cm),
            decay_1_momentum*cos_theta_cm
            },m1);
        fourvector _p;
        _p.setVectM({
            -decay_1_momentum*std::sin(std::acos(cos_theta_cm))*std::cos(phi_cm),
            -decay_1_momentum*std::sin(std::acos(cos_theta_cm))*std::sin(phi_cm),
            -decay_1_momentum*cos_theta_cm
            },m23);
        const fourvector p23 = _p;

        p2.setVectM({
            decay_23_momentum*std::sin(std::acos(cos_theta_23))*std::cos(phi_23),
            decay_23_momentum*std::sin(std::acos(cos_theta_23))*std::sin(phi_23),
            decay_23_momentum*cos_theta_23
            },m2);
        p3.setVectM({
            -decay_23_momentum*std::sin(std::acos(cos_theta_23))*std::cos(phi_23),
            -decay_23_momentum*std::sin(std::acos(cos_theta_23))*std::sin(phi_23),
            -decay_23_momentum*cos_theta_23
            },m3);
        p2.boost(p23.boostVector());
        p3.boost(p23.boostVector());
#else
        p1.set_mass_momentum_theta_phi(m1,decay_1_momentum,std::acos(cos_theta_cm),phi_cm);
        const fourvector p23 = fourvector{}.set_mass_momentum_theta_phi(m23,-decay_1_momentum,std::acos(cos_theta_cm),phi_cm);

        p2.set_mass_momentum_theta_phi(m2,decay_23_momentum,std::acos(cos_theta_23),phi_23);
        p3.set_mass_momentum_theta_phi(m3,-decay_23_momentum,std::acos(cos_theta_23),phi_23);
        p2.boost(p23.get_boost_vector());
        p3.boost(p23.get_boost_vector());
#endif

        if(dm.threebody_dalitz_reweighter.is_enabled()) {
          const double invmass2_12 = (p1+p2).m2();
          const double invmass2_13 = (p1+p3).m2();
          const double dalitz_weight = dm.threebody_dalitz_reweighter(invmass2_12,invmass2_13);
          decpar->multiply_decay_weight(dalitz_weight);
        }
        break;
      }

#ifdef EXPOSE_PHYSICS_VECTORS
      p1.boost(decpar->decay_momentum().boostVector());
      p2.boost(decpar->decay_momentum().boostVector());
      p3.boost(decpar->decay_momentum().boostVector());
#else
      p1.boost(decpar->decay_momentum().get_boost_vector());
      p2.boost(decpar->decay_momentum().get_boost_vector());
      p3.boost(decpar->decay_momentum().get_boost_vector());
#endif

      daughters.push_back({sign_flip * dm.daughters[0].first, p1, dm.daughters[0].second});
      daughters.push_back({sign_flip * dm.daughters[1].first, p2, dm.daughters[1].second});
      daughters.push_back({sign_flip * dm.daughters[2].first, p3, dm.daughters[2].second});

    }
    else {
      throw std::runtime_error("4+ body decays are not implemented");
    }


    for(auto& d : daughters) {
      const int pdg = d.pdg;
      const fourvector& momentum = d.momentum;
      const bool is_final_state = d.final_state;
      decpar->add_child(std::make_unique<decaying_particle_info>(
            decpar, pdg, 
            decpar->decay_position(), momentum,
            is_final_state ? decaying_particle_info::final_state : decaying_particle_info::non_final
            ));
      queue_to_decay.push(decpar->get_children().back().get());
    }

  }


  // stage (1|2).5: choose random pure final state to force that decay into
  // the detector (if option has been enabled; otherwise pure_final_states is empty) 
  if(pure_final_states.size() > 1) {
    const size_t choice = static_cast<size_t>(pure_final_states.size() * rng());
    auto const& keep_this_one = *std::next(pure_final_states.begin(),choice);
    for(auto const& decpar : pure_final_states) {
      if(decpar != keep_this_one) {
        decpar->unset_pre_final_state();
      }
    }
  }

#ifdef DEBUG
  std::cout << "\n\n";
#endif
  // stage 3: generate decay positions
  const bool has_activity_inside_detector = generate_decay_position(top_parent.get(), rng);

  if(has_activity_inside_detector) {
    ret.build_hierarchy(std::move(top_parent));
  }

  return ret;
}

dkgen::core::driver& dkgen::core::driver::add_particle_definition(const dkgen::core::particle_definition& p) {
  particle_content[std::abs(p.pdg())] = p;
  return *this;
}

dkgen::core::driver& dkgen::core::driver::set_geometry(const dkgen::core::geometry& geom) {
  geo = geom;
  return *this;
}


dkgen::core::driver& dkgen::core::driver::set_particle_content(const particle_map& p) {
  particle_content = p;
  return *this;
}
dkgen::core::driver& dkgen::core::driver::set_particle_content(particle_map&& p) {
  particle_content = std::move(p);
  return *this;
}


bool dkgen::core::driver::generate_decay_position(decaying_particle_info_ptr parent,
    random_uniform_0_1_generator rng) const {
  if(!parent->is_decay_set() && !parent->is_final_state()) {

    const particle_definition& p = find_particle(parent->pdg());

    bool decay_set = false;
    if(parent->is_pre_final_state()) {
      // we have to force the decay position to be inside the detector, if possible

      const vector3& direction = parent->production_momentum().vect().unit();
      const vector3& origin = parent->production_position().vect();
      auto detector_points = geo.get_active_volume_intersections_for_beamline_vectors(origin, direction);
      if(detector_points.size() >= 2) {
        const double speed = parent->production_momentum().beta() * config.physical_params().speed_of_light;
        const double tof1 = (detector_points.first() - origin).mag()/speed;
        const double tof2 = (detector_points.last() - origin).mag()/speed;
        const double u = rng();
        // conversion between 0--1 to somewhere along tof1--tof2 with exponential decay
        const double tof = tof1 - p.lifetime() * std::log(1.-u+u*std::exp((tof1-tof2)/p.lifetime()));
        const double weight = std::exp(-tof1/p.lifetime()) - std::exp(-tof2/p.lifetime());
        parent->set_decay_pos_from_tof(tof, config.physical_params().speed_of_light);
        parent->set_decay_weight(weight);
        decay_set = true;
      }
      else {
        // particle cannot decay inside detector
        //  - generate decay without forcing, inside the if(!decay_set) block
      }

    }
    if(!decay_set) {
#ifdef DEBUG
      std::cout << "Setting non-forced decay\n";
#endif
      const double tof = -p.lifetime() * std::log(1. - rng());
      parent->set_decay_pos_from_tof(tof, config.physical_params().speed_of_light);
    }
    for(auto& d : parent->get_children()) { // also need to update the children positions
      d->set_production_position(parent->decay_position());
    }
  }

#ifdef DEBUG
  std::cout << "particle with PDG "<<parent->pdg() << " final "<<parent->is_final_state()<<" pre-final "<<parent->is_pre_final_state()<<'\n';
  std::cout << " - momentum "<<parent->decay_momentum().x()<<","<<parent->decay_momentum().y()<<","<<parent->decay_momentum().z()<<","<<parent->decay_momentum().t()<<"," << '\n';
  std::cout << " - position "<<parent->decay_position().x()<<","<<parent->decay_position().y()<<","<<parent->decay_position().z()<<","<<parent->decay_position().t()<<"," << '\n';
#endif
  // decay position is also set to production position for final states, so should still work for them
  const bool is_inside_detector = parent->is_final_state()
    && geo.is_beamline_vector_in_active_volume(parent->decay_position().vect());

  return std::accumulate(parent->get_children().begin(), parent->get_children().end(), is_inside_detector,
      [this,rng](bool previous_result, const decaying_particle_info::child_t& current_daughter) -> bool {
      return previous_result | generate_decay_position(current_daughter.get(),rng);
      }
      );
}

const dkgen::core::particle_definition& dkgen::core::driver::find_particle(int pdg) const {
  auto p = particle_content.find(std::abs(pdg));
  if (p == particle_content.end()) {
    throw std::runtime_error("Undefined particle requested! PDG code "+std::to_string(pdg));
  }
  return p->second;
};


dkgen::core::driver& dkgen::core::driver::set_config(const dkgen::core::config& conf) {
  if(conf.physical_params().speed_of_light < 0.) {
    throw std::runtime_error("config: system of units has not been fixed yet");
  }
  config = conf;
  return *this;
}
