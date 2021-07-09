#include "GELATO/core/driver.hpp"

#include <queue>
#include <memory>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "GELATO/core/decaying_particle_info.hpp"
#include "GELATO/core/particle.hpp"


//#define DEBUG
#ifdef DEBUG
#include <iostream>
#endif

GELATO::core::particle_history GELATO::core::driver::generate_decays(
    const GELATO::core::particle_info& parent_meson,
    GELATO::core::random_uniform_0_1_generator rng
    ) const {

  particle_history ret;

  // will keep ownership of all children
  std::unique_ptr<decaying_particle_info> top_parent = std::make_unique<decaying_particle_info>(parent_meson);

  // implement queue as std::vector that grows
  // (not too worried about memory size, but worried about timing)
  static std::vector<decaying_particle_info_ptr> queue_to_decay(50);
  queue_to_decay.clear();
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
  queue_to_decay.push_back(top_parent.get());
  auto qiter = queue_to_decay.begin();
  size_t n_decayed = 0; // keeps track of how many decays have been done, to advance qiter
  //while(!queue_to_decay.empty()) {
  while(qiter != queue_to_decay.end()) {
    generate_decay_tree(*qiter, queue_to_decay, pure_final_states, rng);
    // cannot just do ++qiter because queue_to_decay might have reallocated
    qiter = std::next(queue_to_decay.begin(), ++n_decayed);
  }


  // stage 1.5: choose random pure final state to force that decay into
  // the detector (if option has been enabled; otherwise pure_final_states is empty) 
  if(!pure_final_states.empty()) {
    const size_t choice = pure_final_states.size() > 1 ? static_cast<size_t>(pure_final_states.size() * rng()) : 0;
    auto const& keep_this_one = *std::next(pure_final_states.begin(),choice);
    for(auto const& decpar : pure_final_states) {
      if(decpar != keep_this_one) {
        decpar->unset_pre_final_state();
      }
    }
    // cascade state_types up the chain from pre-final state
    decaying_particle_info::state_type state = keep_this_one->get_state_type() + 1;
    decaying_particle_info_ptr pparent = keep_this_one->get_parent();
    while(pparent != nullptr) {
      pparent->set_state_type(state);
      state++;
      pparent = pparent->get_parent();
    }
  }

#ifdef DEBUG
  std::cout << "\n\n";
#endif
  // stage 3: generate decay positions
  const bool has_activity_inside_detector = generate_decay_position(top_parent.get(), rng);

  std::function<double(decaying_particle_info_ptr)> accumulate_weights =
    [&accumulate_weights](decaying_particle_info_ptr particle) {
    return std::accumulate(particle->get_children().begin(), particle->get_children().end(), particle->decay_log_weight(),
        [&accumulate_weights](double a, auto& b){ return a + accumulate_weights(b.get()); });
  };
  const double event_log_weight =
    configuration.get_minimum_event_log_weight() > std::numeric_limits<double>::lowest() ?
    accumulate_weights(top_parent.get()) : 0.;

  if(has_activity_inside_detector && event_log_weight >= configuration.get_minimum_event_log_weight()) {
    ret.build_hierarchy(std::move(top_parent));
  }

  return ret;
}

GELATO::core::driver& GELATO::core::driver::add_particle_definition(const GELATO::core::particle_definition& p) {
  particle_content.push_back(p);
  sort_particles();
  return *this;
}

GELATO::core::driver& GELATO::core::driver::set_geometry(const GELATO::core::geometry& geom) {
  geo = geom;
  return *this;
}


GELATO::core::driver& GELATO::core::driver::set_particle_content(const particle_map& p) {
  particle_content = p;
  sort_particles();
  return *this;
}
GELATO::core::driver& GELATO::core::driver::set_particle_content(particle_map&& p) {
  particle_content = std::move(p);
  sort_particles();
  return *this;
}

void GELATO::core::driver::sort_particles() {
  std::sort(particle_content.begin(), particle_content.end(),
      [](auto& a, auto& b) { return std::abs(a.pdg()) < std::abs(b.pdg()); });
}


bool GELATO::core::driver::generate_decay_position(decaying_particle_info_ptr parent,
    random_uniform_0_1_generator rng) const {
  if(!parent->is_decay_position_set() && !parent->is_final_state()) {

    const particle_definition& p = find_particle(parent->pdg());
    if(p.lifetime() < 0.) {
      return false; // doesn't decay inside detector because it doesn't decay
                    // also won't set decay positions of any children etc. because
                    // theoretically they will be produced at infinity as this doesn't decay
    }

    bool decay_set = false;
    if(p.lifetime() == 0.) {
      // decay straight away
      parent->set_decay_pos_from_tof(0., configuration.physical_params().speed_of_light);
      decay_set = true;
    }
    const double dilated_lifetime = p.lifetime() * parent->production_momentum().gamma();

    const decaying_particle_info::state_type state = parent->get_state_type();
    if(!decay_set && (state == decaying_particle_info::state_types::pre_final_state
          || state == decaying_particle_info::state_types::pre_final_state+1)) {
      // we have to force the decay position to be inside the detector, if possible

      const vector3& direction = parent->production_momentum().vect().mag() > 0. ?
        parent->production_momentum().vect().unit() : vector3{0.,0.,0.};
      const vector3& origin = parent->production_position().vect();
      auto find_pre_final_state_direction = [](decaying_particle_info_ptr p) {
        for(auto const& d : p->get_children()) {
          if(d->is_pre_final_state()) return d->production_momentum().vect().unit();
        }
        throw std::runtime_error("Cannot find daughter direction of pre-final state");
      };
      auto detector_points = 
        state == decaying_particle_info::state_types::pre_final_state
        ? geo.get_active_volume_intersections_for_beamline_vectors(origin, direction)
        : geo.project_detector_onto_first_vector_along_direction_of_second_vector(
            origin, direction, find_pre_final_state_direction(parent));
      if(detector_points.size() >= 2) {
        const double speed = parent->production_momentum().beta() * configuration.physical_params().speed_of_light;
        const double tof1 = (detector_points.first() - origin).mag()/speed;
        const double tof2 = (detector_points.last() - origin).mag()/speed;
        const double u = rng();
        // conversion between 0--1 to somewhere along tof1--tof2 with exponential decay
        // argument of std::log is always positive because 0<u<1, lifetime>0,tof1<tof2 => std::exp<1
        // FIXME: this code will introduce floating-point rounding errors if p.lifetime is too big
        //        might need a switch to a linear approximation
        //        [too big = ~10^300 sec, i.e. much longer than lifetime of unvierse ^ many powers]
        
        const double tof = tof1 - dilated_lifetime * std::log(1.- u + u*std::exp((tof1-tof2)/dilated_lifetime));

        // to get from production point to detector
        const double log_survival_prob = -tof1/dilated_lifetime;

        // probability of not decaying beyond detector
        // = 1 - integral(tof2, infinity)
        // conditional ::==> tof2-tof1
        const double log_prob_doesnt_decay_beyond_detector = std::log(1.-std::exp(-(tof2-tof1)/dilated_lifetime));

        // prob decays inside detector = prob survives to detector * prob doesn't decay beyond detector
        const double log_prob_decays_in_detector = log_survival_prob + log_prob_doesnt_decay_beyond_detector;

        parent->set_decay_pos_from_tof(tof, configuration.physical_params().speed_of_light);
        parent->set_decay_log_weight(log_prob_decays_in_detector);
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
      const double u = rng();
      const double tof = -dilated_lifetime * std::log(u);
      parent->set_decay_pos_from_tof(tof, configuration.physical_params().speed_of_light);
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

const GELATO::core::particle_definition& GELATO::core::driver::find_particle(int pdg) const {
  auto p = std::lower_bound(particle_content.begin(), particle_content.end(), std::abs(pdg),
      [](auto& a, int find_pdg) { return std::abs(a.pdg()) < find_pdg; });
  if (p == particle_content.end() || std::abs(p->pdg()) != std::abs(pdg)) {
    throw std::runtime_error("Undefined particle requested! PDG code "+std::to_string(pdg));
  }
  return *p;
}


GELATO::core::driver& GELATO::core::driver::set_config(const config& conf) {
  if(conf.physical_params().speed_of_light < 0.) {
    throw std::runtime_error("configuration: system of units has not been fixed yet");
  }
  configuration = conf;
  return *this;
}


void
GELATO::core::driver::generate_decay_tree(decaying_particle_info_ptr parent, std::vector<decaying_particle_info_ptr>& queue,
    std::vector<decaying_particle_info_ptr>& pure_final_states,
    random_uniform_0_1_generator rng) const {

  if(parent->is_final_state()) return;

  const particle_definition& parent_particle_info = find_particle(parent->pdg());
  if(parent_particle_info.lifetime() < 0.) return; // stable particle
  double weight;
  auto const& dm = parent_particle_info.generate_weighted_decay_mode(rng, weight);
#ifdef DEBUG
  if(std::abs(parent->pdg())==89 || std::abs(parent->pdg())==91) {
    std::cout << "for decay_mode size: "<<dm.daughters.size()<<" pdgs:";
    for(auto& dmdm: dm.daughters) std::cout << " "<<dmdm.first;
    std::cout<<std::endl;
  }
#endif
  
  parent->multiply_decay_log_weight(std::log(weight));
  if(dm.is_null()) return;

  if(configuration.force_decays_inside_detector() && dm.is_pure_final_state()) {
    parent->set_pre_final_state();
    pure_final_states.push_back(parent);
  }

  auto twobody_decay_momentum_in_com_frame = [](double parent_mass, double m1, double m2) -> double {
    return std::sqrt((parent_mass+m1+m2)*(parent_mass-m1-m2)*(parent_mass+m1-m2)*(parent_mass-m1+m2))/(2.*parent_mass);
  };

  // need to flip signs of daughters if parent decay was antiparticle
  const bool sign_flip = (parent->pdg() != parent_particle_info.pdg()) ? true : false;


  auto add_new_daughter_to_queue = [&parent, &queue](int pdgc, const fourvector& pos,
      fourvector&& mom, auto fs) {
    parent->add_child(std::make_unique<decaying_particle_info>(parent, pdgc, pos, std::move(mom), fs));
    queue.push_back(parent->get_children().back().get());
  };

  if(dm.daughters.size() == 2) {
    auto const& d1 = find_particle(dm.daughters[0].first);
    auto const& d2 = find_particle(dm.daughters[1].first);
    const bool sign_flip1 = (d1.pdg() != dm.daughters[0].first) ? true : false;
    const bool sign_flip2 = (d2.pdg() != dm.daughters[1].first) ? true : false;
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
        },m2);
#else
    p1.set_mass_momentum_theta_phi(m1,decay_momentum,std::acos(cos_theta_cm),phi_cm);
    p2.set_mass_momentum_theta_phi(m2,-decay_momentum,std::acos(cos_theta_cm),phi_cm);
#endif
    

    // need to do the angular reweighting before boosting
    switch(dm.reweighting_type) {
      case decay_mode::reweighter_type::twobody:
        {
          // use parent, unless parent is at rest
          // else use -ve parent of parent (pparent), unless pparent is also at rest
          // else use z-axis
          const vector3 reference = ((parent->decay_momentum().vect().mag() > 0)
              ? parent->decay_momentum().vect().unit()
              : ((parent->get_parent() != nullptr && parent->get_parent()->decay_momentum().vect().mag() > 0)
                ? -(parent->get_parent()->decay_momentum().vect().unit())
                : vector3{0.,0.,1.} // last resort, use z-axis
                )
              );
          const double costh = p1.vect().unit().dot(reference);
          const double dalitz_weight = dm.twobody_dalitz_reweighter(costh);
          parent->multiply_decay_log_weight(std::log(dalitz_weight));
        }
        break;
      default:
        // do nothing
        break;
    }

#ifdef EXPOSE_PHYSICS_VECTORS
    p1.boost(parent->decay_momentum().boostVector());
    p2.boost(parent->decay_momentum().boostVector());
#else
    p1.boost(parent->decay_momentum().get_boost_vector());
    p2.boost(parent->decay_momentum().get_boost_vector());
#endif

    add_new_daughter_to_queue(sign_flip ^ sign_flip1 ? d1.antipdg() : d1.pdg(), parent->decay_position(), std::move(p1),
          dm.daughters[0].second ?
          decaying_particle_info::state_types::final_state : decaying_particle_info::state_types::non_final);
    add_new_daughter_to_queue(sign_flip ^ sign_flip2 ? d2.antipdg() : d2.pdg(), parent->decay_position(), std::move(p2),
          dm.daughters[1].second ?
          decaying_particle_info::state_types::final_state : decaying_particle_info::state_types::non_final);

  }
  else if(dm.daughters.size() == 3) {
    auto const& d1 = find_particle(dm.daughters[0].first);
    auto const& d2 = find_particle(dm.daughters[1].first);
    auto const& d3 = find_particle(dm.daughters[2].first);
    const bool sign_flip1 = (d1.pdg() != dm.daughters[0].first) ? true : false;
    const bool sign_flip2 = (d2.pdg() != dm.daughters[1].first) ? true : false;
    const bool sign_flip3 = (d3.pdg() != dm.daughters[2].first) ? true : false;

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
      break; // break out of while(true)
    }

    // need to do the angular ones before boosting
    switch(dm.reweighting_type) {
      case decay_mode::reweighter_type::threebody:
        {
          const double red_invmass2_12 = (p1+p2).m2() / parent->decay_momentum().m2();
          const double red_invmass2_13 = (p1+p3).m2() / parent->decay_momentum().m2();
          const double dalitz_weight = dm.threebody_dalitz_reweighter(red_invmass2_12,red_invmass2_13);
          parent->multiply_decay_log_weight(std::log(dalitz_weight));
        }
        break;
      case decay_mode::reweighter_type::threebody_angular:
        {
          const double red_invmass2_12 = (p1+p2).m2() / parent->decay_momentum().m2();
          const double red_invmass2_13 = (p1+p3).m2() / parent->decay_momentum().m2();
          // use parent, unless parent is at rest
          // else use -ve parent of parent (pparent), unless pparent is also at rest
          // else use z-axis
          const vector3 reference = ((parent->decay_momentum().vect().mag() > 0)
              ? parent->decay_momentum().vect().unit()
              : ((parent->get_parent() != nullptr && parent->get_parent()->decay_momentum().vect().mag() > 0)
                ? -(parent->get_parent()->decay_momentum().vect().unit())
                : vector3{0.,0.,1.} // last resort, use z-axis
                )
              );
          const double costh_2 = p2.vect().unit().dot(reference);
          const double costh_3 = p3.vect().unit().dot(reference);
#ifdef EXPOSE_PHYSICS_VECTORS
          const double phi_2 = p2.vect().unit().deltaPhi(reference);
          const double phi_23 = p2.vect().unit().deltaPhi(p3.vect().unit());
#else
          const double phi_2 = p2.vect().unit().delta_phi(reference);
          const double phi_23 = p2.vect().unit().delta_phi(p3.vect().unit());
#endif
          const double dalitz_weight = dm.angular_threebody_dalitz_reweighter(
              red_invmass2_12,red_invmass2_13,costh_2,phi_2,costh_3,phi_23);
          parent->multiply_decay_log_weight(std::log(dalitz_weight));
        }
        break;
      default:
        // do nothing
        break;
    }

#ifdef EXPOSE_PHYSICS_VECTORS
    p1.boost(parent->decay_momentum().boostVector());
    p2.boost(parent->decay_momentum().boostVector());
    p3.boost(parent->decay_momentum().boostVector());
#else
    p1.boost(parent->decay_momentum().get_boost_vector());
    p2.boost(parent->decay_momentum().get_boost_vector());
    p3.boost(parent->decay_momentum().get_boost_vector());
#endif

    add_new_daughter_to_queue(sign_flip ^ sign_flip1 ? d1.antipdg() : d1.pdg(), parent->decay_position(), std::move(p1),
          dm.daughters[0].second ?
          decaying_particle_info::state_types::final_state : decaying_particle_info::state_types::non_final);
    add_new_daughter_to_queue(sign_flip ^ sign_flip2 ? d2.antipdg() : d2.pdg(), parent->decay_position(), std::move(p2),
          dm.daughters[1].second ?
          decaying_particle_info::state_types::final_state : decaying_particle_info::state_types::non_final);
    add_new_daughter_to_queue(sign_flip ^ sign_flip3 ? d3.antipdg() : d3.pdg(), parent->decay_position(), std::move(p3),
          dm.daughters[2].second ?
          decaying_particle_info::state_types::final_state : decaying_particle_info::state_types::non_final);
  }
  else {
    throw std::runtime_error("4+ body decays are not implemented");
  }

#ifdef DEBUG
if(std::abs(parent->pdg())==89 || std::abs(parent->pdg())==91) {
  std::cout << "put on parent size: "<<parent->get_children().size()<<" pdgs:";
  for(auto& dmdm: parent->get_children()) std::cout << " "<<dmdm->pdg();
  std::cout<<std::endl;
}
#endif

}
