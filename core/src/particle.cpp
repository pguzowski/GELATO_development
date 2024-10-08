#include "GELATO/core/particle.hpp"

#include <algorithm>
#include <string>

GELATO::core::decay_mode::decay_mode(double br, daughter_vector_t dgt) :
  branching_ratio{br},
  daughters{std::move(dgt)},
  reweighting_type{reweighter_type::none}
{
}

bool GELATO::core::decay_mode::is_pure_final_state() const {
  return std::all_of(daughters.begin(), daughters.end(),[](auto p){return p.second;});
}

GELATO::core::particle_definition::particle_definition(int pdgc, double m, double lt, bool scj) :
  pdg_code{pdgc},
  m{m},
  ltime{lt},
  self_conjugate(scj)
{
  if(m < 0.) {
    throw std::runtime_error("Negative particle mass");
  }
}
      
GELATO::core::particle_definition& GELATO::core::particle_definition::add_decay(const decay_mode& dm) {
  if(sum_branching_ratios.size() > 0) {
    throw std::runtime_error("particle_definition::add_decay(const&): Decay table has been finalised already!");
  }
  if(dm.branching_ratio < 0.) {
    throw std::runtime_error("Negative branching ratio is unphysical");
  }
  decay_table.push_back(dm);
  return *this;
}
GELATO::core::particle_definition& GELATO::core::particle_definition::add_decay(decay_mode&& dm) {
  if(sum_branching_ratios.size() > 0) {
    throw std::runtime_error("particle_definition::add_decay(&&):Decay table has been finalised already!");
  }
  if(dm.branching_ratio < 0.) {
    throw std::runtime_error("Negative branching ratio is unphysical");
  }
  decay_table.push_back(std::move(dm));
  return *this;
}

GELATO::core::particle_definition& GELATO::core::particle_definition::finalise_decay_table() {
  if(sum_branching_ratios.size() > 0) {
    throw std::runtime_error("particle_definition::finalise_decay_table(): Decay table has been finalised already!");
  }
  // sorting by largest branching ratio speeds up code later on when randomly choosing modes
  std::sort(decay_table.begin(), decay_table.end(), [](auto& a, auto& b){ return a.branching_ratio > b.branching_ratio; });
  double rolling_sum_branching_ratios = 0.;
  for(auto const& dm : decay_table) {
    // also check that only 2 and 3 body decays are present
    if(dm.daughters.size() < 2 || dm.daughters.size() > 3) {
      throw std::runtime_error("only 2- and 3-body decays are currently implemented");
    }
    rolling_sum_branching_ratios += dm.branching_ratio;
    sum_branching_ratios.push_back(rolling_sum_branching_ratios);
  }
  return *this;
}

const GELATO::core::decay_mode& GELATO::core::particle_definition::generate_decay_mode(random_uniform_0_1_generator rng) const {
  if(decay_table.empty()) return null_decay;
  if(sum_branching_ratios.empty()) {
    throw std::runtime_error("Decay table has not been finalised");
  }
  if(sum_branching_ratios.back() > 1.0000001) {
    throw std::runtime_error("Branching ratio sum > 1.0000001 : "+std::to_string(sum_branching_ratios.back()-1.));
  }
  if(decay_table.size() == 1 && sum_branching_ratios.back() == 1.) {
    return decay_table.at(0);
  }
  const double x = (sum_branching_ratios.back() > 1. ? sum_branching_ratios.back() : 1.) * rng();
  for(size_t i = 0; i < decay_table.size(); ++i) {
    if(x < sum_branching_ratios.at(i)) {
      return decay_table.at(i);
    }
  }
  return null_decay;
}



const GELATO::core::decay_mode&
GELATO::core::particle_definition::generate_weighted_decay_mode(random_uniform_0_1_generator rng, double& weight) const {
  if(decay_table.empty()) {
    weight = 1.;
    return null_decay;
  }
  if(sum_branching_ratios.empty()) {
    throw std::runtime_error("Decay table has not been finalised");
  }
  weight = sum_branching_ratios.back();
  if(decay_table.size() == 1) {
    return decay_table.back();
  }
  const double x = rng() * weight;
  for(size_t i = 0; i < decay_table.size(); ++i) {
    if(x < sum_branching_ratios.at(i)) {
      return decay_table.at(i);
    }
  }
  // should in theory never get here.
  // But might, because of summing doubles with rounding error,
  // in which case return the last decay mode
  return decay_table.back();
}

GELATO::core::decay_mode GELATO::core::particle_definition::null_decay{};
