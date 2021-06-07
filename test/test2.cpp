#include "dkgen/core/driver.hpp"
#include "dkgen/core/decaying_particle_info.hpp"

#include "dkgen/physics/higgs_portal_scalar.hpp"

#include <iostream>
#include <sstream>
#include <cmath>
#include <random>

int main(int argc, char** argv) {
  size_t n_to_gen = 10;
  double mass = 0.140, theta = 1e-3;
  bool force = false;
  for(int i = 0; i < argc /* as -n option will have argument */; ++i) {
    if(i+1 < argc && std::string(argv[i]) == "-n") {
      n_to_gen = std::atoll(argv[i+1]);
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-m") {
      mass = std::atof(argv[i+1]);
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-t") {
      theta = std::atof(argv[i+1]);
      i++;
      continue;
    }
    if(std::string(argv[i]) == "-f") {
      force = true;
      continue;
    }
  }

  dkgen::core::driver driver;

  dkgen::core::geometry geo({0,0,1.e4},{1e3,1e3,1e3});
  driver.set_geometry(geo);
  dkgen::core::config conf;
  conf.fix_system_of_units(dkgen::core::config::system_of_units::GeV_cm_ns);
  conf.set_force_decays_in_detector(force);
  driver.set_config(conf);

  const double scalar_mass = mass; // GeV
  const double scalar_theta = theta;
  const std::string metadata{(std::ostringstream() << "model_theta=" << scalar_theta).str()};
  //const std::string metadata = [](){ auto s = std::istringstream(); s << "model_theta=" << scalar_theta; return s.str(); }();
  driver.set_particle_content(dkgen::physics::higgs_portal_from_kaons::create_particle_content(scalar_mass, scalar_theta,conf));

  auto& kaon = conf.physical_params().find_particle("kaon_pm");
  const int kpdg = -kaon.pdgcode;
  const double kmass = kaon.mass;
  const double kmom = 0.5;
  
  std::uniform_real_distribution<double> rng;
  std::default_random_engine gen;
  //std::vector<double> moms(n_to_gen);
  size_t i = 0;
  while(i++ < n_to_gen) {
    auto const& res = driver.generate_decays(
        {
          kpdg,
          {0.,0.,0.,0.}, // production position
          {0.,0.,0.,0.}, // decay position
          {0.,0.,kmom,std::sqrt(kmom*kmom + kmass*kmass)}, // momentum
          dkgen::core::decaying_particle_info::state_type::non_final // state
        },
        [&rng, &gen]()->double{return rng(gen);});

    if(res) {
      std::cout << res.build_hepevt_output().build_text(metadata.c_str()) << std::endl;
    }
  }

  return 0;
}
