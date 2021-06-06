#include "dkgen/core/driver.hpp"
#include "dkgen/core/particle.hpp"
#include "dkgen/core/decaying_particle_info.hpp"

#include <iostream>
#include <cmath>
#include <random>

int main(int argc, char** argv) {
  size_t n_to_gen = 10;
  bool to_reweight = false;
  for(int i = 0; i < argc; ++i) {
    // number of primary decays to simulate
    if(i+1 < argc && std::string(argv[i]) == "-n") {
      n_to_gen = std::atoll(argv[i+1]);
      i++;
    }
    // use three body reweighting
    else if(std::string(argv[i]) == "-r") {
      to_reweight = true;
    }
  }

  dkgen::core::driver driver;

  dkgen::core::geometry geo({0,1e4,1.e4},{1e3,1e3,1e3});
  driver.set_geometry(geo);
  dkgen::core::config conf;
  conf.fix_system_of_units(dkgen::core::config::system_of_units::GeV_cm_ns);
  conf.set_force_decays_in_detector(true);
  driver.set_config(conf);

  driver.add_particle_definition(
      dkgen::core::particle_definition{1,100.,1e3}
      .add_decay({1.,{{-2,false},{-3,false}}})
      .finalise_decay_table()
      );
  
  driver.add_particle_definition(
      dkgen::core::particle_definition{-2,80.,.4e3}
      .add_decay({1.,{{11,true},{-11,true}}})
      .finalise_decay_table()
      );
  

  if(to_reweight) {
    // toy reweighter
    dkgen::core::dalitz_function rw{[](dkgen::core::dalitz_function::inv_mass_1_2_squared m12,
        dkgen::core::dalitz_function::inv_mass_1_3_squared m13) -> double {
      return 1.-std::abs((m12-m13)/(m12+m13)); // maximise symmetric invariant masses
    }};
    
    driver.add_particle_definition(
        dkgen::core::particle_definition{3,15.,.2e3}
        .add_decay({1.,{{13,true},{-13,true},{11,true}},rw})
        .finalise_decay_table()
        );
  }
  else {
    driver.add_particle_definition(
        dkgen::core::particle_definition{3,15.,.2e3}
        .add_decay({1.,{{13,true},{-13,true},{-11,true}}}) // no reweight
        .finalise_decay_table()
        );
  }

  driver.add_particle_definition({11,511e-6,-1.});
  driver.add_particle_definition({13,0.109,-1.});
  
  std::uniform_real_distribution<double> rng;
  std::default_random_engine gen;
  //std::vector<double> moms(n_to_gen);
  const double mom = 10.;
  //for(auto mom : moms) {
  size_t i = 0;
  while(i++ < n_to_gen) {
    auto const& res = driver.generate_decays(
        {
        1,
        {0.,0.,0.,0.},
        {0.,0.,0.,0.},
        {0.,0.,mom,std::sqrt(100.*100.+mom*mom)},
        dkgen::core::decaying_particle_info::state_type::non_final
        },
        [&rng, &gen]()->double{return rng(gen);});

    if(res) {
      std::cout << res.build_hepevt_output().build_text() << std::endl;
    }
  }

  return 0;
}
