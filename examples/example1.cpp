#include "GELATO/core/driver.hpp"
#include "GELATO/core/particle.hpp"
#include "GELATO/core/decaying_particle_info.hpp"

#include <iostream>
#include <cmath>
#include <random>

/* generates a toy decay table with the structure:

    1 --- 3 --- 11
       |     |
       |      - 13
       |     |
       |      - (-13)
       |
        - 2 --- 11
             |
              - (-11)

  options:
    -n: number of initial decays to simulate
    -r: whether to use the 3-body reweighter for particle #3
  */

int main(int argc, char** argv) {
  size_t n_to_gen = 10;
  bool to_reweight = false;
  for(int i = 0; i < argc; ++i) {
    // number of primary decays to simulate
    if(i+1 < argc && std::string(argv[i]) == "-n") {
      n_to_gen = std::atoll(argv[i+1]);
      i++;
      continue;
    }
    // use three body reweighting
    if(std::string(argv[i]) == "-r") {
      to_reweight = true;
      continue;
    }
  }

  GELATO::core::driver driver;

  std::unique_ptr<GELATO::core::geometry> geo
    = std::make_unique<GELATO::core::geometry>(GELATO::core::vector3{0,1e4,1.e4},GELATO::core::vector3{1e3,1e3,1e3});
  driver.set_geometry(std::move(geo));
  GELATO::core::config conf;
  conf.fix_system_of_units(GELATO::core::config::system_of_units::GeV_cm_ns);
  conf.set_force_decays_in_detector(true);
  driver.set_config(conf);

  // add 1 -> 2 + 3
  driver.add_particle_definition(
      GELATO::core::particle_definition{1,100.,1e3}
      .add_decay({1.,{{-2,false},{-3,false}}})
      .finalise_decay_table()
      );
  
  // add 2 -> 11 + -11
  driver.add_particle_definition(
      GELATO::core::particle_definition{-2,80.,.4e3}
      .add_decay({1.,{{11,true},{-11,true}}})
      .finalise_decay_table()
      );
  

  // add 3 -> 11 + 13 + (-13), with optional reweighter
  if(to_reweight) {
    // toy reweighter
    GELATO::core::threebody_dalitz_function rw{[](
        GELATO::core::threebody_dalitz_function::reduced_inv_mass_1_2_squared m12,
        GELATO::core::threebody_dalitz_function::reduced_inv_mass_1_3_squared m13
        ) -> double {
      return 1.-std::abs((m12-m13)/(m12+m13)); // maximise symmetric invariant masses
    }};
    
    driver.add_particle_definition(
        GELATO::core::particle_definition{3,15.,.2e3}
        .add_decay(
          GELATO::core::decay_mode{1.,{{13,true},{-13,true},{11,true}}}.set_threebody_dalitz_reweighter(rw)
          )
        .finalise_decay_table()
        );
  }
  else {
    driver.add_particle_definition(
        GELATO::core::particle_definition{3,15.,.2e3}
        .add_decay({1.,{{13,true},{-13,true},{-11,true}}}) // no reweight
        .finalise_decay_table()
        );
  }

  // add stable 11, 13
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
        // initial particle with pdg code 1
        {
        1, // pdg code
        {0.,0.,0.,0.}, // production position
        {0.,0.,0.,0.}, // decay position
        {0.,0.,mom,std::sqrt(100.*100.+mom*mom)}, // decay momentum
        },
        // lambda for the RNG
        [&rng, &gen]()->double{return rng(gen);});

    if(res) {
      std::cout << res.build_hepevt_output().build_text() << std::endl;
    }
  }

  return 0;
}
