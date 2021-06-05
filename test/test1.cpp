#include "driver.hpp"
#include "particle.hpp"
#include "decaying_particle_info.hpp"

#include <iostream>
#include <cmath>
#include <random>

int main(int argc, char** argv) {
  size_t n_to_gen = 10;
  for(size_t i = 0; i < argc-1 /* as -n option will have argument */; ++i) {
    if(std::string(argv[i]) == "-n") {
      n_to_gen = std::atoll(argv[i+1]);
      i++;
    }
  }

  decaygen::driver driver;
  const decaygen::rotation rot = [](){
      decaygen::rotation R;

      // Rotation matrix using the 0,0,0 position for MicroBooNE (beam to det input) // From NuMI flux note
      R.rotate_axes({  0.92103853804025681562,    0.022713504803924120662,  0.38880857519374290021  }, 
                    {  4.6254001262154668408e-05, 0.99829162468141474651,  -0.058427989452906302359 }, 
                    { -0.38947144863934973769,    0.053832413938664107345,  0.91946400794392302291  }); // Also inverts so now to det to beam
      R.invert(); // go back to beam to det
      return R;
      }();

  //decaygen::vector3 vout = rot * decaygen::vector3{0.,0.,1.};
  //std::cout << vout.x() << " , "<<vout.y()<<" , "<<vout.z()<<std::endl;

  decaygen::geometry geo({0,1e4,1.e4},{1e3,1e3,1e3});
  driver.set_geometry(geo);
  decaygen::config conf;
  conf.fix_system_of_units(decaygen::config::GeV_cm_ns);
  conf.set_force_decays_in_detector(true);
  driver.set_config(conf);

  driver.add_particle_definition(
      decaygen::particle{1,100.,1e3}
      .add_decay({1.,{{2,false},{3,false}}})
      .finalise_decay_table()
      );
  
  driver.add_particle_definition(
      decaygen::particle{2,80.,.4e3}
      .add_decay({1.,{{11,true},{-11,true}}})
      .finalise_decay_table()
      );
  
  driver.add_particle_definition(
      decaygen::particle{3,15.,.2e3}
      .add_decay({1.,{{13,true},{-13,true},{11,true}}})
      .finalise_decay_table()
      );

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
        decaygen::decaying_particle_info::non_final
        },
        [&rng, &gen]()->double{return rng(gen);});

    if(res) {
      std::cout << res.build_hepevt_output().build_text() << std::endl;
    }
  }

  return 0;
}
