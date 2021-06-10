#include "dkgen/core/driver.hpp"
#include "dkgen/core/decaying_particle_info.hpp"

#include <gsl/gsl_integration.h>

#include "dkgen/physics/heavy_neutral_leptons.hpp"
//#include "../../development_misc/arXiv_1905_00284.hpp"

#include <iostream>
#include <sstream>
#include <cmath>
#include <random>

int main(int argc, char** argv) {
  size_t n_to_gen = 10;
  long long max_n_to_output = -1;
  double mass = 0.140, Ue4 = 0, Um4 = 0, Ut4 = 0;
  bool force = false, majorana = false, debug = false;
  for(int i = 0; i < argc /* as -n option will have argument */; ++i) {
    if(i+1 < argc && std::string(argv[i]) == "-n") { // number initial particles to gen
      n_to_gen = std::atoll(argv[i+1]);
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-N") { // max number of ouput events
      max_n_to_output = std::atoll(argv[i+1]);
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-m") { // mass
      mass = std::atof(argv[i+1]);
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-E") { // Ue4 angle
      Ue4 = std::atof(argv[i+1]);
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-M") { // Umu4 angle
      Um4 = std::atof(argv[i+1]);
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-T") { // Utau4 angle
      Ut4 = std::atof(argv[i+1]);
      i++;
      continue;
    }
    if(std::string(argv[i]) == "-f") {
      force = true;
      continue;
    }
    if(std::string(argv[i]) == "-j") { // 'j' for Majorana (m/M already taken)
      majorana = true;
      continue;
    }
    if(std::string(argv[i]) == "-d") { 
      debug = true;
      continue;
    }
  }
  if(Ue4 == 0. && Um4 == 0. && Ut4 == 0.) {
    std::cerr << " Must set a mixing angle -E/-M/-T" << std::endl;
    return 1;
  }

  dkgen::core::driver driver;

  dkgen::core::geometry geo({0,0,1.e4},{1e3,1e3,1e3});
  driver.set_geometry(geo);
  dkgen::core::config conf;
  conf.fix_system_of_units(dkgen::core::config::system_of_units::GeV_cm_ns);
  conf.set_force_decays_in_detector(force);
  driver.set_config(conf);

  //const double scalar_mass = mass; // GeV
  //const double scalar_theta = theta;
  //const std::string metadata = "";
  const std::string metadata{
    (std::ostringstream()
     << (majorana ? "type=Majorana" : "type=Dirac")
     <<" U_e4=" << Ue4
     <<" U_m4=" << Um4
     <<" U_t4=" << Ut4
     ).str()
  };
  if(max_n_to_output < 1) max_n_to_output = 10;
  namespace hnl = dkgen::physics::heavy_neutral_leptons;
  //namespace hnl = dkgen::physics::arXiv_1905_00284;
  for(int mm = 1; mm <= max_n_to_output; ++mm) {
    const double delta = 1./max_n_to_output;
    double mass = delta * mm * conf.physical_params().find_particle("kaon_pm").mass;
    //std::cout << mass << std::endl;
    const hnl::model_parameters params{mass, Ue4, Um4, Ut4, false};
    auto const& particles = hnl::create_particle_content(params,conf,
        hnl::all_decay_modes,
        {hnl::production_modes::k_e2}, false
        );
    for(auto& p : particles) {
      if(p.pdg() == 91) {
        std::cout << mass;
        for (auto & d : p.get_decay_table()) {
          std::cout << " (" << d.branching_ratio;
          for(auto& m : d.daughters) {
            std::cout <<" " << m.first;
          }
          std::cout << ")";
        }
        std::cout << std::endl;
      }
    }
  }

  if (false) {

    //const std::string metadata = [](){ auto s = std::istringstream(); s << "model_theta=" << scalar_theta; return s.str(); }();
    const hnl::model_parameters params{mass, Ue4, Um4, Ut4, majorana};
    auto const& particles = hnl::create_particle_content(params,conf /*,
        hnl::all_decay_modes,
        {hnl::production_modes::mu_e}*/
        );

         if(debug) {
           for(auto p : particles) {
             std::cout << "particle "<<p.pdg()<<std::endl;
             for(auto d : p.get_decay_table()) {
               std::cout << "  decay BR="<<d.branching_ratio;
               for(auto m : d.daughters) {
                 std::cout <<" " << m.first;
               }
               std::cout <<std::endl;
             }
           }
           return 0;
         }

    driver.set_particle_content(particles);

    auto& kaon = conf.physical_params().find_particle("kaon_pm");
    const int kpdg = -kaon.pdgcode;
    const double kmass = kaon.mass;
    const double kmom = 0.5;

    std::uniform_real_distribution<double> rng;
    std::default_random_engine gen;
    //std::vector<double> moms(n_to_gen);
    size_t i = 0;
    long long n_output = 0;
    while(i++ < n_to_gen) {
      auto const& res = driver.generate_decays(
          {
          kpdg,
          {0.,0.,0.,0.}, // production position
          {0.,0.,0.,0.}, // decay position
          {0.,0.,kmom,std::sqrt(kmom*kmom + kmass*kmass)}, // momentum
          },
          [&rng, &gen]()->double{return rng(gen);});
      if(res && (max_n_to_output  < 1 || n_output++ < max_n_to_output)) {
        std::cout << res.build_hepevt_output().build_text(metadata.c_str()) << std::endl;
      }
    }
  }

  return 0;
}
