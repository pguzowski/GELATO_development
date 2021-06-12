#include "dkgen/core/driver.hpp"
#include "dkgen/core/decaying_particle_info.hpp"

#include <gsl/gsl_integration.h>

#include "dkgen/physics/heavy_neutral_leptons.hpp"
//#include "../../development_misc/arXiv_1905_00284.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <random>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

int main(int argc, char** argv) {
  size_t n_to_gen = 10;
  long long max_n_to_output = -1;
  double mass = 0.140, Ue4 = 0, Um4 = 0, Ut4 = 0;
  bool force = false, majorana = false, debug = false;
  std::string fluxfn = "";
  double pot_per_flux_file = 100e3;
  std::string rootfn = "out.root";
  bool has_hepout = false;
  std::string hepfn = "";
  std::string decay_mode = "all";
  std::string flux_mode = "all";
  bool is_numi = false;
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
    if(i+1 < argc && std::string(argv[i]) == "-i") { // input flux
      fluxfn = argv[i+1];
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-r") { // output root file
      rootfn = argv[i+1];
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-H") { // output hepevt file
      hepfn = argv[i+1];
      has_hepout = true;
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-p") { // POT in the flux file
      pot_per_flux_file = std::atof(argv[i+1]);
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-D") { // decay modes to simulate
      decay_mode = argv[i+1];
      i++;
      if(decay_mode != "all" && decay_mode != "e_pi" && decay_mode != "mu_pi") {
        std::cerr << "Unrecognised decay mode (-D) "<<decay_mode<<" (must be e_pi or mu_pi or all)"<<std::endl;
        return -1;
      }
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-F") { // production modes to simulate
      flux_mode = argv[i+1];
      i++;
      if(flux_mode != "all" && flux_mode != "mu" && flux_mode != "pi" && flux_mode != "k+-" && flux_mode != "k0") {
        std::cerr << "Unrecognised flux mode (-F) '"<<flux_mode<<"' (must be 'k+-' or 'k0' or 'mu' or 'pi' or 'all')"<<std::endl;
        return -1;
      }
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-g") { // production modes to simulate
      std::string arg = argv[i+1];
      i++;
      if(arg == "numi") {
        is_numi = true;
      }
      else if(arg == "bnb") {
        is_numi = false;
      }
      else {
        std::cerr << " Geometry type "<<arg<<" must be 'bnb' or 'numi' (default bnb)"<<std::endl;
        return -1;
      }
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

  dkgen::core::config conf;
  conf.fix_system_of_units(dkgen::core::config::system_of_units::GeV_cm_ns);
  conf.set_force_decays_in_detector(force);
  driver.set_config(conf);
  
  if(is_numi) {
  /*
   *                NUMI 
   */
  dkgen::core::geometry geo({31387.58422, 3316.402543, 60100.2414},{1.25e2,1.25e2,5e2},
    dkgen::core::rotation{}.rotate_axes({  0.92103853804025681562,    0.022713504803924120662,  0.38880857519374290021  }, 
                                        {  4.6254001262154668408e-05, 0.99829162468141474651,  -0.058427989452906302359 }, 
                                        { -0.38947144863934973769,    0.053832413938664107345,  0.91946400794392302291  }));
  driver.set_geometry(geo);

  } else {
  /*
   *                BNB
   */
  dkgen::core::geometry geo({0., 0., 475e2}, {1.25e2, 1.25e2, 5e2});
  driver.set_geometry(geo);
  }
  
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
  namespace hnl = dkgen::physics::heavy_neutral_leptons;
  if(false) {
    if(max_n_to_output < 1) max_n_to_output = 10;
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
  }

  if (true) {

    const hnl::model_parameters params{mass, Ue4, Um4, Ut4, majorana};
    const std::vector<hnl::decay_modes> decay_modes = (decay_mode=="mu_pi"
        ? std::vector<hnl::decay_modes>{hnl::decay_modes::mu_pi}
        : (decay_mode == "e_pi"
          ? std::vector<hnl::decay_modes>{hnl::decay_modes::e_pi}
          : hnl::all_decay_modes
          )
        );
    const std::vector<hnl::production_modes> prod_modes =
      (flux_mode=="k+-"
       ? hnl::all_kaon_pm_production_modes
       : (flux_mode == "k0"
         ? hnl::all_kaon_0L_production_modes
         : (flux_mode == "pi"
           ? hnl::all_pion_production_modes
           : (flux_mode == "mu"
             ? hnl::all_muon_production_modes
             : hnl::all_production_modes
             )
           )
         )
      );
    auto const& particles = hnl::create_particle_content(params,conf, decay_modes, prod_modes, !debug );

    if(debug) {
      for(auto p : particles) {
        std::cout << "particle "<<p.pdg()<<std::endl;
        for(auto d : p.get_decay_table()) {
          std::cout << "  decay BR="<<d.branching_ratio;
          for(auto m : d.daughters) {
            std::cout <<" " << m.first << (m.second ? " F":" D");
          }
          std::cout <<std::endl;
        }
      }
      return 0;
    }

    driver.set_particle_content(particles);

    const double kaon_pm_br = hnl::kaon_pm_hnl_branching_ratio(params, conf, prod_modes);
    const double kaon_0L_br = hnl::kaon_0L_hnl_branching_ratio(params, conf, prod_modes);
    const double pion_br = hnl::pion_hnl_branching_ratio(params, conf, prod_modes);
    const double muon_br = hnl::muon_hnl_branching_ratio(params, conf, prod_modes);

    std::vector<dkgen::core::particle_info> input_flux;
    if(fluxfn.empty()) {
      auto& kaon = conf.physical_params().find_particle("kaon_pm");
      const int kpdg = -kaon.pdgcode;
      const double kmass = kaon.mass;
      const double kmom = 0.5;
      input_flux.push_back({
          kpdg,
          {0.,0.,0.,0.}, // production position
          {0.,0.,0.,0.}, // decay position
          {0.,0.,kmom,std::sqrt(kmom*kmom + kmass*kmass)}, // momentum
          });
    }
    else {
      std::ifstream fluxf{fluxfn};
      int nread = 0;
      while(!fluxf.eof()) {
        int pdg;
        double wt, vx,vy,vz,vt, px,py,pz,e;
        fluxf>>pdg;
        if(fluxf.eof()) break;
        fluxf>>wt>>vx>>vy>>vz>>vt>>px>>py>>pz>>e;
        nread++;
        if(std::abs(pdg) == 321) {
          if(kaon_pm_br <= 0) continue;
          wt *= kaon_pm_br;
        }
        else if(std::abs(pdg) == 130) {
          if(kaon_0L_br <= 0) continue;
          wt *= kaon_0L_br;
        }
        else if(std::abs(pdg) == 211) {
          if(pion_br <= 0) continue;
          wt *= pion_br;
        }
        else if(std::abs(pdg) == 13) {
          if(muon_br <= 0) continue;
          wt *= muon_br;
        }
        input_flux.push_back({pdg,{vx,vy,vz,vt-1},{vx,vy,vz,vt},{px,py,pz,e},wt});
      }
      std::cout << "Number of flux entries read: "<<nread<<", stored: "<<input_flux.size()<<std::endl;
    }
    if(input_flux.empty()) {
      std::cerr << "no flux!"<< std::endl;
      return -1;
    }

    std::uniform_real_distribution<double> rng;
    std::default_random_engine gen;
    //std::vector<double> moms(n_to_gen);
    size_t i = 0;
    long long n_output = 0;
    auto fluxiter = input_flux.begin();
    int nloops = 0;

    TFile *f = new TFile(rootfn.c_str(),"recreate");
    TTree *t = new TTree("t","");
    TLorentzVector *p1 = new TLorentzVector;
    TLorentzVector *p2 = new TLorentzVector;
    TLorentzVector *p3 = new TLorentzVector;
    int dtype;
    double weight;
    t->Branch("p1",&p1);
    t->Branch("p2",&p2);
    t->Branch("p3",&p3);
    t->Branch("type",&dtype);
    t->Branch("weight",&weight);

    std::ofstream hepout;
    if(has_hepout) {
      hepout.open(hepfn);
    }
    while(i++ < n_to_gen) {
      auto const& res = driver.generate_decays(*fluxiter, [&rng, &gen]()->double{return rng(gen);});
      if(res && (max_n_to_output  < 1 || n_output++ < max_n_to_output)) {
        if(has_hepout) {
          hepout << res.build_hepevt_output().build_text();
        }
        auto const& hepevt = res.build_hepevt_output();
        p3->SetXYZT(0.,0.,0.,0.);
        size_t imoms = 0;
        auto moms = std::vector<TLorentzVector*>{p1,p2,p3};
        std::vector<int> pdgs;
        weight = hepevt.total_weight;
        int number_of_particle_in_list = 1;
        int number_of_HNL = 0;
        for(auto const& p : hepevt.particle_info) {
          if(p.status == 2 && std::abs(p.pdg)>=89 && std::abs(p.pdg) <= 91) {
            number_of_HNL = number_of_particle_in_list;
          }
          else if(p.status == 1 && p.first_mother == number_of_HNL) {
            pdgs.push_back(std::abs(p.pdg));
            moms[imoms++]->SetXYZT(p.momentum.x(),p.momentum.y(),p.momentum.z(),p.momentum.t());
          }
          number_of_particle_in_list++;
        }
        std::sort(pdgs.begin(),pdgs.end());
        dtype=0;
        //std::cerr << pdgs.size() << std::endl;
        if(pdgs.size() == 2) {
          if(pdgs[0]==13 && pdgs[1]==211) dtype = 1; // mu pi
          else if(pdgs[0]==11 && pdgs[1]==211) dtype = 2; // e pi
          else if(pdgs[0]==12 && pdgs[1]==111) dtype = 3; // nu pi0
          else std::cerr << "unknown decay "<<pdgs[0]<<" "<<pdgs[1]<<std::endl;
        }
        else if(pdgs.size() == 3) {
          if(pdgs[0]==11 && pdgs[1]==12 && pdgs[2] ==13) dtype=4; // nu e mu
          else if(pdgs[0]==11 && pdgs[1]==11 && pdgs[2] ==12) dtype=5; // nu e e
          else if(pdgs[0]==12 && pdgs[1]==13 && pdgs[2] ==13) dtype=6; // nu mu mu
          else std::cerr << "unknown decay "<<pdgs[0]<<" "<<pdgs[1]<<" " <<pdgs[2]<<std::endl;
        }
        t->Fill();
        if(n_output >= max_n_to_output) break;
      }
      fluxiter++;
      if(fluxiter==input_flux.end()) {
        fluxiter=input_flux.begin();
        nloops++;
      }
    }
    t->Write();
    f->Close();
    const double tot_pot = nloops * pot_per_flux_file + pot_per_flux_file * std::distance(input_flux.begin(),fluxiter)/input_flux.size();
    if(has_hepout) {
      hepout << "# weight=0 pot_per_flux_file=" << pot_per_flux_file
      << " total_pot="<<tot_pot<<" "<<metadata<<std::endl;
    }
  }

  return 0;
}
