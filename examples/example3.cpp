#include "GELATO/core/driver.hpp"
#include "GELATO/core/decaying_particle_info.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <random>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>



//#include "../../development_misc/arXiv_1905_00284.hpp"
#include "GELATO/physics/heavy_neutral_leptons.hpp"
/*
 * Generates a decay table for a Heavy Neutral Lepton model,
 * based on arXiv:1905.00284, including 2-body and 3-body reweighters
 
  options:
    -n: number of initial mesons to decay (set to 0 for infinite loop, until -N condition is met)
    -N: numer of HNL decays in detector to output
    
    -m: mass of HNL
    -E: U_e4 model paremeter (default 0)
    -M: U_mu4 model paremeter (default 0)
    -T: U_tau4 model paremeter (default 0)
    -j: Majorana HNL (default false/Dirac)

    -r: output ROOT file
    -H: output HEPEVT text file
    
    -g: beamline/detector geometry (bnb or numi or sbnd) 
    
    -i: input flux file
    -p: POT in the flux file
    -F: HNL production modes to use (k+-, k0, pi, mu)
    
    -D: HNL decay modes to use (e_pi, mu_pi, or all)
    -f: force decays inside detector

    -w: unweighted mode (with estimate of max weight)
    -W: unweighted mode (with number of burn samples to estimate max weight)
    
    -d: debug mode (print decay table & exit)

  */


int main(int argc, char** argv) {
  size_t n_to_gen = 10;
  long long max_n_to_output = -1;
  double mass = -1, Ue4 = 0, Um4 = 0, Ut4 = 0;
  bool force = false, majorana = false, debug = false;
  std::string fluxfn = "";
  double pot_per_flux_file = 100e3;
  std::string rootfn = "";
  bool has_hepout = false;
  std::string hepfn = "";
  std::string decay_mode = "all";
  std::string flux_mode = "all";
  std::string geo_type = "bnb";
  const double most_negative_double = std::numeric_limits<double>::lowest();
  double max_weight = most_negative_double;
  size_t unweighted_burn_size = 0;
  for(int i = 0; i < argc; ++i) {
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
      geo_type = argv[i+1];
      i++;
      if(geo_type != "bnb" && geo_type != "numi" && geo_type != "sbnd") {
        std::cerr << " Geometry type "<<geo_type<<" must be 'bnb' or 'numi' or 'sbnd' (default bnb)"<<std::endl;
        return -1;
      }
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-w") { // unweighted mode; estimate of max weight
      max_weight = std::log(std::atof(argv[i+1]));
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-W") { // unweighted mode; number of burn samples to run
      unweighted_burn_size = std::atoll(argv[i+1]);
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
  if(mass < 0 || mass > .5) {
    std::cerr << " Must set a HNL mass (<0.5 GeV)" << std::endl;
    return 1;
  }


  GELATO::core::driver driver;

  GELATO::core::config conf;
  conf.fix_system_of_units(GELATO::core::config::system_of_units::GeV_cm_ns);
  conf.set_force_decays_in_detector(force);
  driver.set_config(conf);
  
  if(geo_type == "numi") {
  /*
   *                NUMI 
   */
  std::unique_ptr<GELATO::core::geometry> geo = std::make_unique<GELATO::core::box_geometry>(GELATO::core::vector3{31387.58422, 3316.402543, 60100.2414},GELATO::core::vector3{1.25e2,1.25e2,5e2},
    GELATO::core::rotation{}.rotate_axes({  0.92103853804025681562,    0.022713504803924120662,  0.38880857519374290021  }, 
                                        {  4.6254001262154668408e-05, 0.99829162468141474651,  -0.058427989452906302359 }, 
                                        { -0.38947144863934973769,    0.053832413938664107345,  0.91946400794392302291  }));
  driver.set_geometry(std::move(geo));

  } else if(geo_type=="bnb") {
  /*
   *                BNB
   */
  std::unique_ptr<GELATO::core::geometry> geo = std::make_unique<GELATO::core::box_geometry>(GELATO::core::vector3{0., 0., 475e2}, GELATO::core::vector3{1.25e2, 1.25e2, 5e2});
  driver.set_geometry(std::move(geo));
  }
  else if(geo_type == "sbnd") {
  std::unique_ptr<GELATO::core::geometry> geo = std::make_unique<GELATO::core::box_geometry>(GELATO::core::vector3{-74., 0., 110e2}, GELATO::core::vector3{2e2, 2e2, 2.5e2}); // sbnd is not directly on-axis, but offset by ~0.75m in x
  driver.set_geometry(std::move(geo));
  }
  
  //const double scalar_mass = mass; // GeV
  //const double scalar_theta = theta;
  //const std::string metadata = "";
  const std::string metadata{
    [=]() {
      std::ostringstream s;
      s
        << (majorana ? "type=Majorana" : "type=Dirac")
        <<" U_e4=" << Ue4
        <<" U_m4=" << Um4
        <<" U_t4=" << Ut4
        ;
      return s.str();}()
  };
  
  namespace hnl = GELATO::physics::heavy_neutral_leptons;
  //namespace hnl = GELATO::physics::arXiv_1905_00284;

  if(false) {

    // I was using this loop to do a branching ratio scan across many masses

    if(max_n_to_output < 1) max_n_to_output = 10;
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
        std::cout << "particle "<<p.pdg()<< " lifetime: ";
        if(p.lifetime() < 0) std::cout << "(stable)"<<std::endl;
        else std::cout << p.lifetime() << std::endl;
        size_t idgt = 0;
        for(auto d : p.get_decay_table()) {
          std::cout << "   "<<(++idgt==p.get_decay_table().size()?"└":"├")<<"─ decay BR="<<d.branching_ratio<< "  ───►";
          for(auto m : d.daughters) {
            std::cout <<" " << m.first << (m.second ? " (Final)":" (Decaying)");
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

    std::vector<GELATO::core::particle_info> input_flux;
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
      if(!fluxf.is_open()) {
        // probably file doesn't exist;
        std::cerr << "File "<<fluxfn<<" couldn't be opened."<<std::endl;
        return -1;
      }
      int nread = 0;
      while(!fluxf.eof() && fluxf.good()) {
        int pdg;
        double wt, vx,vy,vz,vt, px,py,pz,e;
        fluxf>>pdg;
        if(fluxf.eof() || !fluxf.good()) break;
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
        input_flux.push_back({pdg,{vx,vy,vz,vt-1},{vx,vy,vz,vt},{px,py,pz,e},std::log(wt)});
      }
      std::cout << "Number of flux entries read: "<<nread<<", stored: "<<input_flux.size()<<std::endl;
    }
    if(input_flux.empty()) {
      std::cerr << "no flux!"<< std::endl;
      return -1;
    }

    std::uniform_real_distribution<double> rng;
    std::default_random_engine gen;
    size_t i = 0;
    long long n_output = 0;
    auto fluxiter = input_flux.begin();
    int nloops = 0;
    double sum_weight = 0.;
    double pot_burnt = 0.;

    if(max_weight > most_negative_double || unweighted_burn_size > 0) {
      if(max_weight <= most_negative_double) {
        size_t nburnt = 0;
        while(nburnt < unweighted_burn_size) {
          auto const& rrr = driver.generate_decays(*fluxiter, [&rng, &gen]()->double{return rng(gen);});
          if(rrr) {
            auto const& res = rrr.build_hepevt_output();
            //std::cerr << "Burning "<<nburnt <<" with w="<<res.total_weight<<std::endl;
            if(res.total_log_weight > max_weight) {
              max_weight = res.total_log_weight+0.05; // std::log(1.1) ~= 0.05
              //std::cout << "Burnt "<<nburnt << " setting max_weight to "<<max_weight<<std::endl;
            }
            nburnt++;
          }
          fluxiter++;
          if(fluxiter==input_flux.end()) {
            fluxiter=input_flux.begin();
          }
        }
        std::cout << "Burnt "<<nburnt<<" events and found log_max_weight = "<<max_weight<<std::endl;
        pot_burnt = pot_per_flux_file * std::distance(input_flux.begin(),fluxiter)/input_flux.size();
      }
    }

    TFile *rootf = rootfn.empty() ? 0 : new TFile(rootfn.c_str(), "recreate");
    TTree *t = 0;
    TLorentzVector *p0pos = 0;
    TLorentzVector *p0mom = 0;
    TLorentzVector *p1 = 0;
    TLorentzVector *p2 = 0;
    TLorentzVector *p3 = 0;
    TLorentzVector *phnl = 0;
    TLorentzVector *decpos = 0;
    int dtype;
    double weight;

    if(rootf) {
      rootf->cd();
      t = new TTree("t","");
      p0pos = new TLorentzVector;
      p0mom = new TLorentzVector;
      p1 = new TLorentzVector;
      p2 = new TLorentzVector;
      p3 = new TLorentzVector;
      phnl = new TLorentzVector;
      decpos = new TLorentzVector;
      t->Branch("mesmom",&p0mom);
      t->Branch("mespos",&p0pos);
      t->Branch("p1",&p1); // momentum of first daughter
      t->Branch("p2",&p2); // momentum of second daughter
      t->Branch("p3",&p3); // momentum of third daughter
      t->Branch("phnl",&phnl);
      t->Branch("decpos",&decpos);
      t->Branch("type",&dtype);
      t->Branch("weight",&weight);
    }

    std::ofstream hepout;
    if(has_hepout) {
      hepout.open(hepfn);
    }
    while(n_to_gen==0 || i++ < n_to_gen) {
      auto const& res = driver.generate_decays(*fluxiter, [&rng, &gen]()->double{return rng(gen);});
      if(res && (max_n_to_output  < 1 || n_output < max_n_to_output)) {
        auto hepevt = res.build_hepevt_output();
        if(max_weight > most_negative_double) {
          const double weight = hepevt.total_log_weight;
          if(weight > max_weight) {
            std::cerr << "Error! Weight found "<<weight<<" > max_weight "<<max_weight<<"; output may be unrepresentative of truth."<<std::endl;
            break;
          }
          const double u = rng(gen);
          if(weight < std::log(u) + max_weight) continue;
          hepevt.total_log_weight = 0.;
        }
        sum_weight += std::exp(hepevt.total_log_weight);
        if(has_hepout) {
          hepout << hepevt.build_text();
        }
        if(rootf) {
          // reset p3 in case of 2-body decay
          p3->SetXYZT(0.,0.,0.,0.);
          size_t imoms = 0;
          auto moms = std::vector<TLorentzVector*>{p1,p2,p3};
          std::vector<int> pdgs;
          weight = std::exp(hepevt.total_log_weight);
          int number_of_particle_in_list = 1;
          int number_of_HNL = 0;
          for(auto const& p : hepevt.particle_info) {
            if(p.status == 0) {
              p0pos->SetXYZT(p.position.x(),p.position.y(),p.position.z(),p.position.t());
              p0mom->SetXYZT(p.momentum.x(),p.momentum.y(),p.momentum.z(),p.momentum.t());
            }
            else if(p.status == 2 && std::abs(p.pdg)>=89 && std::abs(p.pdg) <= 91) {
              number_of_HNL = number_of_particle_in_list;
              phnl->SetXYZT(p.momentum.x(),p.momentum.y(),p.momentum.z(),p.momentum.t());
            }
            else if(p.status == 1 && p.first_mother == number_of_HNL) {
              decpos->SetXYZT(p.position.x(),p.position.y(),p.position.z(),p.position.t());
              pdgs.push_back(std::abs(p.pdg));
              moms[imoms++]->SetXYZT(p.momentum.x(),p.momentum.y(),p.momentum.z(),p.momentum.t());
            }
            number_of_particle_in_list++;
          }

          // dtypes:
          dtype=0;
          std::sort(pdgs.begin(),pdgs.end());
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
        }
        if(++n_output >= max_n_to_output) break;
      }
      fluxiter++;
      if(fluxiter==input_flux.end()) {
        fluxiter=input_flux.begin();
        nloops++;
      }
    }
    const double tot_pot = (max_weight > most_negative_double ? 1./std::exp(max_weight) : 1.)
      * (pot_per_flux_file * (nloops + std::distance(input_flux.begin(),fluxiter) * 1./input_flux.size()) - pot_burnt);
    if(rootf) {
      t->SetTitle(Form("pot: %f",tot_pot));
      t->Write();
      rootf->Close();
      delete rootf;
      rootf=0;
    }
    if(has_hepout) {
      hepout << "#END total_weight="<<sum_weight<<" pot_per_flux_file=" << pot_per_flux_file
      << " total_pot="<<tot_pot<<" "<<metadata<<std::endl;
    }
  }

  return 0;
}
