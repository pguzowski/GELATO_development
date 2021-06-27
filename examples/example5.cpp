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
#include "GELATO/physics/leptophilic_gauge_boson.hpp"
/*
 * Generates a decay table for a Heavy Neutral Lepton model,
 * based on arXiv:1905.00284, including 2-body and 3-body reweighters
 
  options:
    -n: number of initial mesons to decay (set to 0 for infinite loop, until -N condition is met)
    -N: numer of boson decays in detector to output
    
    -m: mass of boson
    -G: gauge coupling
    -L: gauge symmetry (e-mu, e-tau, mu-tau)

    -r: output ROOT file
    -H: output HEPEVT text file
    
    -g: beamline/detector geometry (bnb or numi or sbnd) 
    
    -i: input flux file in the format file_name:POT
    -F: boson production modes to use (k+-, pi+-, pi0, eta)
    
    -D: boson decay modes to use (e_e, mu_mu, tau_tau, or all)
    -f: force decays inside detector

    -w: unweighted mode (with estimate of max weight)
    -W: unweighted mode (with number of burn samples to estimate max weight)
    
    -d: debug mode (print decay table & exit)

  */


int main(int argc, char** argv) {
  size_t n_to_gen = 10;
  long long max_n_to_output = -1;
  double mass = -1, coupling = -1;
  std::string gaugestr = "";
  namespace lgb = GELATO::physics::leptophilic_gauge_boson;
  lgb::leptophilic_gauges gauge = lgb::leptophilic_gauges::Le_minus_Ltau;
  bool force = false, majorana = false, debug = false;
  std::string fluxfn = "";
  double pot_per_flux_file = 100e3;
  std::string rootfn = "";
  bool has_hepout = false;
  std::string hepfn = "";
  std::string decay_mode = "all";
  std::string flux_mode = "all";
  std::string geo_type = "bnb";
  double max_weight = -1.;
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
    if(i+1 < argc && std::string(argv[i]) == "-G") { //  gauge coupling
      coupling = std::atof(argv[i+1]);
      i++;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-L") { // which gauge
      gaugestr = argv[i+1];
      i++;
      if(gaugestr != "e-mu" && gaugestr != "e-tau" && gaugestr != "mu-tau") {
        std::cerr << "Unrecognised gauge (-L) "<<gaugestr<<" (must be e-mu or e-tau or mu-tau)"<<std::endl;
        return -1;
      }
      if(gaugestr == "e-mu")    gauge = lgb::leptophilic_gauges::Le_minus_Lmu;
      if(gaugestr == "e-tau")   gauge = lgb::leptophilic_gauges::Le_minus_Ltau;
      if(gaugestr == "mu-tau")  gauge = lgb::leptophilic_gauges::Lmu_minus_Ltau;
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-i") { // input flux
      std::string fluxfn_plus_pot = argv[i+1];
      i++;
      const size_t colonpos = fluxfn_plus_pot.find_last_of(":");
      if(colonpos == std::string::npos) {
        std::cerr << "Flux filename must be in the format 'file_name:POT'" << std::endl;
        return -1;
      }
      fluxfn = fluxfn_plus_pot.substr(0,colonpos);
      pot_per_flux_file = std::atof(fluxfn_plus_pot.substr(colonpos+1).c_str());
      if(pot_per_flux_file <= 0) {
        std::cerr << "Couldn't parse POT="<<fluxfn_plus_pot.substr(colonpos+1)<<std::endl;
        return -1;
      }
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
    if(i+1 < argc && std::string(argv[i]) == "-D") { // decay modes to simulate
      decay_mode = argv[i+1];
      i++;
      if(decay_mode != "all" && decay_mode != "e_e" && decay_mode != "mu_mu" && decay_mode != "tau_tau") {
        std::cerr << "Unrecognised decay mode (-D) "<<decay_mode<<" (must be e_e or mu_mu or tau_tau or all)"<<std::endl;
        return -1;
      }
      continue;
    }
    if(i+1 < argc && std::string(argv[i]) == "-F") { // production modes to simulate
      flux_mode = argv[i+1];
      i++;
      if(flux_mode != "all" && flux_mode != "eta" && flux_mode != "pi+-" && flux_mode != "k+-" && flux_mode != "pi0") {
        std::cerr << "Unrecognised flux mode (-F) '"<<flux_mode
          <<"' (must be 'k+-' or 'pi+-' or 'eta' or 'pi0' or 'all')"<<std::endl;
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
      max_weight = std::atof(argv[i+1]);
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
  if(coupling < 0.) {
    std::cerr << " Must set a gauge coupling -G" << std::endl;
    return 1;
  }
  if(gaugestr == "") {
    std::cerr << " Must set a gauge symmetry -L (e-mu, e-tau, mu-tau)" << std::endl;
    return 1;
  }
  if(mass < 0 || mass > .5) {
    std::cerr << " Must set a boson mass (<0.5 GeV)" << std::endl;
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
  GELATO::core::geometry geo({31387.58422, 3316.402543, 60100.2414},{1.25e2,1.25e2,5e2},
    GELATO::core::rotation{}.rotate_axes({  0.92103853804025681562,    0.022713504803924120662,  0.38880857519374290021  }, 
                                        {  4.6254001262154668408e-05, 0.99829162468141474651,  -0.058427989452906302359 }, 
                                        { -0.38947144863934973769,    0.053832413938664107345,  0.91946400794392302291  }));
  driver.set_geometry(geo);

  } else if(geo_type=="bnb") {
  /*
   *                BNB
   */
  GELATO::core::geometry geo({0., 0., 475e2}, {1.25e2, 1.25e2, 5e2});
  driver.set_geometry(geo);
  }
  else if(geo_type == "sbnd") {
  GELATO::core::geometry geo({-74., 0., 110e2}, {2e2, 2e2, 2.5e2}); // sbnd is not directly on-axis, but offset by ~0.75m in x
  driver.set_geometry(geo);
  }
  
  //const double scalar_mass = mass; // GeV
  //const double scalar_theta = theta;
  //const std::string metadata = "";
  const std::string metadata{
    [=]() {
      std::ostringstream s;
      s
        <<" g=" << coupling
        <<" L=" << gaugestr
        ;
      return s.str();}()
  };
  

  if(false) {

    // I was using this loop to do a branching ratio scan across many masses

    if(max_n_to_output < 1) max_n_to_output = 10;
    for(int mm = 1; mm <= max_n_to_output; ++mm) {
      const double delta = 1./max_n_to_output;
      double mass = delta * mm * conf.physical_params().find_particle("kaon_pm").mass;
      //std::cout << mass << std::endl;
      const lgb::model_parameters params{mass, coupling, gauge};
      auto const& particles = lgb::create_particle_content(params,conf,
          lgb::all_decay_modes,
          {lgb::production_modes::k_mu}, false
          );
      for(auto& p : particles) {
        if(p.pdg() == 32) {
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

    const lgb::model_parameters params{mass, coupling, gauge};
    const std::vector<lgb::decay_modes> decay_modes = (decay_mode=="mu_mu"
        ? std::vector<lgb::decay_modes>{lgb::decay_modes::mu_mu}
        : (decay_mode == "e_e"
          ? std::vector<lgb::decay_modes>{lgb::decay_modes::e_e}
          : (decay_mode == "tau_tau"
            ? std::vector<lgb::decay_modes>{lgb::decay_modes::tau_tau}
            : lgb::all_decay_modes
            )
          )
        );
    const std::vector<lgb::production_modes> prod_modes =
      (flux_mode=="k+-"
       ? lgb::all_kaon_pm_production_modes
       : (flux_mode == "pi+-"
         ? lgb::all_pion_pm_production_modes
         : (flux_mode == "pi0"
           ? std::vector<lgb::production_modes>{lgb::production_modes::pi0}
           : (flux_mode == "eta"
             ? std::vector<lgb::production_modes>{lgb::production_modes::eta}
             : lgb::all_production_modes
             )
           )
         )
      );
    auto const& particles = lgb::create_particle_content(params,conf, decay_modes, prod_modes, !debug );

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

    const double kaon_pm_br = lgb::kaon_pm_lgb_branching_ratio(params, conf, prod_modes);
    const double pion_pm_br = lgb::pion_lgb_branching_ratio(params, conf, prod_modes);
    const double pi0_br = lgb::pi0_lgb_branching_ratio(params, conf, prod_modes);
    const double eta_br = lgb::eta_lgb_branching_ratio(params, conf, prod_modes);

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
        else if(std::abs(pdg) == 211) {
          if(pion_pm_br <= 0) continue;
          wt *= pion_pm_br;
        }
        else if(std::abs(pdg) == 111) {
          if(pi0_br <= 0) continue;
          wt *= pi0_br;
        }
        else if(std::abs(pdg) == 221) {
          if(eta_br <= 0) continue;
          wt *= eta_br;
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
    size_t i = 0;
    long long n_output = 0;
    auto fluxiter = input_flux.begin();
    int nloops = 0;
    double sum_weight = 0.;
    double pot_burnt = 0.;

    if(max_weight > 0. || unweighted_burn_size > 0) {
      if(max_weight <= 0.) {
        size_t nburnt = 0;
        while(nburnt < unweighted_burn_size) {
          auto const& rrr = driver.generate_decays(*fluxiter, [&rng, &gen]()->double{return rng(gen);});
          if(rrr) {
            auto const& res = rrr.build_hepevt_output();
            //std::cerr << "Burning "<<nburnt <<" with w="<<res.total_weight<<std::endl;
            if(res.total_weight > max_weight) {
              max_weight = 1.1*res.total_weight;
              //std::cout << "Burnt "<<nburnt << " setting max_weight to "<<max_weight<<std::endl;
            }
            nburnt++;
          }
          fluxiter++;
          if(fluxiter==input_flux.end()) {
            fluxiter=input_flux.begin();
          }
        }
        std::cout << "Burnt "<<nburnt<<" events and found max_weight = "<<max_weight<<std::endl;
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
    TLorentzVector *pboson = 0;
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
      pboson = new TLorentzVector;
      decpos = new TLorentzVector;
      t->Branch("mesmom",&p0mom);
      t->Branch("mespos",&p0pos);
      t->Branch("p1",&p1); // momentum of first daughter
      t->Branch("p2",&p2); // momentum of second daughter
      t->Branch("p3",&p3); // momentum of third daughter
      t->Branch("pboson",&pboson);
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
        if(max_weight > 0) {
          const double weight = hepevt.total_weight;
          if(weight > max_weight) {
            std::cerr << "Error! Weight found "<<weight<<" > max_weight "<<max_weight<<"; output may be unrepresentative of truth."<<std::endl;
            break;
          }
          const double u = rng(gen);
          if(weight < u * max_weight) continue;
          hepevt.total_weight = 1.;
        }
        sum_weight += hepevt.total_weight;
        if(has_hepout) {
          hepout << hepevt.build_text();
        }
        if(rootf) {
          // reset p3 in case of 2-body decay
          p3->SetXYZT(0.,0.,0.,0.);
          size_t imoms = 0;
          auto moms = std::vector<TLorentzVector*>{p1,p2,p3};
          std::vector<int> pdgs;
          weight = hepevt.total_weight;
          int number_of_particle_in_list = 1;
          int number_of_boson = 0;
          for(auto const& p : hepevt.particle_info) {
            if(p.status == 0) {
              p0pos->SetXYZT(p.position.x(),p.position.y(),p.position.z(),p.position.t());
              p0mom->SetXYZT(p.momentum.x(),p.momentum.y(),p.momentum.z(),p.momentum.t());
            }
            else if(p.status == 2 && p.pdg == 32) {
              number_of_boson = number_of_particle_in_list;
              pboson->SetXYZT(p.momentum.x(),p.momentum.y(),p.momentum.z(),p.momentum.t());
            }
            else if(p.status == 1 && p.first_mother == number_of_boson) {
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
            if(pdgs[0]==11 && pdgs[1]==-11) dtype = 1; // mu pi
            else if(pdgs[0]==13 && pdgs[1]==-13) dtype = 2; // e pi
            else if(pdgs[0]==15 && pdgs[1]==-15) dtype = 3; // nu pi0
            else std::cerr << "unknown decay "<<pdgs[0]<<" "<<pdgs[1]<<std::endl;
          }
          else if(pdgs.size() == 3) {
            std::cerr << "unknown decay "<<pdgs[0]<<" "<<pdgs[1]<<" " <<pdgs[2]<<std::endl;
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
    const double tot_pot = (max_weight > 0. ? 1./max_weight : 1.)
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


#if 0


#include "GELATO/physics/leptophilic_gauge_boson.hpp"

#include "GELATO/core/driver.hpp"
#include "GELATO/core/decaying_particle_info.hpp"

#include <iostream>

int main(int , char** ) {
  namespace lgb = GELATO::physics::leptophilic_gauge_boson;

  GELATO::core::config conf;
  conf.fix_system_of_units(GELATO::core::config::system_of_units::GeV_cm_ns);
  //conf.set_force_decays_in_detector(force);
  
  for(int m = 1; m < 400; ++m) {
    const double boson_mass = 1e-3 * m;
    const double gauge_coupling = 1e-3;
    lgb::model_parameters params{boson_mass, gauge_coupling, lgb::leptophilic_gauges::Lmu_minus_Ltau };
    //auto const& particles = lgb::create_particle_content(params, conf, lgb::all_decay_modes, lgb::all_production_modes, false);
    //double br_k_mu = 0., br_k_e = 0., br_pi_mu = 0., br_pi_e = 0.;
    //for(auto const& p : particles) {
    //}
    const double br_k_mu  = lgb::kaon_pm_lgb_branching_ratio(params, conf, {lgb::production_modes::k_mu});
    const double br_k_e   = lgb::kaon_pm_lgb_branching_ratio(params, conf, {lgb::production_modes::k_e});
    const double br_pi_mu = lgb::pion_lgb_branching_ratio(params, conf, {lgb::production_modes::pi_mu});
    const double br_pi_e  = lgb::pion_lgb_branching_ratio(params, conf, {lgb::production_modes::pi_e});
    const double br_pi0  = lgb::pi0_lgb_branching_ratio(params, conf, {lgb::production_modes::pi0});
    const double br_eta  = lgb::eta_lgb_branching_ratio(params, conf, {lgb::production_modes::eta});

    std::cout << boson_mass << " " << br_k_mu <<" " <<br_k_e<<" "<<br_pi_mu<<" "<<br_pi_e<<" "<<br_pi0<<" "<<br_eta <<'\n';
  }
}
#endif
