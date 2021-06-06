#include "dkgen/core/config.hpp"

#include <stdexcept>
#include <limits>

#ifdef XUSING_ROOT
#include <TDatabasePDG.h>
#endif

dkgen::core::config& dkgen::core::config::fix_system_of_units(system_of_units sys) {
  switch(sys) {
    case system_of_units::GeV_cm_ns:
      pparams.speed_of_light = 29.9791932;
      pparams.hbar = 6.582119569e-16;
      pparams.higgs_vev = 246.22;
      pparams.CKM_Vtd = 8.1e-3;
      pparams.CKM_Vts = 39.4e-3;
#ifdef XUSING_ROOT
      {
        TDatabasePDG db{};
        auto kaon_pm = db.GetParticle("K+");
        auto kaon_0L = db.GetParticle("K_L0");
        auto pion_pm = db.GetParticle("pi+");
        auto pion_0 = db.GetParticle("pi0");
        auto elec = db.GetParticle("e-");
        auto muon = db.GetParticle("mu-");
        auto top = db.GetParticle("t");
        const double seconds_to_ns = 1e9; // TDatabasePDG lifetimes are in seconds
        // some lifetimes are not present in the database (set to 0)
        pparams.particles = {
          {"kaon_pm", { kaon_pm->PdgCode(), kaon_pm->Mass(), kaon_pm->Lifetime() * seconds_to_ns }},
          {"kaon_0L", { kaon_0L->PdgCode(), kaon_0L->Mass(), 
                          kaon_0L->Lifetime() > 0 ? kaon_0L->Lifetime() * seconds_to_ns : 51.16 }},
          {"pion_pm", { pion_pm->PdgCode(), pion_pm->Mass(), pion_pm->Lifetime() * seconds_to_ns }},
          {"pion_0",  { pion_0->PdgCode(), pion_0->Mass(),
                          pion_0->Lifetime() > 0 ? pion_0->Lifetime() * seconds_to_ns : 8.43e-8 }},
          {"elec",    { elec->PdgCode(),    elec->Mass(),
                        elec->Lifetime() > 0 ? elec->Lifetime() * seconds_to_ns : std::numeric_limits<double>::infinity() }}, 
          {"muon",    { muon->PdgCode(),    muon->Mass(),
                          muon->Lifetime() > 0 ? muon->Lifetime() * seconds_to_ns : 2197. }},
          {"top",     { top->PdgCode(),     top->Mass(),     top->Lifetime() * seconds_to_ns     }},
        };
      }
#else // not USING_ROOT
      pparams.particles = {
        // name,      pdg,   mass,          lifetime
        {"kaon_pm", { 321, 0.493677,        12.38 }},
        {"kaon_0L", { 130, 0.497648,        51.16 }},
        {"pion_pm", { 211, 0.13957061,      26.03 }},
        {"pion_0",  { 111, 0.1349770,       8.43e-8 }},
        {"elec",    {  11, 0.5109989461e-3, std::numeric_limits<double>::infinity() }}, 
        {"muon",    {  13, 0.1056583745,    2197. }},
        {"top",     {   6, 172.9,           0. }},
      };
#endif
      break;
  }
  return *this;
}

const dkgen::core::config::standard_particle&
dkgen::core::config::physical_parameters::find_particle(const std::string& name) const {
  auto p = particles.find(name);
  if(p == particles.end()) {
    throw std::runtime_error("Error! "+name+" not found in list of particles");
  }
  return p->second;
}

