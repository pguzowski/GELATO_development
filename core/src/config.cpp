#include "GELATO/core/config.hpp"

#include <stdexcept>

#ifdef XUSING_ROOT
#include <TDatabasePDG.h>
#endif

GELATO::core::config& GELATO::core::config::fix_system_of_units(system_of_units sys) {
  switch(sys) {
    case system_of_units::GeV_cm_ns:
      pparams.speed_of_light = 29.9791932;
      pparams.hbar = 6.582119569e-16;
      pparams.higgs_vev = 246.22;
      pparams.CKM_Vtd = 8.1e-3;
      pparams.CKM_Vts = 39.4e-3;
      pparams.CKM_Vud = 0.9737;
      pparams.CKM_Vus = 0.224;
      pparams.gFermi = 1.166379e-5;
      pparams.sin2thW = 0.223;
      pparams.pion_decay_constant = 0.130;
      pparams.particles = {
        // name,      pdg,   mass,          lifetime
        {"kaon_pm", { 321, 0.493677,         12.38 }},
        {"kaon_0L", { 130, 0.497648,         51.16 }},
        {"pion_pm", { 211, 0.13957061,       26.03 }},
        {"pion_0",  { 111, 0.1349770,         8.43e-8 }},
        {"elec",    {  11, 0.5109989461e-3,  -1. }}, 
        {"muon",    {  13, 0.1056583745,      2.197e3 }},
        {"gamma",   {  22, 0.,               -1. }},
        {"nu_e",    {  12, 0.,               -1. }},
        {"nu_mu",   {  14, 0.,               -1. }},
        {"nu_tau",  {  16, 0.,               -1. }},
        {"top",     {   6, 172.9,             0. }},
      };
      break;
  }
  return *this;
}

const GELATO::core::config::standard_particle&
GELATO::core::config::physical_parameters::find_particle(const std::string& name) const {
  auto p = particles.find(name);
  if(p == particles.end()) {
    throw std::runtime_error("Error! "+name+" not found in list of particles");
  }
  return p->second;
}

