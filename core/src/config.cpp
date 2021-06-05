#include "dkgen/core/config.hpp"

#include <stdexcept>
#include <limits>

dkgen::core::config& dkgen::core::config::fix_system_of_units(dkgen::core::config::system_of_units sys) {
  switch(sys) {
    case GeV_cm_ns:
      pparams.speed_of_light = 29.9791932;
      pparams.hbar = 6.582119569e-16;
      pparams.higgs_vev = 246.22;
      pparams.CKM_Vtd = 8.1e-3;
      pparams.CKM_Vts = 39.4e-3;
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

