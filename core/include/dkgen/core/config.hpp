#ifndef __dkgen_core_config_hpp__
#define __dkgen_core_config_hpp__

#include <string>
#include <map>

namespace dkgen {
  namespace core {
    class config {
      public:
        enum class system_of_units { GeV_cm_ns };
        struct standard_particle {
          int pdgcode;
          double mass;
          double lifetime;
        };
        struct physical_parameters {
          double speed_of_light;
          double hbar;
          double higgs_vev;
          double CKM_Vtd;
          double CKM_Vts;
          physical_parameters() : speed_of_light(-1.) {}
          const standard_particle& find_particle(const std::string& name) const;
          private:
          std::map<std::string, standard_particle> particles;
          friend class config;
        };
        const physical_parameters& physical_params() const { return pparams; }
        bool force_decays_inside_detector() const { return force_decays_in_detector; }
        config& fix_system_of_units(system_of_units sys);
        config& set_force_decays_in_detector(bool fdid) { force_decays_in_detector = fdid; return *this; }
      private:
        bool force_decays_in_detector = false;
        physical_parameters pparams;
    };
  }
}

#endif
