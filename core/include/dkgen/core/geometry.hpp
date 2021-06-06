#ifndef __dkgen_core_geometry_hpp__
#define __dkgen_core_geometry_hpp__

#include "dkgen/core/ordered_list_of_vectors.hpp"

namespace dkgen {
  namespace core {
    class geometry {
      public:

        geometry() = default;
        geometry(vector3 centre, vector3 dims, rotation rot_from_detctor_to_beamline_system = {});

        geometry& set_active_volume_centre_in_beamline_system(vector3 centre);
        geometry& set_active_volume_half_dimensions_in_detector_system(vector3 dims);
        geometry& set_active_volume_rotation_from_detector_to_beamline_system(rotation rot);

        bool is_beamline_vector_in_active_volume(const vector3& v) const;
        bool is_detector_vector_in_active_volume(const vector3& v) const;
        ordered_list_of_vectors get_active_volume_intersections_for_beamline_vectors(
            const vector3& origin, const vector3& direction) const;
        ordered_list_of_vectors get_active_volume_intersections_for_detector_vectors(
            const vector3& origin, const vector3& direction) const;
        vector3 rotate_beamline_vector_to_detector_coordinates(const vector3& vin) const;
        std::pair<vector3,vector3>
          rotate_and_translate_beamline_vector_to_detector_coordinates(const vector3& vpos, const vector3& vdir) const;
        vector3 rotate_and_translate_beamline_vector_to_detector_coordinates(const vector3& vin) const;
        vector3 rotate_detector_vector_to_beamline_coordinates(const vector3& vin) const;
        std::pair<vector3, vector3>
          rotate_and_translate_detector_vector_to_beamline_coordinates(const vector3& vpos, const vector3& vdir) const;
        vector3 rotate_and_translate_detector_vector_to_beamline_coordinates(const vector3& vin) const;
      private:
        vector3 active_volume_centre_in_beamline_system;
        vector3 active_volume_half_dimensions_in_detector_system;
        rotation active_volume_rotation_from_detector_to_beamline_system;
        rotation active_volume_rotation_from_beamline_to_detector_system;
        void make_active_volume_rotation_from_beamline_to_detector_system();
    };
  }
}

#endif
