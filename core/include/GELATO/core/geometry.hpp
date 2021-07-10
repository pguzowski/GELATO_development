#ifndef __GELATO_core_geometry_hpp__
#define __GELATO_core_geometry_hpp__

#include "GELATO/core/ordered_list_of_vectors.hpp"

namespace GELATO {
  namespace core {
    class geometry {
      public:

        geometry(vector3 centre, rotation rot_from_detctor_to_beamline_system = {});
        virtual ~geometry() = default;

        geometry& set_active_volume_centre_in_beamline_system(vector3 centre);
        geometry& set_active_volume_rotation_from_detector_to_beamline_system(rotation rot);

        bool is_beamline_vector_in_active_volume(const vector3& v) const;
        virtual bool is_detector_vector_in_active_volume(const vector3& v) const = 0;

        ordered_list_of_vectors get_active_volume_intersections_for_beamline_vectors(
            const vector3& origin, const vector3& direction) const;
        virtual ordered_list_of_vectors get_active_volume_intersections_for_detector_vectors(
            const vector3& origin, const vector3& direction) const = 0;

        virtual ordered_list_of_vectors project_detector_onto_first_vector_along_direction_of_second_vector(
            const vector3& origin1, const vector3& direction1, const vector3& direction2
            ) const = 0;
        
        vector3 rotate_beamline_vector_to_detector_coordinates(const vector3& vin) const;
        std::pair<vector3,vector3>
          rotate_and_translate_beamline_vector_to_detector_coordinates(const vector3& vpos, const vector3& vdir) const;
        vector3 rotate_and_translate_beamline_vector_to_detector_coordinates(const vector3& vin) const;
        
        vector3 rotate_detector_vector_to_beamline_coordinates(const vector3& vin) const;
        std::pair<vector3, vector3>
          rotate_and_translate_detector_vector_to_beamline_coordinates(const vector3& vpos, const vector3& vdir) const;
        vector3 rotate_and_translate_detector_vector_to_beamline_coordinates(const vector3& vin) const;
      protected:
        vector3 active_volume_centre_in_beamline_system;
        rotation active_volume_rotation_from_detector_to_beamline_system;
        rotation active_volume_rotation_from_beamline_to_detector_system;
        void make_active_volume_rotation_from_beamline_to_detector_system();
    };

    class box_geometry : public geometry {
      public:
        box_geometry(vector3 centre, vector3 half_dims, rotation rot_from_detctor_to_beamline_system = {});

        box_geometry& set_active_volume_half_dimensions_in_detector_system(vector3 dims);

        bool is_detector_vector_in_active_volume(const vector3& v) const;

        ordered_list_of_vectors get_active_volume_intersections_for_detector_vectors(
            const vector3& origin, const vector3& direction) const;

        ordered_list_of_vectors project_detector_onto_first_vector_along_direction_of_second_vector(
            const vector3& origin1, const vector3& direction1, const vector3& direction2
            ) const;
      private:
        vector3 active_volume_half_dimensions_in_detector_system;
    };

    
    class sphere_geometry : public geometry {
      public:
        sphere_geometry(vector3 centre, double rad);

        bool is_detector_vector_in_active_volume(const vector3& v) const;

        ordered_list_of_vectors get_active_volume_intersections_for_detector_vectors(
            const vector3& origin, const vector3& direction) const;

        ordered_list_of_vectors project_detector_onto_first_vector_along_direction_of_second_vector(
            const vector3& origin1, const vector3& direction1, const vector3& direction2
            ) const;
      private:
        double radius;
    };

    /*
    // cylinder with length along local z coordinate, and radius in local x-y plane
    class cylinder_geometry : public geometry {
      public:
        cylinder_geometry(vector3 centre, double rad, double len, rotation rot_from_detctor_to_beamline_system = {});

        bool is_detector_vector_in_active_volume(const vector3& v) const;

        ordered_list_of_vectors get_active_volume_intersections_for_detector_vectors(
            const vector3& origin, const vector3& direction) const;

        ordered_list_of_vectors project_detector_onto_first_vector_along_direction_of_second_vector(
            const vector3& origin1, const vector3& direction1, const vector3& direction2
            ) const;
      private:
        double radius;
        double length;
    };
    */
  }
}

#endif
