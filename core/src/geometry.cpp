#include "GELATO/core/geometry.hpp"

#include <utility>
#include <vector>
#include <cmath>

#undef DEBUG
#ifdef DEBUG
#include <iostream>
#endif



///////// BASE class: geometry /////////////////////////////////////////////////////////////////////////////

GELATO::core::geometry::geometry(vector3 centre, rotation rot_from_detctor_to_beamline_system)
  : active_volume_centre_in_beamline_system{std::move(centre)},
  active_volume_rotation_from_detector_to_beamline_system{std::move(rot_from_detctor_to_beamline_system)}
{
  make_active_volume_rotation_from_beamline_to_detector_system();
}

GELATO::core::geometry& GELATO::core::geometry::set_active_volume_centre_in_beamline_system(vector3 centre) {
  active_volume_centre_in_beamline_system = std::move(centre);
  return *this;
}

GELATO::core::geometry& GELATO::core::geometry::set_active_volume_rotation_from_detector_to_beamline_system(rotation rot) {
  active_volume_rotation_from_detector_to_beamline_system = std::move(rot);
  make_active_volume_rotation_from_beamline_to_detector_system();
  return *this;
}

void GELATO::core::geometry::make_active_volume_rotation_from_beamline_to_detector_system() {
  active_volume_rotation_from_beamline_to_detector_system =
#ifdef EXPOSE_PHYSICS_VECTORS
    active_volume_rotation_from_detector_to_beamline_system.inverse();
#else
    active_volume_rotation_from_detector_to_beamline_system.get_inverse();
#endif
}

bool GELATO::core::geometry::is_beamline_vector_in_active_volume(const vector3& v) const {
  return is_detector_vector_in_active_volume(rotate_and_translate_beamline_vector_to_detector_coordinates(v));
}

GELATO::core::ordered_list_of_vectors GELATO::core::geometry::get_active_volume_intersections_for_beamline_vectors(
    const vector3& origin, const vector3& direction) const {
#ifdef DEBUG
  std::cout << __FILE__<<" "<<__LINE__<<" "
    <<origin.x()<<","<<origin.y()<<","<<origin.z()<<" "
    <<direction.x()<<","<<direction.y()<<","<<direction.z()<<std::endl;
#endif
  auto const& vv = rotate_and_translate_beamline_vector_to_detector_coordinates(origin, direction);
#ifdef DEBUG
  std::cout << __FILE__<<" "<<__LINE__<<" "
    <<vv.first.x()<<","<<vv.first.y()<<","<<vv.first.z()<<" "
    <<vv.second.x()<<","<<vv.second.y()<<","<<vv.second.z()<<std::endl;
#endif
  auto r = get_active_volume_intersections_for_detector_vectors(vv.first, vv.second);
  // need to rotate vectors back to beamline coordinates
  auto transform = [this](vector3& v) {
#ifdef DEBUG
  std::cout << __FILE__<<" "<<__LINE__<<" "
    <<v.x()<<","<<v.y()<<","<<v.z();
#endif
    v = rotate_and_translate_detector_vector_to_beamline_coordinates(v);
#ifdef DEBUG
  std::cout <<" --> "
    <<v.x()<<","<<v.y()<<","<<v.z()<<" "<<std::endl;
#endif
  };
  r.apply_transformation(transform);
  return r;
}

GELATO::core::vector3 GELATO::core::geometry::rotate_beamline_vector_to_detector_coordinates(const vector3& vin) const {
  return active_volume_rotation_from_beamline_to_detector_system * vin;
}

std::pair<GELATO::core::vector3, GELATO::core::vector3>
GELATO::core::geometry::rotate_and_translate_beamline_vector_to_detector_coordinates(
    const vector3& vpos, const vector3& vdir) const {
  vector3 new_vdir = active_volume_rotation_from_beamline_to_detector_system * vdir;
  vector3 new_vpos =
    active_volume_rotation_from_beamline_to_detector_system * (vpos - active_volume_centre_in_beamline_system);
  return std::make_pair<>(std::move(new_vpos), std::move(new_vdir));
}

GELATO::core::vector3 GELATO::core::geometry::rotate_and_translate_beamline_vector_to_detector_coordinates(const vector3& vin) const {
  return active_volume_rotation_from_beamline_to_detector_system * (vin - active_volume_centre_in_beamline_system);
}

GELATO::core::vector3 GELATO::core::geometry::rotate_detector_vector_to_beamline_coordinates(const vector3& vin) const {
  return active_volume_rotation_from_detector_to_beamline_system * vin;
}

std::pair<GELATO::core::vector3, GELATO::core::vector3>
GELATO::core::geometry::rotate_and_translate_detector_vector_to_beamline_coordinates(
    const vector3& vpos, const vector3& vdir) const {
  vector3 new_vdir = active_volume_rotation_from_detector_to_beamline_system * vdir;
  vector3 new_vpos = active_volume_rotation_from_detector_to_beamline_system * vpos + active_volume_centre_in_beamline_system;
  return std::make_pair<>(std::move(new_vpos), std::move(new_vdir));
}

GELATO::core::vector3
GELATO::core::geometry::rotate_and_translate_detector_vector_to_beamline_coordinates(const vector3& vin) const {
  return active_volume_rotation_from_detector_to_beamline_system * vin + active_volume_centre_in_beamline_system;
}



///////// DERIVED class: box_geometry //////////////////////////////////////////////////////////////////////

GELATO::core::box_geometry::box_geometry(vector3 centre, vector3 half_dims, rotation rot_from_detctor_to_beamline_system)
  : GELATO::core::geometry{centre, rot_from_detctor_to_beamline_system},
  active_volume_half_dimensions_in_detector_system{std::move(half_dims)}
{
}

GELATO::core::box_geometry& GELATO::core::box_geometry::set_active_volume_half_dimensions_in_detector_system(vector3 dims) {
  active_volume_half_dimensions_in_detector_system = std::move(dims);
  return *this;
}

bool GELATO::core::box_geometry::is_detector_vector_in_active_volume(const vector3& v) const {
  return std::abs(v.x()) < active_volume_half_dimensions_in_detector_system.x()
    && std::abs(v.y()) < active_volume_half_dimensions_in_detector_system.y()
    && std::abs(v.z()) < active_volume_half_dimensions_in_detector_system.z();
}

GELATO::core::ordered_list_of_vectors GELATO::core::box_geometry::get_active_volume_intersections_for_detector_vectors(
    const vector3& origin, const vector3& direction) const {
  ordered_list_of_vectors ret;
  // work out all coordinates that interesect the 6 planes of the bounding box
  // first add case where origin is inside the detector
  if(is_detector_vector_in_active_volume(origin)) {
    ret.add(0.,origin);
    if(!(0. < direction.mag())) {
      ret.add(0.,origin);
    }
  }
  if(0. < direction.mag()) {
    auto get_plane_intersections = [](
        double direction_coordinate, double plane_coordinate, double origin_coordinate,
        double other_origin_coordinate_1, double other_origin_coordinate_2,
        double other_direction_coordinate_1, double other_direction_coordinate_2,
        double dimension_1, double dimension_2
        ) {
      std::vector<double> alphas;
      // if parallel to the plane, there is no intersection
      if(std::abs(direction_coordinate) > 0.) {
        // check both sides
        for(int plane_offset: {-1, +1}) {
          const double alpha = (plane_offset * plane_coordinate - origin_coordinate) / direction_coordinate;
          if(alpha < 0.) {
            // going in the wrong direction
            continue;
          }
          const double intersection_coordinate_1 = alpha * other_direction_coordinate_1 + other_origin_coordinate_1;
          const double intersection_coordinate_2 = alpha * other_direction_coordinate_2 + other_origin_coordinate_2;
          if(std::abs(intersection_coordinate_1) < dimension_1 && std::abs(intersection_coordinate_2) < dimension_2) {
            alphas.push_back(alpha);
          }
        }
      }
      return alphas;
    };
    for(auto plane: {'x', 'y', 'z'}) {
      switch(plane) {
        case 'x':
          {
            auto alphas = get_plane_intersections(direction.x(), active_volume_half_dimensions_in_detector_system.x(), origin.x(),
                origin.y(), origin.z(),
                direction.y(), direction.z(),
                active_volume_half_dimensions_in_detector_system.y(), active_volume_half_dimensions_in_detector_system.z());
            for(auto a : alphas) {
              ret.add(a, a * direction + origin);
            }
          }
          break;
        case 'y':
          {
            auto alphas = get_plane_intersections(direction.y(), active_volume_half_dimensions_in_detector_system.y(), origin.y(),
                origin.x(), origin.z(),
                direction.x(), direction.z(),
                active_volume_half_dimensions_in_detector_system.x(), active_volume_half_dimensions_in_detector_system.z());
            for(auto a : alphas) {
              ret.add(a, a * direction + origin);
            }
          }
          break;
        case 'z':
          {
            auto alphas = get_plane_intersections(direction.z(), active_volume_half_dimensions_in_detector_system.z(), origin.z(),
                origin.x(), origin.y(),
                direction.x(), direction.y(),
                active_volume_half_dimensions_in_detector_system.x(), active_volume_half_dimensions_in_detector_system.y());
            for(auto a : alphas) {
              ret.add(a, a * direction + origin);
            }
          }
          break;
      }
    }
  }
  return ret;
}

GELATO::core::ordered_list_of_vectors
GELATO::core::box_geometry::project_detector_onto_first_vector_along_direction_of_second_vector(
  const vector3& origin1, const vector3& direction1, const vector3& direction2
) const {
  ordered_list_of_vectors ret;
  const double dotp = direction1.unit().dot(direction2.unit());

  if(dotp <= -1.) {
    // direction2 is antiparallel to direction1; return empty list
    return ret;
  }
  if(dotp >= 1.) {
    // direction2 is parallel to direction1; no plane possible
    // generate from 0 to max perpendicular distance
    ret.add(0.,origin1);
  }
  bool has_origin = false;
  for(auto sideX: {-1., 1.}) {
    for(auto sideY: {-1., 1.}) {
      for(auto sideZ: {-1., 1.}) {
        const vector3& corner = rotate_and_translate_detector_vector_to_beamline_coordinates({
            sideX * active_volume_half_dimensions_in_detector_system.x(),
            sideY * active_volume_half_dimensions_in_detector_system.y(),
            sideZ * active_volume_half_dimensions_in_detector_system.z(),
            });
        const vector3& corner_from_origin = corner - origin1;
        if(dotp > -1. && dotp < 1.) {
          const double dist = (corner_from_origin.dot(direction2.unit()) - corner_from_origin.dot(direction1.unit()))
            / (dotp*dotp - 1.);
          if(dist > 0.) {
            ret.add(dist, origin1 + dist * direction1.unit());
          }
          else if(!has_origin) {
            ret.add(0., origin1);
            has_origin = true;
          }
        }
        else {
          // direction2 is parallel to direction1; no plane possible
          // use perpendicular projection onto vector1
          const double dist = corner_from_origin.dot(direction1.unit());
          if(dist > 0.) {
            ret.add(dist, origin1 + dist * direction1.unit());
          }
        }
      }
    }
  }
  return ret;
}



///////// DERIVED class: sphere_geometry ///////////////////////////////////////////////////////////////////


///////// DERIVED class: cylinder_geometry /////////////////////////////////////////////////////////////////

