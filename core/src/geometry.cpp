#include "geometry.hpp"

#include <utility>
#include <vector>

#undef DEBUG
#ifdef DEBUG
#include <iostream>
#endif

decaygen::geometry::geometry(decaygen::vector3 centre, decaygen::vector3 dims, decaygen::rotation rot_from_detctor_to_beamline_system)
  : active_volume_centre_in_beamline_system{std::move(centre)},
  active_volume_half_dimensions_in_detector_system{std::move(dims)},
  active_volume_rotation_from_detector_to_beamline_system{std::move(rot_from_detctor_to_beamline_system)}
{
  make_active_volume_rotation_from_beamline_to_detector_system();
}

decaygen::geometry& decaygen::geometry::set_active_volume_centre_in_beamline_system(decaygen::vector3 centre) {
  active_volume_centre_in_beamline_system = std::move(centre);
  return *this;
}

decaygen::geometry& decaygen::geometry::set_active_volume_half_dimensions_in_detector_system(decaygen::vector3 dims) {
  active_volume_half_dimensions_in_detector_system = std::move(dims);
  return *this;
}

decaygen::geometry& decaygen::geometry::set_active_volume_rotation_from_detector_to_beamline_system(decaygen::rotation rot) {
  active_volume_rotation_from_detector_to_beamline_system = std::move(rot);
  make_active_volume_rotation_from_beamline_to_detector_system();
  return *this;
}

void decaygen::geometry::make_active_volume_rotation_from_beamline_to_detector_system() {
  active_volume_rotation_from_beamline_to_detector_system = active_volume_rotation_from_detector_to_beamline_system.get_inverse();
}

bool decaygen::geometry::is_beamline_vector_in_active_volume(const vector3& v) const {
  return is_detector_vector_in_active_volume(rotate_and_translate_beamline_vector_to_detector_coordinates(v));
}

bool decaygen::geometry::is_detector_vector_in_active_volume(const vector3& v) const {
  return std::abs(v.x()) < active_volume_half_dimensions_in_detector_system.x()
    && std::abs(v.y()) < active_volume_half_dimensions_in_detector_system.y()
    && std::abs(v.z()) < active_volume_half_dimensions_in_detector_system.z();
}

decaygen::ordered_list_of_vectors decaygen::geometry::get_active_volume_intersections_for_beamline_vectors(
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
  auto transform = [this](decaygen::vector3& v) {
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

decaygen::ordered_list_of_vectors decaygen::geometry::get_active_volume_intersections_for_detector_vectors(
    const vector3& origin, const vector3& direction) const {
  decaygen::ordered_list_of_vectors ret;
  // work out all coordinates that interesect the 6 planes of the bounding box
  // first add case where origin is inside the detector
  if(is_detector_vector_in_active_volume(origin)) {
    ret.add(0.,origin);
  }
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
  return ret;
}

decaygen::vector3 decaygen::geometry::rotate_beamline_vector_to_detector_coordinates(const decaygen::vector3& vin) const {
  return active_volume_rotation_from_beamline_to_detector_system * vin;
}

std::pair<decaygen::vector3, decaygen::vector3>
decaygen::geometry::rotate_and_translate_beamline_vector_to_detector_coordinates(
    const decaygen::vector3& vpos, const decaygen::vector3& vdir) const {
  decaygen::vector3 new_vdir = active_volume_rotation_from_beamline_to_detector_system * vdir;
  decaygen::vector3 new_vpos =
    active_volume_rotation_from_beamline_to_detector_system * (vpos - active_volume_centre_in_beamline_system);
  return std::make_pair<>(std::move(new_vpos), std::move(new_vdir));
}

decaygen::vector3 decaygen::geometry::rotate_and_translate_beamline_vector_to_detector_coordinates(const decaygen::vector3& vin) const {
  return active_volume_rotation_from_beamline_to_detector_system * (vin - active_volume_centre_in_beamline_system);
}

decaygen::vector3 decaygen::geometry::rotate_detector_vector_to_beamline_coordinates(const decaygen::vector3& vin) const {
  return active_volume_rotation_from_detector_to_beamline_system * vin;
}

std::pair<decaygen::vector3, decaygen::vector3>
decaygen::geometry::rotate_and_translate_detector_vector_to_beamline_coordinates(
    const decaygen::vector3& vpos, const decaygen::vector3& vdir) const {
  decaygen::vector3 new_vdir = active_volume_rotation_from_detector_to_beamline_system * vdir;
  decaygen::vector3 new_vpos = active_volume_rotation_from_detector_to_beamline_system * vpos + active_volume_centre_in_beamline_system;
  return std::make_pair<>(std::move(new_vpos), std::move(new_vdir));
}

decaygen::vector3 decaygen::geometry::rotate_and_translate_detector_vector_to_beamline_coordinates(const decaygen::vector3& vin) const {
  return active_volume_rotation_from_detector_to_beamline_system * vin + active_volume_centre_in_beamline_system;
}
