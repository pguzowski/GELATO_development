#ifndef __GELATO_core_vectors_hpp__
#define __GELATO_core_vectors_hpp__

#ifdef EXPOSE_PHYSICS_VECTORS

#ifndef USING_CLHEP
#error "Need to set USING_CLHEP flag"
#endif

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/Rotation.h>
namespace GELATO {
  namespace core {
    using vector3 = CLHEP::Hep3Vector;
    using fourvector = CLHEP::HepLorentzVector;
    using rotation = CLHEP::HepRotation;
  }
}

#else // EXPOSE_PHYSICS_VECTORS


// hide away the physics implementation from the end user
// allows compiling with different library implementations
// without polluting end user namespace with included libraries

#include <memory>
//#include <experimental/propagate_const>

#include "GELATO/core/physics_vector_sizes.hpp"

namespace GELATO {
  namespace core {
    class vector3 {
      public:
        vector3();
        vector3(double x, double y, double z);
        ~vector3();
        vector3(const vector3&);
        vector3& operator=(const vector3&);
        /* No need for move semantics (all data stored locally in impl_storage)
        vector3(vector3&&);
        vector3& operator=(vector3&&);
        */

        vector3 unit() const;
        vector3 operator*(double m) const;
        double mag() const;
        vector3& set_mag(double m);

        double x() const;
        double y() const;
        double z() const;
        vector3& set_xyz(double x, double y, double z);

        vector3 operator-() const; // unary negation
        vector3 operator-(const vector3& v) const;
        vector3 operator+(const vector3& v) const;

        double dot(const vector3& v) const;
        double delta_phi(const vector3& v) const;

        vector3 cross(const vector3& v) const;
      private:
        struct _vector3_impl_;
        //using impl_ptr = std::experimental::propagate_const<std::unique_ptr<_vector3_impl_>>;
        //using impl_ptr = std::experimental::propagate_const<_vector3_impl_*>;
        //impl_ptr impl;
        static constexpr size_t c_impl_size = vectors::c_vector3_impl_size;
        static constexpr size_t c_impl_alignment = std::alignment_of<double>::value;
        typedef std::aligned_storage<c_impl_size, c_impl_alignment>::type aligned_storage_type;
        aligned_storage_type impl_storage;
        friend class fourvector;
        friend class rotation;
        _vector3_impl_* impl;
        //inline const _vector3_impl_* impl() const { return (_vector3_impl_*)&impl_storage; }
        //inline _vector3_impl_* impl() { return (_vector3_impl_*)&impl_storage; }
    };

    class fourvector {
      public:
        fourvector();
        ~fourvector();
        fourvector(double x, double y, double z, double t);
        fourvector(const fourvector&);
        fourvector& operator=(const fourvector&);
        /* No need for move semantics (all data stored locally in impl_storage)
        fourvector(fourvector&&);
        fourvector& operator=(fourvector&&);
        */

        double x() const;
        double y() const;
        double z() const;
        double t() const;
        fourvector& set_xyzt(double x, double y, double z, double t);
        fourvector& set_t(double t);

        double px() const;
        double py() const;
        double pz() const;
        double e() const;
        fourvector& set_mass_momentum_theta_phi(double m, double mom, double t, double p);

        double m() const;
        double m2() const;

        vector3 vect() const;
        fourvector& set_vec_time(const vector3& vec, double time);

        fourvector& boost(const vector3& boost_vector);
        vector3 get_boost_vector() const;

        double beta() const;
        double gamma() const;

        fourvector operator+(const fourvector& v2) const;
        fourvector operator-(const fourvector& v2) const;
      private:
        struct _fourvector_impl_;
        //using impl_ptr = std::experimental::propagate_const<std::unique_ptr<_fourvector_impl_>>;
        //using impl_ptr = std::experimental::propagate_const<_fourvector_impl_*>;
        //impl_ptr impl;
        static constexpr size_t c_impl_size = vectors::c_fourvector_impl_size;
        static constexpr size_t c_impl_alignment = std::alignment_of<double>::value;
        typedef std::aligned_storage<c_impl_size, c_impl_alignment>::type aligned_storage_type;
        aligned_storage_type impl_storage;
        _fourvector_impl_ *impl;
        //inline const _fourvector_impl_* impl() const { return reinterpret_cast<const _fourvector_impl_*>(&impl_storage); }
        //inline _fourvector_impl_* impl() { return reinterpret_cast<_fourvector_impl_*>(&impl_storage); }
    };

    class rotation {
      public:
        rotation();
        ~rotation();
        rotation(const rotation&);
        rotation& operator=(const rotation&);
        /* No need for move semantics (all data stored locally in impl_storage)
        rotation& operator=(rotation&&);
        rotation(rotation&&);
        */

        vector3 operator*(const vector3& v) const;

        rotation& rotate_axes(const vector3& new_x_axis, const vector3& new_y_axis, const vector3& new_z_axis);
        rotation& rotate(double angle, const vector3& axis);

        rotation& invert();
        rotation get_inverse() const;
      private:
        struct _rotation_impl_;
        //using impl_ptr = std::experimental::propagate_const<std::unique_ptr<_rotation_impl_>>;
        //using impl_ptr = std::experimental::propagate_const<_rotation_impl_*>;
        //impl_ptr impl;
        static constexpr size_t c_impl_size = vectors::c_rotation_impl_size;
        static constexpr size_t c_impl_alignment = std::alignment_of<double>::value;
        typedef std::aligned_storage<c_impl_size, c_impl_alignment>::type aligned_storage_type;
        aligned_storage_type impl_storage;
        _rotation_impl_ *impl;
        //inline const _rotation_impl_* impl() const { return reinterpret_cast<const _rotation_impl_*>(&impl_storage); }
        //inline _rotation_impl_* impl() { return reinterpret_cast<_rotation_impl_*>(&impl_storage); }
    };

  }
}

// needs to be outside namespaces to overload global operator*
inline GELATO::core::vector3 operator*(double a, const GELATO::core::vector3& b) {
  return b*a;
}
#endif // EXPOSE_PHYSICS_VECTORS



#endif // __GELATO_core_vectors_hpp__

