#ifndef __dkgen_core_vectors_hpp__
#define __dkgen_core_vectors_hpp__

#ifdef EXPOSE_PHYSICS_VECTORS

#ifndef USING_CLHEP
#error "Need to set USE_CLHEP flag"
#endif

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/Rotation.h>
namespace dkgen {
  namespace core {
    using vector3 = CLHEP::Hep3Vector;
    using fourvector = CLHEP::HepLorentzVector;
//#define set_t setT
//#define set_vec_time(v,t) set(v,t) 
//#define set_mass_momentum_theta_phi(m,mm,t,p) setVectM(vector3{mm*std::sin(t)*std::cos(p),mm*std::sin(t)*std::sin(p),mm*std::cos(t)},m)
//#define get_boost_vector boostVector
//#define get_inverse inverse
    using rotation = CLHEP::HepRotation;
  }
}

#else // EXPOSE_PHYSICS_VECTORS


// hide away the physics implementation from the end user
// allows compiling with different library implementations
// without polluting end user namespace with included libraries
// However performance is 3x worse due to memory allocations
// need to investigate speedups

#include <memory>
#include <experimental/propagate_const>

namespace dkgen {
  namespace core {
    struct _vector3_impl_;
    class vector3 {
      public:
        using impl_ptr = std::experimental::propagate_const<std::unique_ptr<_vector3_impl_>>;

        vector3();
        vector3(double x, double y, double z);
        ~vector3();
        vector3(const vector3&);
        vector3(vector3&&);
        vector3& operator=(const vector3&);
        vector3& operator=(vector3&&);

        vector3 unit() const;
        vector3 operator*(double m) const;
        double mag() const;
        vector3& set_mag(double m);

        double x() const;
        double y() const;
        double z() const;
        vector3& set_xyz(double x, double y, double z);

        vector3 operator-(const vector3& v) const;
        vector3 operator+(const vector3& v) const;
      private:
        impl_ptr impl;
        friend class fourvector;
        friend class rotation;
    };

    struct _fourvector_impl_;
    class fourvector {
      public:
        using impl_ptr = std::experimental::propagate_const<std::unique_ptr<_fourvector_impl_>>;

        fourvector();
        ~fourvector();
        fourvector(double x, double y, double z, double t);
        fourvector(const fourvector&);
        fourvector(fourvector&&);
        fourvector& operator=(const fourvector&);
        fourvector& operator=(fourvector&&);

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

        fourvector operator+(const fourvector& v2) const;
        fourvector operator-(const fourvector& v2) const;
      private:
        impl_ptr impl;
    };

    struct _rotation_impl_;
    class rotation {
      public:
        using impl_ptr = std::experimental::propagate_const<std::unique_ptr<_rotation_impl_>>;

        rotation();
        ~rotation();
        rotation(const rotation&);
        rotation(rotation&&);
        rotation& operator=(const rotation&);
        rotation& operator=(rotation&&);

        vector3 operator*(const vector3& v) const;

        rotation& rotate_axes(const vector3& new_x_axis, const vector3& new_y_axis, const vector3& new_z_axis);

        rotation& invert();
        rotation get_inverse() const;
      private:
        impl_ptr impl;
    };

  }
}

// needs to be outside namespaces to overload global operator*
dkgen::core::vector3 operator*(double a, const dkgen::core::vector3& b);
#endif // EXPOSE_PHYSICS_VECTORS



#endif // __dkgen_core_vectors_hpp__

