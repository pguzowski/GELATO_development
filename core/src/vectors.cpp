#include "dkgen/core/vectors.hpp"

// need to define USING_CLHEP or USING_ROOT or USING_XXXXX (tbd) on the compile command line

#ifdef USING_CLHEP
#if defined(USING_ROOT) || defined(USING_XXXX)
#error can only have one of USING_CLHEP or USING_ROOT or USING_XXXX set
#endif
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/Rotation.h>
using VECTOR3_CLASS = CLHEP::Hep3Vector;
using VECTOR4_CLASS = CLHEP::HepLorentzVector;
using ROTATION_CLASS = CLHEP::HepRotation;
#define PHYSICS_IS_SET
#endif

#ifdef USING_ROOT
#if defined(USING_CLHEP) || defined(USING_XXXX)
#error can only have one of USING_CLHEP or USING_ROOT or USING_XXXX set
#endif
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRotation.h>
using VECTOR3_CLASS = TVector3;
using VECTOR4_CLASS = TLorentzVector;
using ROTATION_CLASS = TRotation;
#define PHYSICS_IS_SET
#endif


#ifdef USING_XXXX
#error "XXXX" is unimplemented
#endif


#ifndef PHYSICS_IS_SET
#error must define exactly one of USING_CLHEP, USING_ROOT, or USING_XXXX
#endif



namespace dkgen {
  namespace core {
    struct _vector3_impl_ {
      VECTOR3_CLASS v;
    };
    struct _fourvector_impl_ {
      VECTOR4_CLASS v;
    };
    struct _rotation_impl_ {
      ROTATION_CLASS r;
    };
  }
}


dkgen::core::vector3::vector3() : impl{std::make_unique<_vector3_impl_>()} {}

dkgen::core::vector3::vector3(double x, double y, double z) : impl{std::make_unique<_vector3_impl_>()} {
  set_xyz(x,y,z);
}

dkgen::core::vector3::~vector3() { }

dkgen::core::vector3::vector3(const vector3& v) : impl{std::make_unique<_vector3_impl_>(*(v.impl))} {}

dkgen::core::vector3::vector3(vector3&& v) : impl{std::move(v.impl)} {
}

dkgen::core::vector3& dkgen::core::vector3::operator=(const dkgen::core::vector3& v) {
  impl->v = v.impl->v;
  return *this;
}

dkgen::core::vector3& dkgen::core::vector3::operator=(dkgen::core::vector3&& v) {
  impl = std::move(v.impl);
  return *this;
}

dkgen::core::vector3 dkgen::core::vector3::unit() const {
  if(mag() <= 0.) {
    throw std::runtime_error("vector3: unit() requested for zero or negative length vector");
  }
  vector3 ret(*this);
  ret.set_mag(1.);
  return ret;
}

dkgen::core::vector3 dkgen::core::vector3::operator*(double a) const {
  vector3 ret(*this);
  ret.set_mag(ret.mag()*a);
  return ret;
}

dkgen::core::vector3 operator*(double a, const dkgen::core::vector3& b) {
  return b*a;
}

double dkgen::core::vector3::mag() const {
#ifdef USING_CLHEP
  return impl->v.mag();
#elif defined(USING_ROOT)
  return impl->v.Mag();
#endif
}

dkgen::core::vector3& dkgen::core::vector3::set_mag(double m) {
#ifdef USING_CLHEP
  impl->v.setMag(m);
#elif defined(USING_ROOT)
  impl->v.SetMag(m);
#endif
  return *this;
}

dkgen::core::vector3& dkgen::core::vector3::set_xyz(double x, double y, double z) {
#ifdef USING_CLHEP
  impl->v.set(x,y,z);
#elif defined(USING_ROOT)
  impl->v.SetXYZ(x,y,z);
#endif
  return *this;
}

double dkgen::core::vector3::x() const {
#ifdef USING_CLHEP
  return impl->v.x();
#elif defined USING_ROOT
  return impl->v.X();
#endif
}

double dkgen::core::vector3::y() const {
#ifdef USING_CLHEP
  return impl->v.y();
#elif defined USING_ROOT
  return impl->v.Y();
#endif
}

double dkgen::core::vector3::z() const {
#ifdef USING_CLHEP
  return impl->v.z();
#elif defined USING_ROOT
  return impl->v.Z();
#endif
}

dkgen::core::vector3 dkgen::core::vector3::operator-(const dkgen::core::vector3& v2) const {
  vector3 ret(*this);
  ret.impl->v -= v2.impl->v;
  return ret;
}

dkgen::core::vector3 dkgen::core::vector3::operator+(const dkgen::core::vector3& v2) const {
  vector3 ret(*this);
  ret.impl->v += v2.impl->v;
  return ret;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////





dkgen::core::fourvector::fourvector() : impl{std::make_unique<_fourvector_impl_>()} {}

dkgen::core::fourvector::fourvector(double x, double y, double z, double t) : impl{std::make_unique<_fourvector_impl_>()} {
  set_xyzt(x,y,z,t);
}

dkgen::core::fourvector::~fourvector() { }

dkgen::core::fourvector::fourvector(const dkgen::core::fourvector& v) : impl{std::make_unique<_fourvector_impl_>(*(v.impl))} {}

dkgen::core::fourvector::fourvector(dkgen::core::fourvector&& v) : impl{std::move(v.impl)} {
}

dkgen::core::fourvector& dkgen::core::fourvector::operator=(const dkgen::core::fourvector& v) {
  impl->v = v.impl->v;
  return *this;
}

dkgen::core::fourvector& dkgen::core::fourvector::operator=(dkgen::core::fourvector&& v) {
  impl = std::move(v.impl);
  return *this;
}

dkgen::core::fourvector& dkgen::core::fourvector::set_xyzt(double x, double y, double z, double t) {
#ifdef USING_CLHEP
  impl->v.set(x,y,z,t);
#elif defined USING_ROOT
  impl->v.SetXYZT(x,y,z,t);
#endif
  return *this;
}

dkgen::core::fourvector& dkgen::core::fourvector::set_t(double t) {
#ifdef USING_CLHEP
  impl->v.setT(t);
#elif defined USING_ROOT
  impl->v.SetT(t);
#endif
  return *this;
}

double dkgen::core::fourvector::x() const {
#ifdef USING_CLHEP
  return impl->v.x();
#elif defined USING_ROOT
  return impl->v.X();
#endif
}

double dkgen::core::fourvector::y() const {
#ifdef USING_CLHEP
  return impl->v.y();
#elif defined USING_ROOT
  return impl->v.Y();
#endif
}

double dkgen::core::fourvector::z() const {
#ifdef USING_CLHEP
  return impl->v.z();
#elif defined USING_ROOT
  return impl->v.Z();
#endif
}

double dkgen::core::fourvector::t() const {
#ifdef USING_CLHEP
  return impl->v.t();
#elif defined USING_ROOT
  return impl->v.T();
#endif
}

double dkgen::core::fourvector::px() const {
#ifdef USING_CLHEP
  return impl->v.px();
#elif defined USING_ROOT
  return impl->v.Px();
#endif
}

double dkgen::core::fourvector::py() const {
#ifdef USING_CLHEP
  return impl->v.py();
#elif defined USING_ROOT
  return impl->v.Py();
#endif
}

double dkgen::core::fourvector::pz() const {
#ifdef USING_CLHEP
  return impl->v.pz();
#elif defined USING_ROOT
  return impl->v.Pz();
#endif
}

double dkgen::core::fourvector::e() const {
#ifdef USING_CLHEP
  return impl->v.e();
#elif defined USING_ROOT
  return impl->v.E();
#endif
}

dkgen::core::fourvector& dkgen::core::fourvector::set_mass_momentum_theta_phi(double m, double mom, double t, double p) {
#ifdef USING_CLHEP
  impl->v.setVectM({mom*std::sin(t)*std::cos(p), mom*std::sin(t)*std::sin(p), mom*std::cos(t)},m);
#elif defined USING_ROOT
  impl->v.SetVectM({mom*std::sin(t)*std::cos(p), mom*std::sin(t)*std::sin(p), mom*std::cos(t)},m);
#endif
  return *this;
}

double dkgen::core::fourvector::m() const {
#ifdef USING_CLHEP
  return impl->v.m();
#elif defined USING_ROOT
  return impl->v.M();
#endif
}

double dkgen::core::fourvector::m2() const {
#ifdef USING_CLHEP
  return impl->v.m2();
#elif defined USING_ROOT
  return impl->v.M2();
#endif
}

dkgen::core::vector3 dkgen::core::fourvector::vect() const {
  vector3 ret;
#ifdef USING_CLHEP
  ret.impl->v = impl->v.vect();
#elif defined USING_ROOT
  ret.impl->v = impl->v.Vect();
#endif
  return ret;
}

dkgen::core::fourvector& dkgen::core::fourvector::set_vec_time(const dkgen::core::vector3& vec, double time)  {
#ifdef USING_CLHEP
  impl->v.setVect(vec.impl->v);
  impl->v.setT(time);
#elif defined USING_ROOT
  impl->v.SetVect(vec.impl->v);
  impl->v.SetT(time);
#endif
  return *this;
}


dkgen::core::vector3 dkgen::core::fourvector::get_boost_vector() const {
  vector3 ret;
#ifdef USING_CLHEP
  ret.impl->v = impl->v.boostVector();
#elif defined USING_ROOT
  ret.impl->v = impl->v.BoostVector();
#endif
  return ret;
}

dkgen::core::fourvector& dkgen::core::fourvector::boost(const dkgen::core::vector3& bv) {
#ifdef USING_CLHEP
  impl->v.boost(bv.impl->v);
#elif defined USING_ROOT
  impl->v.Boost(bv.impl->v);
#endif
  return *this;
}

double dkgen::core::fourvector::beta() const {
#ifdef USING_CLHEP
  return impl->v.beta();
#elif defined USING_ROOT
  return impl->v.Beta();
#endif
}


dkgen::core::fourvector dkgen::core::fourvector::operator-(const dkgen::core::fourvector& v2) const {
  fourvector ret(*this);
  ret.impl->v -= v2.impl->v;
  return ret;
}

dkgen::core::fourvector dkgen::core::fourvector::operator+(const dkgen::core::fourvector& v2) const {
  fourvector ret(*this);
  ret.impl->v += v2.impl->v;
  return ret;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







dkgen::core::rotation::rotation() : impl{std::make_unique<_rotation_impl_>()} {}

dkgen::core::rotation::~rotation() { }

dkgen::core::rotation::rotation(const dkgen::core::rotation& r) : impl{std::make_unique<_rotation_impl_>(*(r.impl))} {}

dkgen::core::rotation::rotation(dkgen::core::rotation&& r) : impl{std::move(r.impl)} {
}

dkgen::core::rotation& dkgen::core::rotation::operator=(const dkgen::core::rotation& r) {
  impl->r = r.impl->r;
  return *this;
}

dkgen::core::rotation& dkgen::core::rotation::operator=(dkgen::core::rotation&& r) {
  impl = std::move(r.impl);
  return *this;
}

dkgen::core::vector3 dkgen::core::rotation::operator*(const dkgen::core::vector3& vin) const {
  dkgen::core::vector3 ret;
  ret.impl->v = impl->r * vin.impl->v;
  return ret;
}

dkgen::core::rotation& dkgen::core::rotation::invert() { 
#ifdef USING_CLHEP
  impl->r.invert();
#elif defined USING_ROOT
  impl->r.Invert();
#endif
  return *this;
}

dkgen::core::rotation dkgen::core::rotation::get_inverse()  const { 
  rotation ret;
#ifdef USING_CLHEP
  ret.impl->r = impl->r.inverse();
#elif defined USING_ROOT
  ret.impl->r = impl->r.Inverse();
#endif
  return ret;
}

dkgen::core::rotation& dkgen::core::rotation::rotate_axes(const dkgen::core::vector3& new_x, const dkgen::core::vector3& new_y, const dkgen::core::vector3& new_z) {
#ifdef USING_CLHEP
  impl->r.rotateAxes(new_x.impl->v, new_y.impl->v, new_z.impl->v);
#elif defined USING_ROOT
  impl->r.RotateAxes(new_x.impl->v, new_y.impl->v, new_z.impl->v);
#endif
  return *this;
}



#undef USING_CLHEP
#undef USING_ROOT
#undef USING_XXXX
