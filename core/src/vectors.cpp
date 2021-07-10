#ifdef EXPOSE_PHYSICS_VECTORS

// nothing to do here

#else // EXPOSE_PHYSICS_VECTORS
#include "GELATO/core/vectors.hpp"

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



namespace GELATO {
  namespace core {
    struct vector3::_vector3_impl_ : public VECTOR3_CLASS {
      public:
        //using VECTOR3_CLASS::operator=;
        _vector3_impl_& operator=(const VECTOR3_CLASS& v) { 
          VECTOR3_CLASS::operator=(v);
          return *this;
        }
    };
    struct fourvector::_fourvector_impl_ : public VECTOR4_CLASS {
      public:
    };
    struct rotation::_rotation_impl_ : public ROTATION_CLASS {
      public:
        //using ROTATION_CLASS::operator=;
        _rotation_impl_& operator=(const ROTATION_CLASS& v) { 
          ROTATION_CLASS::operator=(v);
          return *this;
        }
    };
  }
}


GELATO::core::vector3::vector3() {
  static_assert(sizeof(_vector3_impl_) <= c_impl_size, "Size of vector3 implementation too big");
  static_assert(alignof(_vector3_impl_) == c_impl_alignment, "Wrong alignment of vector3 implementation");
  new (&impl_storage) _vector3_impl_;
  impl = reinterpret_cast<_vector3_impl_*>(&impl_storage);
}

GELATO::core::vector3::vector3(double x, double y, double z)  {
  static_assert(sizeof(_vector3_impl_) <= c_impl_size, "Size of vector3 implementation too big");
  static_assert(alignof(_vector3_impl_) == c_impl_alignment, "Wrong alignment of vector3 implementation");
  new (&impl_storage) _vector3_impl_ ;
  impl = reinterpret_cast<_vector3_impl_*>(&impl_storage);
  set_xyz(x,y,z);
}

GELATO::core::vector3::~vector3() {
  impl->~_vector3_impl_();
}

GELATO::core::vector3::vector3(const vector3& v)  {
  static_assert(sizeof(_vector3_impl_) <= c_impl_size, "Size of vector3 implementation too big");
  static_assert(alignof(_vector3_impl_) == c_impl_alignment, "Wrong alignment of vector3 implementation");
  new (&impl_storage) _vector3_impl_(*(v.impl));
  impl = reinterpret_cast<_vector3_impl_*>(&impl_storage);
}

GELATO::core::vector3& GELATO::core::vector3::operator=(const vector3& v) {
  *impl = *v.impl;
  return *this;
}

/* No need for move semantics (all data stored locally in impl_storage)

// local storage, just copy
GELATO::core::vector3::vector3(vector3&& v) {
  static_assert(sizeof(_vector3_impl_) <= c_impl_size, "Size of vector3 implementation too big");
  static_assert(alignof(_vector3_impl_) == c_impl_alignment, "Wrong alignment of vector3 implementation");
  new (&impl_storage) _vector3_impl_(*(v.impl)) ;
  impl = reinterpret_cast<_vector3_impl_*>(&impl_storage);
}
// local storage, just copy
GELATO::core::vector3& GELATO::core::vector3::operator=(vector3&& v) {
  *impl = *v.impl;
  return *this;
}

*/

GELATO::core::vector3 GELATO::core::vector3::unit() const {
  /*
  if(mag() <= 0.) {
    throw std::runtime_error("vector3: unit() requested for zero or negative length vector");
  }
  */
  vector3 ret(*this);
  ret.set_mag(1.);
  return ret;
}

GELATO::core::vector3 GELATO::core::vector3::operator*(double a) const {
  vector3 ret(*this);
  if(ret.mag() > 0.) {
    ret.set_mag(ret.mag()*a);
  }
  return ret;
}

/*
GELATO::core::vector3 operator*(double a, const GELATO::core::vector3& b) {
  return b*a;
}
*/

double GELATO::core::vector3::mag() const {
#ifdef USING_CLHEP
  return impl->mag();
#elif defined(USING_ROOT)
  return impl->Mag();
#endif
}

GELATO::core::vector3& GELATO::core::vector3::set_mag(double m) {
#ifdef USING_CLHEP
  impl->setMag(m);
#elif defined(USING_ROOT)
  impl->SetMag(m);
#endif
  return *this;
}

GELATO::core::vector3& GELATO::core::vector3::set_xyz(double x, double y, double z) {
#ifdef USING_CLHEP
  impl->set(x,y,z);
#elif defined(USING_ROOT)
  impl->SetXYZ(x,y,z);
#endif
  return *this;
}

double GELATO::core::vector3::x() const {
#ifdef USING_CLHEP
  return impl->x();
#elif defined USING_ROOT
  return impl->X();
#endif
}

double GELATO::core::vector3::y() const {
#ifdef USING_CLHEP
  return impl->y();
#elif defined USING_ROOT
  return impl->Y();
#endif
}

double GELATO::core::vector3::z() const {
#ifdef USING_CLHEP
  return impl->z();
#elif defined USING_ROOT
  return impl->Z();
#endif
}

GELATO::core::vector3 GELATO::core::vector3::operator-() const {
  vector3 ret(*this);
  *ret.impl = -*ret.impl;
  return ret;
}

GELATO::core::vector3 GELATO::core::vector3::operator-(const vector3& v2) const {
  vector3 ret(*this);
  *ret.impl -= *v2.impl;
  return ret;
}

GELATO::core::vector3 GELATO::core::vector3::operator+(const vector3& v2) const {
  vector3 ret(*this);
  *ret.impl += *v2.impl;
  return ret;
}

double GELATO::core::vector3::dot(const vector3& v) const {
#ifdef USING_CLHEP
  return impl->dot(*v.impl);
#elif defined USING_ROOT
  return impl->Dot(*v.impl);
#endif
}

double GELATO::core::vector3::delta_phi(const vector3& v) const {
#ifdef USING_CLHEP
  return impl->deltaPhi(*v.impl);
#elif defined USING_ROOT
  return impl->DeltaPhi(*v.impl);
#endif
}

GELATO::core::vector3 GELATO::core::vector3::cross(const vector3& v) const {
  vector3 ret(*this);
#ifdef USING_CLHEP
  ret.impl->cross(*v.impl);
#elif defined USING_ROOT
  ret.impl->Cross(*v.impl);
#endif
  return ret;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////





GELATO::core::fourvector::fourvector() {
  static_assert(sizeof(_fourvector_impl_) <= c_impl_size, "Size of fourvector implementation too big");
  static_assert(alignof(_fourvector_impl_) == c_impl_alignment, "Wrong alignment of fourvector implementation");
  new (&impl_storage) _fourvector_impl_;
  impl = reinterpret_cast<_fourvector_impl_*>(&impl_storage);
}

GELATO::core::fourvector::fourvector(double x, double y, double z, double t) {
  static_assert(sizeof(_fourvector_impl_) <= c_impl_size, "Size of fourvector implementation too big");
  static_assert(alignof(_fourvector_impl_) == c_impl_alignment, "Wrong alignment of fourvector implementation");
  new (&impl_storage) _fourvector_impl_;
  impl = reinterpret_cast<_fourvector_impl_*>(&impl_storage);
  set_xyzt(x,y,z,t);
}

GELATO::core::fourvector::~fourvector() {
  impl->~_fourvector_impl_();
}

GELATO::core::fourvector::fourvector(const fourvector& v) {
  static_assert(sizeof(_fourvector_impl_) <= c_impl_size, "Size of fourvector implementation too big");
  static_assert(alignof(_fourvector_impl_) == c_impl_alignment, "Wrong alignment of fourvector implementation");
  new (&impl_storage) _fourvector_impl_(*(v.impl));
  impl = reinterpret_cast<_fourvector_impl_*>(&impl_storage);
}

GELATO::core::fourvector& GELATO::core::fourvector::operator=(const fourvector& v) {
  *impl = *v.impl;
  return *this;
}

/* No need for move semantics (all data stored locally in impl_storage)

// local storage, just copy
GELATO::core::fourvector::fourvector(fourvector&& v) {
  static_assert(sizeof(_fourvector_impl_) <= c_impl_size, "Size of fourvector implementation too big");
  static_assert(alignof(_fourvector_impl_) == c_impl_alignment, "Wrong alignment of fourvector implementation");
  new (&impl_storage) _fourvector_impl_(*(v.impl));
  impl = reinterpret_cast<_fourvector_impl_*>(&impl_storage);
}

GELATO::core::fourvector& GELATO::core::fourvector::operator=(fourvector&& v) {
  *impl = *v.impl;
  return *this;
}
*/

GELATO::core::fourvector& GELATO::core::fourvector::set_xyzt(double x, double y, double z, double t) {
#ifdef USING_CLHEP
  impl->set(x,y,z,t);
#elif defined USING_ROOT
  impl->SetXYZT(x,y,z,t);
#endif
  return *this;
}

GELATO::core::fourvector& GELATO::core::fourvector::set_t(double t) {
#ifdef USING_CLHEP
  impl->setT(t);
#elif defined USING_ROOT
  impl->SetT(t);
#endif
  return *this;
}

double GELATO::core::fourvector::x() const {
#ifdef USING_CLHEP
  return impl->x();
#elif defined USING_ROOT
  return impl->X();
#endif
}

double GELATO::core::fourvector::y() const {
#ifdef USING_CLHEP
  return impl->y();
#elif defined USING_ROOT
  return impl->Y();
#endif
}

double GELATO::core::fourvector::z() const {
#ifdef USING_CLHEP
  return impl->z();
#elif defined USING_ROOT
  return impl->Z();
#endif
}

double GELATO::core::fourvector::t() const {
#ifdef USING_CLHEP
  return impl->t();
#elif defined USING_ROOT
  return impl->T();
#endif
}

double GELATO::core::fourvector::px() const {
#ifdef USING_CLHEP
  return impl->px();
#elif defined USING_ROOT
  return impl->Px();
#endif
}

double GELATO::core::fourvector::py() const {
#ifdef USING_CLHEP
  return impl->py();
#elif defined USING_ROOT
  return impl->Py();
#endif
}

double GELATO::core::fourvector::pz() const {
#ifdef USING_CLHEP
  return impl->pz();
#elif defined USING_ROOT
  return impl->Pz();
#endif
}

double GELATO::core::fourvector::e() const {
#ifdef USING_CLHEP
  return impl->e();
#elif defined USING_ROOT
  return impl->E();
#endif
}

GELATO::core::fourvector& GELATO::core::fourvector::set_mass_momentum_theta_phi(double m, double mom, double t, double p) {
#ifdef USING_CLHEP
  impl->setVectM({mom*std::sin(t)*std::cos(p), mom*std::sin(t)*std::sin(p), mom*std::cos(t)},m);
#elif defined USING_ROOT
  impl->SetVectM({mom*std::sin(t)*std::cos(p), mom*std::sin(t)*std::sin(p), mom*std::cos(t)},m);
#endif
  return *this;
}

double GELATO::core::fourvector::m() const {
#ifdef USING_CLHEP
  return impl->m();
#elif defined USING_ROOT
  return impl->M();
#endif
}

double GELATO::core::fourvector::m2() const {
#ifdef USING_CLHEP
  return impl->m2();
#elif defined USING_ROOT
  return impl->M2();
#endif
}

GELATO::core::vector3 GELATO::core::fourvector::vect() const {
  vector3 ret;
#ifdef USING_CLHEP
  *ret.impl = impl->vect();
#elif defined USING_ROOT
  *ret.impl = impl->Vect();
#endif
  return ret;
}

GELATO::core::fourvector& GELATO::core::fourvector::set_vec_time(const vector3& vec, double time)  {
#ifdef USING_CLHEP
  impl->setVect(*vec.impl);
  impl->setT(time);
#elif defined USING_ROOT
  impl->SetVect(*vec.impl);
  impl->SetT(time);
#endif
  return *this;
}


GELATO::core::vector3 GELATO::core::fourvector::get_boost_vector() const {
  vector3 ret;
#ifdef USING_CLHEP
  *ret.impl = impl->boostVector();
#elif defined USING_ROOT
  *ret.impl = impl->BoostVector();
#endif
  return ret;
}

GELATO::core::fourvector& GELATO::core::fourvector::boost(const vector3& bv) {
#ifdef USING_CLHEP
  impl->boost(*bv.impl);
#elif defined USING_ROOT
  impl->Boost(*bv.impl);
#endif
  return *this;
}

double GELATO::core::fourvector::beta() const {
#ifdef USING_CLHEP
  return impl->beta();
#elif defined USING_ROOT
  return impl->Beta();
#endif
}

double GELATO::core::fourvector::gamma() const {
#ifdef USING_CLHEP
  return impl->gamma();
#elif defined USING_ROOT
  return impl->Gamma();
#endif
}


GELATO::core::fourvector GELATO::core::fourvector::operator-(const fourvector& v2) const {
  fourvector ret(*this);
  *ret.impl -= *v2.impl;
  return ret;
}

GELATO::core::fourvector GELATO::core::fourvector::operator+(const fourvector& v2) const {
  fourvector ret(*this);
  *ret.impl += *v2.impl;
  return ret;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







GELATO::core::rotation::rotation() {
  static_assert(sizeof(_rotation_impl_) <= c_impl_size, "Size of rotation implementation too big");
  static_assert(alignof(_rotation_impl_) == c_impl_alignment, "Wrong alignment of rotation implementation");
  new (&impl_storage) _rotation_impl_;
  impl = reinterpret_cast<_rotation_impl_*>(&impl_storage);
}

GELATO::core::rotation::~rotation() {
  impl->~_rotation_impl_();
}

GELATO::core::rotation::rotation(const GELATO::core::rotation& r) {
  static_assert(sizeof(_rotation_impl_) <= c_impl_size, "Size of rotation implementation too big");
  static_assert(alignof(_rotation_impl_) == c_impl_alignment, "Wrong alignment of rotation implementation");
  new (&impl_storage) _rotation_impl_(*(r.impl));
  impl = reinterpret_cast<_rotation_impl_*>(&impl_storage);
}

GELATO::core::rotation& GELATO::core::rotation::operator=(const rotation& r) {
  *impl = *r.impl;
  return *this;
}

/* No need for move semantics (all data stored locally in impl_storage)

// local storage, just copy
GELATO::core::rotation::rotation(GELATO::core::rotation&& r) {
  static_assert(sizeof(_rotation_impl_) <= c_impl_size, "Size of rotation implementation too big");
  static_assert(alignof(_rotation_impl_) == c_impl_alignment, "Wrong alignment of rotation implementation");
  new (&impl_storage) _rotation_impl_(*(r.impl));
  impl = reinterpret_cast<_rotation_impl_*>(&impl_storage);
}

// local storage, just copy
GELATO::core::rotation& GELATO::core::rotation::operator=(rotation&& r) {
  *impl = *r.impl;
  return *this;
}
*/

GELATO::core::vector3 GELATO::core::rotation::operator*(const vector3& vin) const {
  GELATO::core::vector3 ret;
  *ret.impl = *impl * *vin.impl;
  return ret;
}

GELATO::core::rotation& GELATO::core::rotation::invert() { 
#ifdef USING_CLHEP
  impl->invert();
#elif defined USING_ROOT
  impl->Invert();
#endif
  return *this;
}

GELATO::core::rotation GELATO::core::rotation::get_inverse()  const { 
  rotation ret;
#ifdef USING_CLHEP
  *ret.impl = impl->inverse();
#elif defined USING_ROOT
  *ret.impl = impl->Inverse();
#endif
  return ret;
}

GELATO::core::rotation& GELATO::core::rotation::rotate_axes(const vector3& new_x, const vector3& new_y, const vector3& new_z) {
#ifdef USING_CLHEP
  impl->rotateAxes(*new_x.impl, *new_y.impl, *new_z.impl);
#elif defined USING_ROOT
  impl->RotateAxes(*new_x.impl, *new_y.impl, *new_z.impl);
#endif
  return *this;
}



#undef USING_CLHEP
#undef USING_ROOT
#undef USING_XXXX


#endif // defined EXPOSE_PHYSICS_VECTORS
