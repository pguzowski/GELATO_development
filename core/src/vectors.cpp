#ifdef EXPOSE_PHYSICS_VECTORS

#else // EXPOSE_PHYSICS_VECTORS
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
    struct vector3::_vector3_impl_ : public VECTOR3_CLASS {
      public:
        using VECTOR3_CLASS::operator=;
    };
    //  VECTOR3_CLASS v;
    //};
    struct fourvector::_fourvector_impl_ : public VECTOR4_CLASS {
      public:
        //using VECTOR4_CLASS::operator=;
    };
    //  VECTOR4_CLASS v;
    //};
    struct rotation::_rotation_impl_ : public ROTATION_CLASS {
      public:
        using ROTATION_CLASS::operator=;
    };
    //  ROTATION_CLASS r;
    //};
  }
}


dkgen::core::vector3::vector3() : impl{
  //std::make_unique<_vector3_impl_>()
  //std::unique_ptr<_vector3_impl_>(new (&impl_storage) _vector3_impl_ )
  new (&impl_storage) _vector3_impl_ 
} {
  static_assert(sizeof(_vector3_impl_) <= c_impl_size, "Size of vector3 implementation too big");
  static_assert(alignof(_vector3_impl_) == c_impl_alignment, "Wrong alignment of vector3 implementation");
}

dkgen::core::vector3::vector3(double x, double y, double z) : impl{
  //std::make_unique<_vector3_impl_>()
  //std::unique_ptr<_vector3_impl_>(new (&impl_storage) _vector3_impl_ )
  new (&impl_storage) _vector3_impl_ 
} {
  //static_assert(sizeof(_vector3_impl_) <= sizeof(impl_storage), "Size of vector3 implementation too big");
  static_assert(sizeof(_vector3_impl_) <= c_impl_size, "Size of vector3 implementation too big");
  static_assert(alignof(_vector3_impl_) == c_impl_alignment, "Wrong alignment of vector3 implementation");
  set_xyz(x,y,z);
}

dkgen::core::vector3::~vector3() {
  // no need to delete local storage, and VECTOR3_CLASS doesn't leak memory in the implementations that have been tried
  //delete impl.get();
}

dkgen::core::vector3::vector3(const vector3& v) : impl{
  //std::make_unique<_vector3_impl_>(*(v.impl))
  //std::unique_ptr<_vector3_impl_>(new (&impl_storage) _vector3_impl_(*(v.impl)) )
  new (&impl_storage) _vector3_impl_(*(v.impl)) 
} {
  static_assert(sizeof(_vector3_impl_) <= c_impl_size, "Size of vector3 implementation too big");
  static_assert(alignof(_vector3_impl_) == c_impl_alignment, "Wrong alignment of vector3 implementation");
}

dkgen::core::vector3::vector3(vector3&& v) : impl{
  //std::move(v.impl)
  //std::make_unique<_vector3_impl_>(*(v.impl))
  //std::unique_ptr<_vector3_impl_>(new (&impl_storage) _vector3_impl_(*(v.impl)) )
  new (&impl_storage) _vector3_impl_(*(v.impl)) 
} {
  static_assert(sizeof(_vector3_impl_) <= c_impl_size, "Size of vector3 implementation too big");
  static_assert(alignof(_vector3_impl_) == c_impl_alignment, "Wrong alignment of vector3 implementation");
}

dkgen::core::vector3& dkgen::core::vector3::operator=(const dkgen::core::vector3& v) {
  *impl = *v.impl;
  return *this;
}

dkgen::core::vector3& dkgen::core::vector3::operator=(dkgen::core::vector3&& v) {
  //impl = std::move(v.impl);
  *impl = *v.impl;
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
  return impl->mag();
#elif defined(USING_ROOT)
  return impl->Mag();
#endif
}

dkgen::core::vector3& dkgen::core::vector3::set_mag(double m) {
#ifdef USING_CLHEP
  impl->setMag(m);
#elif defined(USING_ROOT)
  impl->SetMag(m);
#endif
  return *this;
}

dkgen::core::vector3& dkgen::core::vector3::set_xyz(double x, double y, double z) {
#ifdef USING_CLHEP
  impl->set(x,y,z);
#elif defined(USING_ROOT)
  impl->SetXYZ(x,y,z);
#endif
  return *this;
}

double dkgen::core::vector3::x() const {
#ifdef USING_CLHEP
  return impl->x();
#elif defined USING_ROOT
  return impl->X();
#endif
}

double dkgen::core::vector3::y() const {
#ifdef USING_CLHEP
  return impl->y();
#elif defined USING_ROOT
  return impl->Y();
#endif
}

double dkgen::core::vector3::z() const {
#ifdef USING_CLHEP
  return impl->z();
#elif defined USING_ROOT
  return impl->Z();
#endif
}

dkgen::core::vector3 dkgen::core::vector3::operator-(const dkgen::core::vector3& v2) const {
  vector3 ret(*this);
  *ret.impl -= *v2.impl;
  return ret;
}

dkgen::core::vector3 dkgen::core::vector3::operator+(const dkgen::core::vector3& v2) const {
  vector3 ret(*this);
  *ret.impl += *v2.impl;
  return ret;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////





dkgen::core::fourvector::fourvector() : impl{
  //std::make_unique<_fourvector_impl_>()
  new (&impl_storage) _fourvector_impl_
} {
  static_assert(sizeof(_fourvector_impl_) <= c_impl_size, "Size of fourvector implementation too big");
  static_assert(alignof(_fourvector_impl_) == c_impl_alignment, "Wrong alignment of fourvector implementation");
}

dkgen::core::fourvector::fourvector(double x, double y, double z, double t) : impl{
  //std::make_unique<_fourvector_impl_>()
  new (&impl_storage) _fourvector_impl_
} {
  static_assert(sizeof(_fourvector_impl_) <= c_impl_size, "Size of fourvector implementation too big");
  static_assert(alignof(_fourvector_impl_) == c_impl_alignment, "Wrong alignment of fourvector implementation");
  set_xyzt(x,y,z,t);
}

dkgen::core::fourvector::~fourvector() { }

dkgen::core::fourvector::fourvector(const dkgen::core::fourvector& v) : impl{
  //std::make_unique<_fourvector_impl_>(*(v.impl))
  new (&impl_storage) _fourvector_impl_(*(v.impl))
} {
  static_assert(sizeof(_fourvector_impl_) <= c_impl_size, "Size of fourvector implementation too big");
  static_assert(alignof(_fourvector_impl_) == c_impl_alignment, "Wrong alignment of fourvector implementation");
}

dkgen::core::fourvector::fourvector(dkgen::core::fourvector&& v) : impl{
  //std::move(v.impl)
  new (&impl_storage) _fourvector_impl_(*(v.impl))
} {
  static_assert(sizeof(_fourvector_impl_) <= c_impl_size, "Size of fourvector implementation too big");
  static_assert(alignof(_fourvector_impl_) == c_impl_alignment, "Wrong alignment of fourvector implementation");
}

dkgen::core::fourvector& dkgen::core::fourvector::operator=(const dkgen::core::fourvector& v) {
  *impl = *v.impl;
  return *this;
}

dkgen::core::fourvector& dkgen::core::fourvector::operator=(dkgen::core::fourvector&& v) {
  //impl = std::move(v.impl);
  *impl = *v.impl;
  return *this;
}

dkgen::core::fourvector& dkgen::core::fourvector::set_xyzt(double x, double y, double z, double t) {
#ifdef USING_CLHEP
  impl->set(x,y,z,t);
#elif defined USING_ROOT
  impl->SetXYZT(x,y,z,t);
#endif
  return *this;
}

dkgen::core::fourvector& dkgen::core::fourvector::set_t(double t) {
#ifdef USING_CLHEP
  impl->setT(t);
#elif defined USING_ROOT
  impl->SetT(t);
#endif
  return *this;
}

double dkgen::core::fourvector::x() const {
#ifdef USING_CLHEP
  return impl->x();
#elif defined USING_ROOT
  return impl->X();
#endif
}

double dkgen::core::fourvector::y() const {
#ifdef USING_CLHEP
  return impl->y();
#elif defined USING_ROOT
  return impl->Y();
#endif
}

double dkgen::core::fourvector::z() const {
#ifdef USING_CLHEP
  return impl->z();
#elif defined USING_ROOT
  return impl->Z();
#endif
}

double dkgen::core::fourvector::t() const {
#ifdef USING_CLHEP
  return impl->t();
#elif defined USING_ROOT
  return impl->T();
#endif
}

double dkgen::core::fourvector::px() const {
#ifdef USING_CLHEP
  return impl->px();
#elif defined USING_ROOT
  return impl->Px();
#endif
}

double dkgen::core::fourvector::py() const {
#ifdef USING_CLHEP
  return impl->py();
#elif defined USING_ROOT
  return impl->Py();
#endif
}

double dkgen::core::fourvector::pz() const {
#ifdef USING_CLHEP
  return impl->pz();
#elif defined USING_ROOT
  return impl->Pz();
#endif
}

double dkgen::core::fourvector::e() const {
#ifdef USING_CLHEP
  return impl->e();
#elif defined USING_ROOT
  return impl->E();
#endif
}

dkgen::core::fourvector& dkgen::core::fourvector::set_mass_momentum_theta_phi(double m, double mom, double t, double p) {
#ifdef USING_CLHEP
  impl->setVectM({mom*std::sin(t)*std::cos(p), mom*std::sin(t)*std::sin(p), mom*std::cos(t)},m);
#elif defined USING_ROOT
  impl->SetVectM({mom*std::sin(t)*std::cos(p), mom*std::sin(t)*std::sin(p), mom*std::cos(t)},m);
#endif
  return *this;
}

double dkgen::core::fourvector::m() const {
#ifdef USING_CLHEP
  return impl->m();
#elif defined USING_ROOT
  return impl->M();
#endif
}

double dkgen::core::fourvector::m2() const {
#ifdef USING_CLHEP
  return impl->m2();
#elif defined USING_ROOT
  return impl->M2();
#endif
}

dkgen::core::vector3 dkgen::core::fourvector::vect() const {
  vector3 ret;
#ifdef USING_CLHEP
  *ret.impl = impl->vect();
#elif defined USING_ROOT
  *ret.impl = impl->Vect();
#endif
  return ret;
}

dkgen::core::fourvector& dkgen::core::fourvector::set_vec_time(const dkgen::core::vector3& vec, double time)  {
#ifdef USING_CLHEP
  impl->setVect(*vec.impl);
  impl->setT(time);
#elif defined USING_ROOT
  impl->SetVect(*vec.impl);
  impl->SetT(time);
#endif
  return *this;
}


dkgen::core::vector3 dkgen::core::fourvector::get_boost_vector() const {
  vector3 ret;
#ifdef USING_CLHEP
  *ret.impl = impl->boostVector();
#elif defined USING_ROOT
  *ret.impl = impl->BoostVector();
#endif
  return ret;
}

dkgen::core::fourvector& dkgen::core::fourvector::boost(const dkgen::core::vector3& bv) {
#ifdef USING_CLHEP
  impl->boost(*bv.impl);
#elif defined USING_ROOT
  impl->Boost(*bv.impl);
#endif
  return *this;
}

double dkgen::core::fourvector::beta() const {
#ifdef USING_CLHEP
  return impl->beta();
#elif defined USING_ROOT
  return impl->Beta();
#endif
}


dkgen::core::fourvector dkgen::core::fourvector::operator-(const dkgen::core::fourvector& v2) const {
  fourvector ret(*this);
  *ret.impl -= *v2.impl;
  return ret;
}

dkgen::core::fourvector dkgen::core::fourvector::operator+(const dkgen::core::fourvector& v2) const {
  fourvector ret(*this);
  *ret.impl += *v2.impl;
  return ret;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







dkgen::core::rotation::rotation() : impl{
  //std::make_unique<_rotation_impl_>()
  new (&impl_storage) _rotation_impl_
} {
  static_assert(sizeof(_rotation_impl_) <= c_impl_size, "Size of rotation implementation too big");
  static_assert(alignof(_rotation_impl_) == c_impl_alignment, "Wrong alignment of rotation implementation");
}

dkgen::core::rotation::~rotation() { }

dkgen::core::rotation::rotation(const dkgen::core::rotation& r) : impl{
  //std::make_unique<_rotation_impl_>(*(r.impl))
  new (&impl_storage) _rotation_impl_(*(r.impl))
} {
  static_assert(sizeof(_rotation_impl_) <= c_impl_size, "Size of rotation implementation too big");
  static_assert(alignof(_rotation_impl_) == c_impl_alignment, "Wrong alignment of rotation implementation");
}

dkgen::core::rotation::rotation(dkgen::core::rotation&& r) : impl{
  //std::move(r.impl)
  new (&impl_storage) _rotation_impl_(*(r.impl))
} {
  static_assert(sizeof(_rotation_impl_) <= c_impl_size, "Size of rotation implementation too big");
  static_assert(alignof(_rotation_impl_) == c_impl_alignment, "Wrong alignment of rotation implementation");
}

dkgen::core::rotation& dkgen::core::rotation::operator=(const dkgen::core::rotation& r) {
  *impl = *r.impl;
  return *this;
}

dkgen::core::rotation& dkgen::core::rotation::operator=(dkgen::core::rotation&& r) {
  //impl = std::move(r.impl);
  *impl = *r.impl;
  return *this;
}

dkgen::core::vector3 dkgen::core::rotation::operator*(const dkgen::core::vector3& vin) const {
  dkgen::core::vector3 ret;
  *ret.impl = *impl * *vin.impl;
  return ret;
}

dkgen::core::rotation& dkgen::core::rotation::invert() { 
#ifdef USING_CLHEP
  impl->invert();
#elif defined USING_ROOT
  impl->Invert();
#endif
  return *this;
}

dkgen::core::rotation dkgen::core::rotation::get_inverse()  const { 
  rotation ret;
#ifdef USING_CLHEP
  *ret.impl = impl->inverse();
#elif defined USING_ROOT
  *ret.impl = impl->Inverse();
#endif
  return ret;
}

dkgen::core::rotation& dkgen::core::rotation::rotate_axes(const dkgen::core::vector3& new_x, const dkgen::core::vector3& new_y, const dkgen::core::vector3& new_z) {
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
