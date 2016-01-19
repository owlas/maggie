// SingleMNPMasterEquation.hpp - Class to represent the master
// equation for a single uniaxial magnetic nanoparticle. Derived from
// LangevunEquation base.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#ifndef MNPME_H
#define MNPME_H

#include<ODE.hpp>
#include<Constants.hpp>
#include<Field.hpp>
#include<KramersTrig.hpp>
#include<stdexcept>
#include<functional>
#include <boost/multi_array.hpp>
using array_f = boost::multi_array<float,1>;

class SingleMNPMasterEquation : public ODE<float>
{
public:
  // constructor for the general case
  SingleMNPMasterEquation( float anis,
			   float temp,
			   const Field &field,
			   float fieldangle,
			   float radius,
			   float alpha,
			   float Ms,
			   float gamma=Constants::GYROMAG );

  // delegated constructor for aligned case
  SingleMNPMasterEquation( float anis,
			   float temp,
			   const Field &field,
			   float radius );



  // compute the difference in energy between the two states of the
  // system
  float ediff( float t );

  // compute the angle of the minima
  float state_rotation( float t ) const;

  // Compute the drift of the master equation
  void computeDrift( array_f& out, const array_f& in, const float t ) const;

private:
  float k;
  float T;
  float r;
  float psi;
  float a;
  float M;
  float gam;
  float v;
  const Field *fieldPtr;
};

#endif
