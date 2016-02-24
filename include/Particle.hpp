// Particle.hpp
// Head file for the particle class. Represents a single magnetic
// nanoparticle with its parameters.
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//
#ifndef PARTICLE_H
#define PARTICLE_H

#define _USE_MATH_DEFINES
#include<cmath>
using std::sqrt;
using std::pow;
using std::fabs;
using std::exp;

#include<Constants.hpp>

#include<boost/multi_array.hpp>
using array_d = boost::multi_array<double, 1>;

#include<stdexcept>
using std::invalid_argument;

class Particle
{
public:
  // constructor
  Particle( double gamma, double alpha, double Ms, double D, double K,
	    array_d uea );

  // getters for particle properties
  double getGamma() const;
  double getAlpha() const;
  double getMs() const;
  double getD() const;
  double getK() const;
  array_d getUea() const;
  double getV() const;

  // setters for particle properties
  void setGamma( double );
  void setAlpha( double );
  void setMs( double );
  void setSize( double );
  void setK( double );
  void setUea( array_d );

  // Compute the energy barriers for the case of an applied field
  // parallel to the uniaxial anisotropy axis and with an intensity of
  // less that H_k. i.e. h = H/H_k < 1
  array_d alignedEnergyBarriers( double happ ) const;

  // Compute the Neel-Arrhenius transition rates at a given temperature
  array_d neelTransitionRates( double temp, double barrier1,
                               double barrier2 )     const;


private:
  array_d uea;
  double gamma;
  double alpha;
  double ms;
  double d;
  double k;
  double v;
};
#endif
