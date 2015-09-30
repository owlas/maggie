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
using array_f = boost::multi_array<float, 1>;

#include<stdexcept>
using std::invalid_argument;

class Particle
{
public:
  // constructor
  Particle( float gamma, float alpha, float Ms, float D, float K,
	    array_f uea );

  // getters for particle properties
  float getGamma() const;
  float getAlpha() const;
  float getMs() const;
  float getD() const;
  float getK() const;
  array_f getUea() const;
  float getV() const;

  // setters for particle properties
  void setGamma( float );
  void setAlpha( float );
  void setMs( float );
  void setSize( float );
  void setK( float );
  void setUea( array_f );

  // Compute the energy barriers for the case of an applied field
  // parallel to the uniaxial anisotropy axis and with an intensity of
  // less that H_k. i.e. h = H/H_k < 1
  array_f alignedEnergyBarriers( float happ ) const;

  // Compute the Neel-Arrhenius transition rates at a given temperature
  array_f neelTransitionRates( float temp, float barrier1,
                               float barrier2 )     const;
  

private:
  array_f uea;
  float gamma;
  float alpha;
  float ms;
  float d;
  float k;
  float v;
};
#endif
