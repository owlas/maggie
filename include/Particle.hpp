// Particle.hpp
// Head file for the particle class. Represents a single magnetic
// nanoparticle with its parameters.
// 
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
// 
#ifndef PARTICLE_H
#define PARTICLE_H

#include<boost/multi_array.hpp>
typedef boost::multi_array<float, 1> array_f;

class Particle
{
public:
  // constructor
  Particle( float gamma, float alpha, float Ms, float D, float K,
	    array_f uea );

  // getters for particle properties
  float getGamma();
  float getAlpha();
  float getMs();
  float getD();
  float getK();
  array_f getUea();
  float getV();

  // setters for particle properties
  void setGamma( float );
  void setAlpha( float );
  void setMs( float );
  void setSize( float );
  void setK( float );
  void setUea( array_f );

  // Compute the Neel-Brown transition rates over the barrier
  

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
