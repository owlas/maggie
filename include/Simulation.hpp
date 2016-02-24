// Simulation.hpp
//
// Header file for a magnetic nanoparticle simulation. Current
// behaviour allows the user to run a master equation integration to
// determine the relaxation of a magnetic particle cluster.
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//
#ifndef SIM_H
#define SIM_H

#include<ParticleCluster.hpp>
#include <Heun.hpp>
#include <StocLLG.hpp>
#include <Utility.hpp>
#include<stdexcept>
#include<vector>
#include<boost/multi_array.hpp>
using array_f = boost::multi_array<float,1>;
using matrix_f = boost::multi_array<float,2>;
#include<boost/random.hpp>
using boost::mt19937;

class Simulation
{
public:
  // constructor
  Simulation( const ParticleCluster geom, const matrix_f init_state,
              const float h, const unsigned int n, const float temp,
              const array_f field );

  // Do simulation
  int run();

  // save results to hard drive
  int save();

  // approximate energy barriers for single particle flips
  matrix_f computeEnergyBarriers() const;

  // compute the transition matrix
  matrix_f arrheniusMatrix() const;

  // setters
  void setSimLength( const unsigned int );
  void setTimeStep( const float );
  void setState( const matrix_f );
  void setTemp( const float );
  void setField( const array_f );

  // getters
  ParticleCluster getGeometry() const;
  float getTimeStep() const;
  unsigned int getSimLength() const;
  matrix_f getState() const;
  float getTemp() const;
  array_f getField() const;

private:
  const ParticleCluster geom;
  float dt;
  unsigned int N;
  float T;
  array_f h;
  matrix_f state;
  std::vector<float> stability;
};
#endif
