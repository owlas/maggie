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
using array_d = boost::multi_array<double,1>;
using matrix_d = boost::multi_array<double,2>;
#include<boost/random.hpp>
using boost::mt19937;

class Simulation
{
public:
  // constructor
  Simulation( const ParticleCluster geom, const matrix_d init_state,
              const double h, const unsigned int n, const double temp,
              const array_d field );

  // Do simulation
  int run();

  // save results to hard drive
  int save();

  // approximate energy barriers for single particle flips
  matrix_d computeEnergyBarriers() const;

  // compute the transition matrix_d
  matrix_d arrheniusMatrix_D() const;

  // setters
  void setSimLength( const unsigned int );
  void setTimeStep( const double );
  void setState( const matrix_d );
  void setTemp( const double );
  void setField( const array_d );

  // getters
  ParticleCluster getGeometry() const;
  double getTimeStep() const;
  unsigned int getSimLength() const;
  matrix_d getState() const;
  double getTemp() const;
  array_d getField() const;

private:
  const ParticleCluster geom;
  double dt;
  unsigned int N;
  double T;
  array_d  h;
  matrix_d state;
  std::vector<double> stability;
};
#endif
