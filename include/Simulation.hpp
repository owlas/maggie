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
#include<boost/multi_array.hpp>
using array_d = boost::multi_array<double,1>;
using boost::extents;
using bidx = boost::multi_array_types::index;

#include<vector>
using ad_vec = std::vector<array_d>;

#include<boost/random.hpp>
using boost::mt19937;

#include<cmath>
using std::sin; using std::pow;
using std::sqrt; using std::exp;
using std::acos;

class Simulation
{
public:
  // constructor
  Simulation( const ParticleCluster geom, const ad_vec init_state,
              const double h, const unsigned int n, const double temp,
              const array_d field );

  // Do simulation for a single system and save the mag data
    int runFull();

    // Do a simulation of an ensemble of particles
    int runEnsemble( unsigned int );

  // save results to hard drive
  int save();

  // approximate energy barriers for single particle flips
  matrix_d computeEnergyBarriers() const;

  // compute the transition matrix_d
  matrix_d arrheniusMatrix_D() const;

    // computes a random state from the equilibrium distribution
    array_d equilibriumState();

  // setters
  void setSimLength( const unsigned int );
  void setTimeStep( const double );
  void setState( const ad_vec );
  void setTemp( const double );
  void setField( const array_d );

  // getters
  ParticleCluster getGeometry() const;
  double getTimeStep() const;
  unsigned int getSimLength() const;
  const ad_vec& getState() const;
  double getTemp() const;
  array_d getField() const;

private:
  const ParticleCluster geom;
  double dt;
  unsigned int N;
  double T;
  array_d  h;
  ad_vec state;
  std::vector<double> stability;
    mt19937 equilibrium_rng;
};
#endif
