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
#include "MagneticStateController.hpp"
#include <Heun.hpp>
#include "Euler.hpp"
#include <StocLLG.hpp>
#include <Utility.hpp>
#include<stdexcept>
#include<boost/multi_array.hpp>
using array_d = boost::multi_array<double,1>;
using boost::extents;
using bidx = boost::multi_array_types::index;

#include "types.hpp"
using namespace maggie;

#include <omp.h>

#include<memory>

#include<vector>
using ad_vec = std::vector<moment>;

#include<boost/random.hpp>
using boost::mt19937;

#include<cmath>
using std::sin; using std::pow;
using std::sqrt; using std::exp;
using std::acos;

//#include<omp.h>

#include<sstream>

class Simulation
{
public:
    // constructor
    Simulation( const ParticleCluster geom, const std::vector<moment>,
                const double stepsize, const unsigned int n, const temperature t,
                const field happ , const unsigned long int seed=5091);

    // Initialise the simulation
    int init();

  // Do simulation for a single system and save the mag data
    int runFull();

    // Do a simulation of an ensemble of particles
    int runEnsemble( unsigned int );

    // Do a simulation of an ensemble
    int runFullEnsemble( unsigned int );

    // Do a simulation of a single particle until it switches
    int runFPT( const int N_ensemble, const bool alignup = false );

  // save results to hard drive
  int save();

  // approximate energy barriers for single particle flips
  matrix_d computeEnergyBarriers() const;

  // compute the transition matrix_d
  matrix_d arrheniusMatrix_D() const;

    // computes a random state from the equilibrium distribution
  std::vector<maggie::moment> equilibriumState();

    // compute the switch times
    int runResidence( const unsigned int N_switches );

  // setters
  void setSimLength( const unsigned int );
  void setTimeStep( const double );
  void setState( const std::vector<moment> );
  void setTemp( const temperature );
  void setField( const field );

  // getters
  ParticleCluster getGeometry() const;
  double getTimeStep() const;
  unsigned int getSimLength() const;
  const std::vector<moment>& getState() const;
  temperature getTemp() const;
  field getField() const;

private:
    const ParticleCluster geom;
    double dt;
    unsigned int Nsteps;
    temperature T;
    field h;
    const size_t Nparticles;
    std::vector<moment> state;
    std::vector<stability> sigmas;
    mt19937 equilibrium_rng;
    MagneticStateController<Euler<StocLLG>>::unique_ptr simulationController;
    const unsigned long int integrators_seed;
    double reduced_time_factor;
};
#endif
