// Heun.hpp
// Header for the Heun integration scheme
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef HEUN_H
#define HEUN_H

#include<Integrator.hpp>
#include<cmath>
using std::sqrt;
#include<boost/random.hpp>
using boost::variate_generator;
using boost::mt19937;
using boost::normal_distribution;
#include<boost/multi_array.hpp>
using array_f = boost::multi_array<float,1>;
using matrix_f = boost::multi_array<float,2>;

class Heun : public Integrator<SDE>
{

public:
  // Constructor
  Heun( const SDE &sde, const array_f& init_state,
	const float time, const float dt, mt19937 &rng );

  // Step the integrator once
  virtual void step();

  // Manual wiener process mode
  // allow user to control internal wiener increments themselves
  void setManualWienerMode( const bool );
  void setWienerIncrements( const array_f );

private:
  const float h;
    const int dim, wDim;
  array_f dw;
  array_f xPred;
  array_f tmp1;
  array_f tmp1Up;
  matrix_f tmp2;
  matrix_f tmp2Up;
  normal_distribution<float> dist;
  variate_generator<mt19937 &, normal_distribution<float> > gen;
  bool manualWiener{ false };
};
#endif
