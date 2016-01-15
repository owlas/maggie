// Milstein.hpp
// Header for the Milstein numerical solver for SDEs
// 
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
// 
#ifndef MILST_H
#define MILST_H

#include <Integrator.hpp>
#include <boost/multi_array.hpp>
using array_f = boost::multi_array<float,1>;
using matrix_f = boost::multi_array<float,2>;
using array3_f = boost::multi_array<float,3>;

#include <boost/random.hpp>
using boost::variate_generator;
using boost::mt19937;
using boost::normal_distribution;
#include <StocLLG.hpp>
#include <cmath>
using std::sqrt; using std::pow;
#include <algorithm>

class Milstein : public Integrator<SDE>
{
public:
  // Constructor
  Milstein( const SDE &le, const array_f& init_state,
	    const float time, const float dt,
	    mt19937 &rng_1, mt19937 &rng_2 ); 

  // One step of the integrator
  virtual void step();

  // Fourier series approximation of double integral
  // may only be applicable for stratonovich
  void setP( const int pset );
  int getP() const;

  // Manual wiener process mode
  // allow user to control internal wiener increments themselves
  void setManualWienerMode( const bool );
  void setWienerIncrements( const array_f );

private:
  const float h;
  const int dim, wDim;
  normal_distribution<float> dist_1;
  normal_distribution<float> dist_2;
  variate_generator<mt19937&, normal_distribution<float> > gen_1;
  variate_generator<mt19937&, normal_distribution<float> > gen_2;
  array_f next_state;
  array_f dw;
  array_f dw2;
  array_f tmp1;
  matrix_f tmp2;
  matrix_f tmp22;
  array3_f tmp3;
  unsigned int p; // Fourier series truncation
  array_f mu;
  matrix_f eta;
  matrix_f zeta;
  bool manualWiener{ false };
};
#endif
