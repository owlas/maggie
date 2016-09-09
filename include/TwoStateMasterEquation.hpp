// TwoStateMasterEquation.hpp
// Implementation of the differential equation for a two state master
// equation.
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef TWOSTATEME_H
#define TWOSTATEME_H

#include<ODE.hpp>
#include<functional>
#include<stdexcept>
#include<iostream>
using std::cout;
using std::endl;
using std::invalid_argument;

#include <boost/multi_array.hpp>
using array2 = std::array<double,2>;

class TwoStateMasterEquation : public ODE<2>
{
public:

  // constructor for fixed rates
  TwoStateMasterEquation( const double rate1, const double rate2 );

  // constructor for variable rates
  TwoStateMasterEquation( const std::function<double( double )> rateFunc1,
                          const std::function<double( double)> rateFunc2 );

  // Compute the drift vector for the current state
    void computeDrift( ODE<2>::array& out, const ODE<2>::array& in, const double t )
      const;

  // get rate one and two at a set time
  double getRate1( const double t=0) const;
  double getRate2( const double t=0) const;

private:
  const std::function<double( double )> w1;
  const std::function<double( double )> w2;
};

#endif
