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
using array_f = boost::multi_array<float,1>;

class TwoStateMasterEquation : public ODE<float>
{
public:

  // constructor for fixed rates
  TwoStateMasterEquation( const float rate1, const float rate2 );

  // constructor for variable rates
  TwoStateMasterEquation( const std::function<float( float )> rateFunc1,
                          const std::function<float( float )> rateFunc2 );

  // Compute the drift vecctor for the current state
  virtual void computeDrift( array_f& out, const array_f& in, const float t )
    const;

  // get rate one and two at a set time
  float getRate1( const float t=0) const;
  float getRate2( const float t=0) const;

private:
  const std::function<float( float )> w1;
  const std::function<float( float )> w2;
};

#endif
