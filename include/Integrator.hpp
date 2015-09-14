// Integrator.hpp
// Header for the absract base class for integrators
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include<LangevinEquation.hpp>
#include<boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;

class Integrator
{
 public:

  Integrator( LangevinEquation &ld, array_f&, float time=0 );

  // Get the current state
  array_f& getState();
  // Set the integrator state
  void setState( array_f );

  // Get the current time
  float getTime();
  // Set the integrator time
  void setTime( float );

  // Compute an integration step
  virtual void step() const = 0;

  // Reset the integrator to initial conditions
  void reset();

 private:  
  array_f state;
  array_f initial_state;
  float t;
  float initial_t;
  LangevinEquation *lEq;
};  
  

#endif

