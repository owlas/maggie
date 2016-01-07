// Integrator.hpp
// Header for the absract base class for integrators
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include<SDE.hpp>
#include<boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;
#include<iostream>
using std::cout;
using std::endl;
#include<Integrator.hpp>
#include<stdexcept>
using std::invalid_argument;

template<class T>
class Integrator
{
 public:

  Integrator( const T& equation, const array_f&,
	      const float time=0 );

  // Get the current state
  array_f getState() const;
  // Set the integrator state
  void setState( const array_f );

  // Get the current time
  float getTime() const;
  // Set the integrator time
  void setTime( const float );

  // Get the pointer to the associated Langevin Equation
  const T& getLE() const;

  // Compute an integration step
  virtual void step() = 0;

  // Reset the integrator to initial conditions
  // or a new initial condition
  void reset();
  void reset( const array_f );

protected:
  array_f state;

 private:  
  const T* lEq;
  array_f initial_state;
  float t;
  float initial_t;
};
  
#include "tpp/Integrator.tpp"
#endif
