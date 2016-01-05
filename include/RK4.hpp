// RK4.hpp
// Header for the RK4 Integrator
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef RK4_H
#define RK4_H
#include<Integrator.hpp>

class RK4 : public Integrator
{

public:
  //Constructor
  RK4( const LangevinEquation &le, const array_f& init_state,
       const float init_time, const float dt );

  // Step the integrator once
  virtual void step();

private:
  const float h;
  const float dim;
  array_f k1;
  array_f k2;
  array_f k3;
  array_f k4;
  array_f tmp;
};

#endif
