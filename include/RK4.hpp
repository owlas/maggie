// RK4.hpp
// Header for the RK4 Integrator
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<Integrator.hpp>

class RK4 : public Integrator
{

public:
  //Constructor
  RK4( LangevinEquation &le, array_f& init_state, float init_time,
       float dt );

  // Step the integrator once
  virtual void step();

private:
  float h;
  float dim;
  array_f k1;
  array_f k2;
  array_f k3;
  array_f k4;
  array_f tmp;
};
