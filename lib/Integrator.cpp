// Integrator.cpp
// Implementation for the abstract Integrator class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<iostream>
using std::cout;
using std::endl;

#include<LangevinEquation.hpp>
#include<boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;

#include<Integrator.hpp>
#include<stdexcept>
using std::invalid_argument;

// Constructor
Integrator::Integrator( const LangevinEquation& le, const array_f& init_state,
			const float time )
  : state( boost::extents[le.getDim()] )
  , lEq( &le )
  , initial_state( boost::extents[le.getDim()] )
{
  initial_t = time;
  setTime( time );
  initial_state = init_state;
  setState( init_state );
}

// Get state
array_f Integrator::getState() const { return state; }
// Set State
void Integrator::setState( const array_f s )
{
if( int( s.shape()[0] ) != lEq->getDim() )
    {
      throw invalid_argument( "Error: state dimension does not match"
			      " Langevin equation dimension");
    }
  else
    {
      state = s;
    }
}

// get time
float Integrator::getTime() const { return t; }
// set time
void Integrator::setTime( const float time ) { t = time; }

// get the pointer to the integrator
const LangevinEquation &Integrator::getLE() const { return *lEq; }

// reset the integrator to the initial condition
void Integrator::reset()
{
  state = initial_state;
  t = initial_t;
}

// reset the integrator with a new inital condition
void Integrator::reset( const array_f new_init )
{
  initial_state = new_init;
  reset();
}
