// Field.cpp - implementation of time-varying external field shapes
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<Field.hpp>

float Field::sin( float h, float f, float t )
{
  return h*std::sin( 2*M_PI*f*t );
}
float Field::square( float h, float f, float t )
{
  return h*( int( t*f*2 )%2 ? -1 : 1 );
}
float Field::cos( float h, float f, float t )
{
  return h*std::cos( 2*M_PI*f*t );
}
