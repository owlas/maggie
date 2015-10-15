// Field.cpp - implementation of time-varying external field shapes
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<Field.hpp>

// Abstract class to hold a field calculation
Field::Field( float field )
{
  setH( field );
}
void Field::setH( float field ) { h=field; }
float Field::getH() const { return h; }

// Abstrsct class for periodic time-varying fields
FieldPeriodic::FieldPeriodic( float field, float freq )
        : Field( field )
{
  setF( freq );
}
void FieldPeriodic::setF( float freq ) { f=freq; }
float FieldPeriodic::getF() const { return f; }

// Sine field
FieldACSine::FieldACSine( float field, float freq )
  : FieldPeriodic( field, freq ) {}
float FieldACSine::getField( float t ) const
{
  return getH()*std::sin( 2*M_PI*getF()*t );
}

// Cosine field
FieldACCosine::FieldACCosine( float field, float freq )
  : FieldPeriodic( field, freq ) {}
float FieldACCosine::getField( float t ) const
{
  return getH()*std::cos( 2*M_PI*getF()*t );
}

// Square wave
FieldACSquare::FieldACSquare( float field, float freq )
  : FieldPeriodic( field, freq ) {}
float FieldACSquare::getField( float t ) const
{
  return getH()*( int( t*getF()*2 )%2 ? -1 : 1 );
}

// Ramped square wave
FieldACRamp::FieldACRamp( float field, float freq,
			  float riseTime )
  : FieldPeriodic( field, freq )
{
  setRT( riseTime );
}
float FieldACRamp::getRT() const { return rt; }
void FieldACRamp::setRT( float r )
{
  // rise time must be 0-50%
  if( ( r>0.5 ) or ( r<0.0 ) )
    throw std::invalid_argument("The rise time must be between"
				" 0.0 - 0.5 (proportion of"
				" cyclical period)" );
  // set the rise time
  rt = r;
  // compute the ramping gradient
  float T=1/getF();
  ramp = 2*getH()/( T*rt );
  // compute the max interval
  left=T/2*rt;
  right=T/2*( 1-rt );
}
  
float FieldACRamp::getField( float t ) const
{
  float tcyc = std::fmod( t, 1/( getF() ) );
  float T2 = 1/( 2*getF() );
  return tcyc<left ? tcyc*ramp :
    tcyc<right ? getH() :
    tcyc<( T2+left ) ? getH() - ramp*( tcyc-right ) :
    tcyc<( T2+right ) ? -getH() :
    -getH() + ramp*( tcyc-T2-right );
}
