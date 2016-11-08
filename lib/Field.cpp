// Field.cpp - implementation of time-varying external field shapes
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<Field.hpp>

// Abstract class to hold a field calculation
Field::Field( double field )
{
  setH( field );
}
void Field::setH( double field ) { h=field; }
double Field::getH() const { return h; }

// The most basic field is simply a constant
FieldConstant::FieldConstant( double field )
  : Field( field ) {}
double FieldConstant::getField( double ) const
{ return getH(); }

// Abstrsct class for periodic time-varying fields
FieldPeriodic::FieldPeriodic( double field, double freq )
        : Field( field )
{
  setF( freq );
}
void FieldPeriodic::setF( double freq ) { f=freq; }
double FieldPeriodic::getF() const { return f; }

// Sine field
FieldACSine::FieldACSine( double field, double freq )
  : FieldPeriodic( field, freq ) {}
double FieldACSine::getField( double t ) const
{
  return getH()*std::sin( 2*M_PI*getF()*t );
}

// Cosine field
FieldACCosine::FieldACCosine( double field, double freq )
  : FieldPeriodic( field, freq ) {}
double FieldACCosine::getField( double t ) const
{
  return getH()*std::cos( 2*M_PI*getF()*t );
}

// Square wave
FieldACSquare::FieldACSquare( double field, double freq )
  : FieldPeriodic( field, freq ) {}
double FieldACSquare::getField( double t ) const
{
  return getH()*( int( t*getF()*2 )%2 ? -1 : 1 );
}

// Ramped square wave
FieldACRamp::FieldACRamp( double field, double freq,
			  double riseTime )
  : FieldPeriodic( field, freq )
{
  setRT( riseTime );
}
double FieldACRamp::getRT() const { return rt; }
void FieldACRamp::setRT( double r )
{
  // rise time must be 0-50%
  if( ( r>0.5 ) or ( r<0.0 ) )
    throw std::invalid_argument("The rise time must be between"
				" 0.0 - 0.5 (proportion of"
				" cyclical period)" );
  // set the rise time
  rt = r;
  // compute the ramping gradient
  double T=1/getF();
  ramp = 2*getH()/( T*rt );
  // compute the max interval
  left=T/2*rt;
  right=T/2*( 1-rt );
}

double FieldACRamp::getField( double t ) const
{
  double tcyc = std::fmod( t, 1/( getF() ) );
  double T2 = 1/( 2*getF() );
  return tcyc<left ? tcyc*ramp :
    tcyc<right ? getH() :
    tcyc<( T2+left ) ? getH() - ramp*( tcyc-right ) :
    tcyc<( T2+right ) ? -getH() :
    -getH() + ramp*( tcyc-T2-right );
}

FieldFourierSquare::FieldFourierSquare(
    const double amplitude, const double frequency,
    const unsigned int _n_components
    )
    : FieldPeriodic( amplitude, frequency )
    , n_components( _n_components ) {}

double FieldFourierSquare::getField( double t ) const
{
    double field = 0;
    for( unsigned int k=1; k<n_components+1; ++k )
    {
        field += std::sin( 2*M_PI*( 2*k - 1 )*getF()*t )
            / ( 2*k-1 );
    }
    field *= 4 / M_PI * getH();
    return field;
}
