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
