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

// Sine field
FieldACSine::FieldACSine( float field, float freq )
  : Field( field )
{
  setF( freq );
}
void FieldACSine::setF( float freq ) { f=freq; }
float FieldACSine::getF() const { return f; }
float FieldACSine::getField( float t ) const
{
  return getH()*std::sin( 2*M_PI*getF()*t );
}

// Sine field
FieldACCosine::FieldACCosine( float field, float freq )
  : Field( field )
{
  setF( freq );
}
void FieldACCosine::setF( float freq ) { f=freq; }
float FieldACCosine::getF() const { return f; }
float FieldACCosine::getField( float t ) const
{
  return getH()*std::cos( 2*M_PI*getF()*t );
}

// Sine field
FieldACSquare::FieldACSquare( float field, float freq )
  : Field( field )
{
  setF( freq );
}
void FieldACSquare::setF( float freq ) { f=freq; }
float FieldACSquare::getF() const { return f; }
float FieldACSquare::getField( float t ) const
{
  return getH()*( int( t*getF()*2 )%2 ? -1 : 1 );
}
