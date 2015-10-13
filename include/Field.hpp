// Field.hpp - Functions to compute different shaped fields.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#ifndef FIELD_H
#define FIELD_H

#define _USE_MATH_DEFINES
#include<cmath>


class Field
{
public:
  Field( float );

  void setH( float );
  float getH() const;
  virtual float getField( float ) const = 0;
  
private:
  float h;
};

class FieldPeriodic : public Field
{
public:
  FieldPeriodic( float h, float f );
  void setF( float );
  float getF() const;
  virtual float getField( float ) const = 0;
private:
  float f;
};

class FieldACSine : public FieldPeriodic
{
public:
  FieldACSine( float h, float f );
  virtual float getField( float ) const;
};

class FieldACCosine : public FieldPeriodic
{
public:
  FieldACCosine( float h, float f );
  virtual float getField( float ) const;
};

class FieldACSquare : public FieldPeriodic
{
public:
  FieldACSquare( float h, float f );
  virtual float getField( float ) const;
};

#endif
