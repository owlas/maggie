// Field.hpp - Functions to compute different shaped fields.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#ifndef FIELD_H
#define FIELD_H

#define _USE_MATH_DEFINES
#include<cmath>
#include<stdexcept>

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

class FieldACRamp : public FieldPeriodic
{
public:
  FieldACRamp( float h, float f, float rTime );
  virtual float getField( float ) const;
  float getRT() const;
  void setRT( float );
private:
  float rt;
  float left;
  float right;
  float ramp;
};

#endif
