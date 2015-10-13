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

class FieldACSine : public Field
{
public:
  FieldACSine( float h, float f );
  void setF( float );
  float getF() const;
  virtual float getField( float ) const;

private:
  float f;
};

class FieldACCosine : public Field
{
public:
  FieldACCosine( float h, float f );
  void setF( float );
  float getF() const;
  virtual float getField( float ) const;

private:
  float f;
};

class FieldACSquare : public Field
{
public:
  FieldACSquare( float h, float f );
  void setF( float );
  float getF() const;
  virtual float getField( float ) const;

private:
  float f;
};

#endif
