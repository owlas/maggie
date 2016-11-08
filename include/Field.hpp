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
  Field( double );

  void setH( double );
  double getH() const;
  virtual double getField( double ) const = 0;

private:
  double h;
};

class FieldConstant : public Field
{
public:
  FieldConstant( double h );
  virtual double getField( double ) const;
};

class FieldPeriodic : public Field
{
public:
  FieldPeriodic( double h, double f );
  void setF( double );
  double getF() const;
  virtual double getField( double ) const = 0;
private:
  double f;
};

class FieldACSine : public FieldPeriodic
{
public:
  FieldACSine( double h, double f );
  virtual double getField( double ) const;
};

class FieldACCosine : public FieldPeriodic
{
public:
  FieldACCosine( double h, double f );
  virtual double getField( double ) const;
};

class FieldACSquare : public FieldPeriodic
{
public:
  FieldACSquare( double h, double f );
  virtual double getField( double ) const;
};

class FieldACRamp : public FieldPeriodic
{
public:
  FieldACRamp( double h, double f, double rTime );
  virtual double getField( double ) const;
  double getRT() const;
  void setRT( double );
private:
  double rt;
  double left;
  double right;
  double ramp;
};

class FieldFourierSquare : public FieldPeriodic
{
public:
    FieldFourierSquare( const double amplitude,
                        const double frequency,
                        const unsigned int n_components);
    virtual double getField( double ) const;
private:
    unsigned int n_components;
};

#endif
