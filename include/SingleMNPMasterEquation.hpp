// SingleMNPMasterEquation.hpp - Class to represent the master
// equation for a single uniaxial magnetic nanoparticle. Derived from
// LangevunEquation base.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#ifndef MNPME_H
#define MNPME_H

#include<ODE.hpp>
#include<Constants.hpp>
#include<Field.hpp>
#include<KramersTrig.hpp>
#include<stdexcept>
#include<functional>
#include <utility>
using dpair = std::pair<double, double>;
#include <cmath>
using std::abs;

class SingleMNPMasterEquation : public ODE<2>
{
public:
    // constructor for the general case
    SingleMNPMasterEquation( double anis,
                             double temp,
                             const Field &field,
                             double fieldangle,
                             double radius,
                             double alpha,
                             double Ms,
                             double gamma=Constants::GYROMAG );

    // delegated constructor for aligned case
    SingleMNPMasterEquation( double anis,
                             double temp,
                             const Field &field,
                             double radius );



  // compute the difference in energy between the two states of the
  // system
  double ediff( double t );

  // compute the angle of the minima
  dpair state_rotation( const double t ) const;

  // Compute the drift of the master equation
    void computeDrift( ODE<2>::array& out, const ODE<2>::array& in, const double t ) const;

private:
  double k;
  double T;
  double r;
  double psi;
  double a;
  double M;
  double gam;
  double v;
  const Field *fieldPtr;
};

#endif
