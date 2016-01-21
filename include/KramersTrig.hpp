// KramersTrig.hpp - header file for KramersTrig functions
// see implementation for details.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#ifndef KRAM_H
#define KRAM_H

#define _USE_MATH_DEFINES

#include<maggie.hpp>
#include<cmath>
using std::sin; using std::cos;
using std::pow; using std::sqrt;
using std::acos;

namespace KramersTrig
{
  double k_ebar_1( double s, double h, double psi );
  double k_ebar_2( double s, double h, double psi );
  double k_angular_0( double s, double h, double psi,
                     double gamma, double Ms );
  double k_angular_1( double s, double h, double psi,
                     double gamma, double Ms );
  double k_angular_2( double s, double h, double psi,
                     double gamma, double Ms );
  double k_coef_1( double s, double h, double psi );
  double k_coef_2( double s, double h, double psi );
  double k_taun( double gamma, double Ms, double alpha, double V,
                double T );
  double k_damped_angular_0( double s, double h, double psi, double gamma,
                            double Ms, double alpha, double V, double T );
  double ihd_rate_1( double s, double h, double psi, double gamma, double Ms,
                    double alpha, double V, double T );
  double ihd_rate_2( double s, double h, double psi, double gamma, double Ms,
                    double alpha, double V, double T );
  double theta_max( double h, double psi );
  double theta_min1( double h, double psi );
  double theta_min2( double h, double psi );
}

#endif
