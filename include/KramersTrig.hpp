// KramersTrig.hpp - header file for KramersTrig functions
// see implementation for details.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#define _USE_MATH_DEFINES

#include<maggie.hpp>
#include<cmath>
using std::sin; using std::cos;
using std::pow; using std::sqrt;
using std::acos;

namespace KramersTrig
{
  float k_ebar_1( float s, float h, float psi );
  float k_ebar_2( float s, float h, float psi );
  float k_angular_0( float s, float h, float psi,
                     float gamma, float Ms );
  float k_angular_1( float s, float h, float psi,
                     float gamma, float Ms );
  float k_angular_2( float s, float h, float psi,
                     float gamma, float Ms );
  float k_coef_1( float s, float h, float psi );
  float k_coef_2( float s, float h, float psi );
  float k_taun( float gamma, float Ms, float alpha, float V,
                float T );
  float k_damped_angular_0( float s, float h, float psi, float gamma,
                            float Ms, float alpha, float V, float T );
  float ihd_rate_1( float s, float h, float psi, float gamma, float Ms,
                    float alpha, float V, float T );
  float ihd_rate_2( float s, float h, float psi, float gamma, float Ms,
                    float alpha, float V, float T );
  float theta_max( float h, float psi );
  float theta_min1( float h, float psi );
  float theta_min2( float h, float psi );
}
