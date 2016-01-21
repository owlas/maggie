// KramersTrig.cpp - a namespace containing functions to compute the
// transition rates for non-axially symmetric Stoner-Wohlfarth like
// magnetic nanoparticles.
//
// Functions are taken from the equations presented in Kalmykov, JAP,
// 96, 1138 (2004) - DOI: 10.1063/1.1760839
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<KramersTrig.hpp>

// Taylor expansion for the first activation energy
double KramersTrig::k_ebar_1( double s, double h, double psi )
{
  double res;
  res = 1;
  res += 2*h*( cos( psi ) - sin( psi ) );
  res += pow( h,2 );
  res += 0.5*pow( h,3 )*sin( 2*psi )*( cos( psi )-sin( psi ) );
  res += 0.5*pow( h,4 )*pow( sin( 2*psi ), 2 );
  res -= ( 1.0/32.0 )*pow( h,5 )*sin( 2*psi )*
    ( 7*sin( psi )+3*sin( 3*psi )-7*cos( psi )+3*cos( 3*psi ) );
  res += 0.5*pow( h,6 )*pow( sin( 2*psi ), 2 );
  res*= s;
  return res;
}

// Taylor expansion for the second activation energy
double KramersTrig::k_ebar_2( double s, double h, double psi )
{
  double res;
  res = 1;
  res -= 2*h*( cos( psi ) + sin( psi ) );
  res += pow( h,2 );
  res += 0.5*pow( h,3 )*sin( 2*psi )*( cos( psi )+sin( psi ) );
  res += 0.5*pow( h,4 )*pow( sin( 2*psi ), 2 );
  res += ( 1.0/32.0 )*pow( h,5 )*sin( 2*psi )*
    ( 7*sin( psi )+3*sin( 3*psi )+7*cos( psi )-3*cos( 3*psi ) );
  res += 0.5*pow( h,6 )*pow( sin( 2*psi ), 2 );
  res *= s;
  return res;
}

// Taylor expansion for the angular terms [zeroth]
double KramersTrig::k_angular_0( double s, double h, double psi,
                                double gamma, double Ms )
{
  double res;
  res = 2;
  res -= h*sin( psi );
  res -= pow( h,2 )*( 3+cos( 2*psi ) )/8.0;
  res -= pow( h,3 )*( 19+17*cos( 2*psi ) )*sin( psi )/16.0;
  res -= pow( h,4 )*( 351+28*cos( 2*psi )-283*cos( 4*psi ) )/512.0;
  /*res -= pow( h,5 )*( 2021+1332*cos( 2*psi )-663*cos( 4*psi ) )
   *sin( psi )/1024.0;*/
  res *= s*( gamma/Ms )*sqrt( h*sin( psi ) );
  return res;
}

// Taylor expansion for the angular terms [first]
double KramersTrig::k_angular_1( double s, double h, double psi,
                                double gamma, double Ms )
{
  double res;
  res = 2;
  res += 2*h*cos( psi );
  res -= pow( h,2 )*pow( sin( psi ),2 );
  res += 3*pow( h,3 )*cos( psi )*pow( sin( psi ),2 );
  res -= pow( h,4 )*( 21+19*cos( 2*psi ) )*pow( sin( psi ),2 )/8.0;
  res *= s*( gamma/Ms );
  return res;
}

// Taylor expansion for the angular terms [second]
double KramersTrig::k_angular_2( double s, double h, double psi,
                                double gamma, double Ms )
{
  double res;
  res = 2;
  res -= 2*h*cos( psi );
  res -= pow( h,2 )*pow( sin( psi ),2 );
  res -= 3*pow( h,3 )*cos( psi )*pow( sin( psi ),2 );
  res -= pow( h,4 )*( 21+19*cos( 2*psi ) )*pow( sin( psi ),2 )/8.0;
  res *= s*gamma/Ms;
  return res;
}


// Taylor expansion for the coefficient 1
double KramersTrig::k_coef_1( double s, double h, double psi )
{
  double res;
  res = 1;
  res += pow( h,2 )*pow( cos( psi ),2 )
    *( 0.5+h*sin( psi )+3*pow( h,2 )*( 5-3*cos( 2*psi ) )/16.0 );
  res *= 2*s*h*sin( psi );
  return res;
}

// Taylor expansion for the coefficient 2
double KramersTrig::k_coef_2( double s, double h, double psi )
{
  double res;
  res = -1;
  res += h*sin( psi );
  res += pow( h,2 )*pow( cos( psi ),2 );
  res += 5*pow( h,3 )*sin( 2*psi )*cos( psi )/4.0;
  res += pow( h,4 )*pow( sin( 2*psi ),2 );
  res *= 2*s;
  return res;
}

// Neel relaxation time constant
double KramersTrig::k_taun( double gamma, double Ms, double alpha, double V,
              double T )
{
    double beta = ( V/Constants::KB_d ) * ( 1/T );
    double res = beta * ( Ms/gamma );
    res *= ( 1+pow( alpha,2 ) )/( 2*alpha );
  return res;
}

// Damped saddle angular frequency
double KramersTrig::k_damped_angular_0( double s, double h, double psi, double gamma,
                                       double Ms, double alpha, double V, double T )
{
  double taun = KramersTrig::k_taun( gamma, Ms, alpha, V, T );
  double c1 = KramersTrig::k_coef_1( s, h, psi );
  double c2 = KramersTrig::k_coef_2( s, h, psi );
  double res = -c1-c2+sqrt( pow( c2-c1, 2 ) -4*( 1/pow( alpha,2 ) )*c1*c2 );
  res *= 0.25;
  res /= taun;
  return res;
}

// Kalmykov's expression for the IHD escape rate with oblique field
double KramersTrig::ihd_rate_1( double s, double h, double psi, double gamma, double Ms,
                               double alpha, double V, double T )
{
  double W0 = KramersTrig::k_damped_angular_0( s, h, psi, gamma, Ms, alpha, V, T );
  double w1 = KramersTrig::k_angular_1( s, h, psi, gamma, Ms );
  double w0 = KramersTrig::k_angular_0( s, h, psi, gamma, Ms );
  double bDelV = KramersTrig::k_ebar_1( s, h, psi );

  return W0*w1/( 2*M_PI*w0 )*exp( -bDelV );
}

// Kalmykov's expression for the IHD escape rate with oblique field
double KramersTrig::ihd_rate_2( double s, double h, double psi, double gamma, double Ms,
                               double alpha, double V, double T )
{
  double W0 = KramersTrig::k_damped_angular_0( s, h, psi, gamma, Ms, alpha, V, T );
  double w2 = KramersTrig::k_angular_2( s, h, psi, gamma, Ms );
  double w0 = KramersTrig::k_angular_0( s, h, psi, gamma, Ms );
  double bDelV = KramersTrig::k_ebar_2( s, h, psi );

  return W0*w2/( 2*M_PI*w0 )*exp( -bDelV );
}

// Kalmykov's expressions for the energy maxima and minima
double KramersTrig::theta_max( double h, double psi )
{
  return acos( -h*cos(psi)
	       -pow( h,2 )*sin(psi)*cos(psi)
	       *( 1+h*sin(psi)
		  + pow( h,2 )/4.0*( 3-cos(2*psi) )
		  + pow( h,3 )/2.0*sin(psi)*( 2+cos(2*psi) )
		  + pow( h,4 )/64.0*(73-20*cos(2*psi)
				     -29*cos(4*psi) ) ) );
}

double KramersTrig::theta_min1( double h, double psi )
{
  return acos( 1 - pow( h,2 )/2.0*pow( sin(psi),2 )
	       *( 1-2*h*cos(psi)
		  + pow( h,2 )/8.0*( 13+11*cos(2*psi) )
		  - pow( h,3 )*( 3+cos(2*psi)*cos(psi) )
		  + pow( h,4 )/64.0*( 184+156*cos(2*psi)
				      -19*cos(4*psi) ) ) );
}

double KramersTrig::theta_min2( double h, double psi )
{
  return acos( -1 + pow( h,2 )/2.0*pow( sin(psi),2 )
	       *( 1+2*h*cos(psi)
		  + pow( h,2 )/8.0*( 13+11*cos(2*psi) )
		  + pow( h,3 )*( 3+cos(2*psi) )*cos(psi)
		  + pow( h,4 )/64.0*( 184+156*cos(2*psi)
				      -19*cos(4*psi) ) ) );
}
