// SingleMNPMasterEquation.cpp - Implementation of the master equation
// for a single uniaxial magnetic nanoparticle.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<SingleMNPMasterEquation.hpp>

SingleMNPMasterEquation::
SingleMNPMasterEquation( float anis,
			 float temp,
			 const Field &field,
			 float fieldangle,
			 float radius,
			 float alpha,
			 float Ms,
			 float gamma )
  : ODE( 2 )
  , k( anis )
  , T( temp )
  , r( radius )
  , psi( fieldangle )
  , a( alpha )
  , M( Ms )
  , gam( gamma )
{
  fieldPtr = &field;
  v = 4.0/3.0*M_PI*std::pow( r,3 );
}

// Delegate constructor for the aligned case
SingleMNPMasterEquation::
SingleMNPMasterEquation( float anis,
			 float temp,
			 const Field &field,
			 float radius )
  : SingleMNPMasterEquation( anis, temp, field, 0, radius,
			     0.1, 100 ) {}

void SingleMNPMasterEquation::computeDrift( array_f& out,
					    array_f& in,
					    float t )
{
  float rate1, rate2;
  float h = fieldPtr->getField( t );
  if( h*sin( psi )<0.03 )
    {
      float ebar1 = k*v*std::pow( 1-h, 2 );
      float ebar2 = k*v*std::pow( 1+h, 2 );
      rate1 = 1/( 1e-10*std::exp( ebar1/( Constants::KB*T ) ) );
      rate2 = 1/( 1e-10*std::exp( ebar2/( Constants::KB*T ) ) );
    }
  else
    {
    float sigma = k*v/( Constants::KB*T );
    rate1 = KramersTrig::ihd_rate_1( sigma, h, psi, gam, M,
				     a, v, T );
    rate2 = KramersTrig::ihd_rate_2( sigma, h, psi, gam, M,
				     a, v, T );
    }

  out[0] = -( rate1+rate2 )*in[0] + rate2;
  out[1] = -( rate1+rate2 )*in[1] + rate1;
}

float SingleMNPMasterEquation::ediff( float t )
{
  float h = fieldPtr->getField( t );
  if( h*sin( psi ) < 0.03 )
    return -4*k*v*h;
  else
    {
    float beta = v/( Constants::KB*T );
    float sigma = k*beta;
    // Multiply by K_b*T because of Eq.6 in Kalmykov paper
    // energy diff functions are reduced by this.
    return
      Constants::KB*T*(   KramersTrig::k_ebar_1( sigma, h,
						 psi )
			- KramersTrig::k_ebar_2( sigma, h,
						 psi ) );
    }
}

float SingleMNPMasterEquation::state_rotation( float t ) const 
{
  float h = fieldPtr->getField( t );
  if( h*sin( psi ) < 0.03 )
    return 0;
  else
    return KramersTrig::theta_min1( h, psi );
}
