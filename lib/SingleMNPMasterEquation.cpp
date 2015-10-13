// SingleMNPMasterEquation.cpp - Implementation of the master equation
// for a single uniaxial magnetic nanoparticle.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<SingleMNPMasterEquation.hpp>

SingleMNPMasterEquation::
SingleMNPMasterEquation(float temperature,
                        float anisotropy_constant,
                        float particle_radius,
                        const Field &field,
                        float field_angle )
  : LangevinEquation( 2 )
  , k( anisotropy_constant )
  , T( temperature )
  , r( particle_radius )
  , psi( field_angle )
{
  // get the field pointer
  fieldPtr = &field;
  // set volume of particle
  v = 4.0/3.0*M_PI*std::pow( r,3 );
}

void SingleMNPMasterEquation::computeDrift( array_f& out, array_f& in, float t )
{
  float rate1, rate2;

  if( fieldPtr->getField( t )*sin( psi )<0.03 )


    {
      float ebar1 = k*v*std::pow( 1-fieldPtr->getField( t ), 2 );
      float ebar2 = k*v*std::pow( 1+fieldPtr->getField( t ), 2 );
      rate1 = 1/( TAU0*std::exp( ebar1/( Constants::KB*T ) ) );
      rate2 = 1/( TAU0*std::exp( ebar2/( Constants::KB*T ) ) );
    }
  /*else
    {
    float sigma = k*v/( Constants::KB*T );
    rate1 = KramersTrig::ihd_rate_1( sigma, fieldPtr->getField( t ), psi );
    rate2 = KramersTrig::ihd_rate_2( sigma, fieldPtr->getField( t ), psi );
    }*/

  out[0] = -( rate1+rate2 )*in[0] + rate2;
  out[1] = -( rate1+rate2 )*in[1] + rate1;
}

float SingleMNPMasterEquation::ediff( float t )
{
  if( fieldPtr->getField( t )*sin( psi ) < 0.03 )
    return -4*k*v*fieldPtr->getField( t );
  else
    return 0;
  /*else
    {
    float beta = v/( Constants::KB*T );
    float sigma = k*beta;
    // Multiply by K_b*T because of Eq.6 in Kalmykov paper
    // energy diff functions are reduced by this.
    return Constants::KB*T*( KrmaersTrig::k_ebar_1( sigma,
    fieldPtr->getField( t ),
    psi )
    - KramersTrig::k_ebar_2( sigma,
    fieldPtr->getField( t ),
    psi ) );
    }*/
}

