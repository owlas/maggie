// SingleMNPMasterEquation.cpp - Implementation of the master equation
// for a single uniaxial magnetic nanoparticle.
//
// Oliver Laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<SingleMNPMasterEquation.hpp>
#include<iostream>

SingleMNPMasterEquation::
SingleMNPMasterEquation( float anis,
                         float temp,
                         const Field &field,
                         float fieldangle,
                         float radius,
                         float alpha,
                         float Ms,
                         float gamma )
    : ODE<float>( 2 )
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
                                            const array_f& in,
                                            const float t ) const
{
    float rate1, rate2;
    float h = fieldPtr->getField( t );
    /* For weak fields and low angles approximate as aligned
       and use the Neel-Arrhenius laws */
    if( abs( h )*sin( psi )<0.03 )
    {
        float ebar1 = k*v*std::pow( 1-h, 2 );
        float ebar2 = k*v*std::pow( 1+h, 2 );
        rate2 = 1/( 1e-10*std::exp( ebar1/( Constants::KB*T ) ) );
        rate1 = 1/( 1e-10*std::exp( ebar2/( Constants::KB*T ) ) );
    }
    /* Otherwise use the Kramers theory approach
       This does not work for negative fields, so we exploit the symmetry to
       solve the problem */
    else
    {
        float sigma = k*v/( Constants::KB*T );
        rate1 = KramersTrig::ihd_rate_1( sigma, abs( h ), psi, gam, M,
                                         a, v, T );
        rate2 = KramersTrig::ihd_rate_2( sigma, abs( h ), psi, gam, M,
                                             a, v, T );
        if( h<0 )
        {
            float temp = rate1;
            rate1 = rate2;
            rate2 = temp;
        }
    }

    out[0] = -( rate1+rate2 )*in[0] + rate2;
    out[1] = -( rate1+rate2 )*in[1] + rate1;
}

float SingleMNPMasterEquation::ediff( float t )
{
    float h = fieldPtr->getField( t );
    if( abs( h )*sin( psi ) < 0.03 )
        return -4*k*v*h;
    else
    {
        float beta = v/( Constants::KB*T );
        float sigma = k*beta;
        // Multiply by K_b*T because of Eq.6 in Kalmykov paper
        // energy diff functions are reduced by this.
        return
            Constants::KB*T*(   - KramersTrig::k_ebar_1( sigma, h,
                                                         psi )
                                + KramersTrig::k_ebar_2( sigma, h,
                                                         psi ) );
    }
}

pair SingleMNPMasterEquation::state_rotation( float t ) const
{
    float t1, t2;
    float h = fieldPtr->getField( t );
    if( abs( h )*sin( psi ) < 0.03 )
    {
        t1=1;
        t2=-1;
    }
    else
    {
        t1 = KramersTrig::theta_min1( h, psi );
        t2 = KramersTrig::theta_min2( h, psi );
    }
    pair ans = { t1, t2 };
    return ans;
}
