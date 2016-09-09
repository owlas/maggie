// SDES.cpp -- implementation of stochastic differential equations as a
// Langevin Equation object.
//
// Oliver Laslett <O.Laslett@soton.ac.uk>
// 2015
//

#include "../include/SDEs.hpp"

OH::OH( const double t, const double m, const double s )
    : SDE<1,1>()
    , theta( t )
    , mu( m )
    , sigma( s )
{
    // empty
}
double OH::getTheta() const { return theta; }
double OH::getMu() const { return mu; }
double OH::getSigma() const { return sigma; }

void OH::computeDrift( OH::array& out, const OH::array& in, const double ) const
{
    out[0] = theta*( mu-in[0] );
}

void OH::computeDiffusion( OH::matrix& out, const OH::array&, const double ) const
{
    out[0][0] = sigma;
}

Wiener::Wiener()
    : SDE<1,1>()
{
    // empty
}

void Wiener::computeDrift( Wiener::array& out, const Wiener::array&, const double ) const
{
    out[0] = 0;
}

void Wiener::computeDiffusion( Wiener::matrix& out, const Wiener::array&, const double )
    const
{
    out[0][0] = 1;
}

ODEConstantDrift::ODEConstantDrift( const double drift )
    : ODE<1>()
    , a( drift )
{
    //empty
}

void ODEConstantDrift::computeDrift( ODEConstantDrift::array& out,
                                     const ODEConstantDrift::array&,
				     const double ) const
{
    out[0] = a;
}

SDE_AXpB::SDE_AXpB( const double drift, const double diff )
    : SDE<1,1>()
    , a( drift )
    , b( diff )
{}

void SDE_AXpB::computeDrift( SDE_AXpB::array& out,
                             const SDE_AXpB::array& in,
                             const double )
    const
{
    out[0] = a*in[0];
}

void SDE_AXpB::computeDiffusion( SDE_AXpB::matrix& out,
                                 const SDE_AXpB::array&,
                                 const double )
    const
{
    out[0][0] = b;
}

SDE_AXpBX::SDE_AXpBX( const double drift, const double diff )
    : SDE<1,1>()
    , a( drift )
    , b( diff )
{}

void SDE_AXpBX::computeDrift( SDE_AXpBX::array& out,
                              const SDE_AXpBX::array& in,
                              const double )
    const
{
    out[0] = a*in[0];
}

void SDE_AXpBX::computeDiffusion( SDE_AXpBX::matrix& out,
                                  const SDE_AXpBX::array& in,
				  const double ) const
{
    out[0][0] = b*in[0];
}

MultiplicativeConstantNoise::MultiplicativeConstantNoise( const double drift,
							  const double diffusion )
    : SDE<1,1>()
    , a( drift )
    , b( diffusion )
{
    //empty
}

void MultiplicativeConstantNoise::computeDrift( MultiplicativeConstantNoise::array& out,
						const MultiplicativeConstantNoise::array& in,
						const double ) const
{
    out[0] = a*in[0];
}

void MultiplicativeConstantNoise::computeDiffusion( MultiplicativeConstantNoise::matrix& out,
						    const MultiplicativeConstantNoise::array& in,
						    const double ) const
{
    out[0][0] = b*in[0];
}

void MultiplicativeConstantNoise::computeDiffusionDerivatives(
    MultiplicativeConstantNoise::array3& out,
    const MultiplicativeConstantNoise::array&,
    const double ) const
{
    out[0][0][0] = b;
}
