// SDES.cpp -- implementation of stochastic differential equations as a
// Langevin Equation object.
// 
// Oliver Laslett <O.Laslett@soton.ac.uk>
// 2015
// 

#include <SDEs.hpp>

OH::OH( const float t, const float m, const float s )
  : SDE( 1,1 )
  , theta( t )
  , mu( m )
  , sigma( s )
{
  // empty
}
float OH::getTheta() const { return theta; }
float OH::getMu() const { return mu; }
float OH::getSigma() const { return sigma; }

void OH::computeDrift( array_f& out, const array_f& in, const float ) const
{
  out[0] = theta*( mu-in[0] );
}

void OH::computeDiffusion( matrix_f& out, const array_f&, const float ) const
{
  out[0][0] = sigma;
}

Wiener::Wiener()
  : SDE( 1,1 )
{
  // empty
}

void Wiener::computeDrift( array_f& out, const array_f&, const float ) const
{
  out[0] = 0;
}

void Wiener::computeDiffusion( matrix_f& out, const array_f&, const float )
  const
{
  out[0][0] = 1;
}

ODEConstantDrift::ODEConstantDrift( const float drift )
  : ODE( 1 )
  , a( drift )
{
  //empty
}

void ODEConstantDrift::computeDrift( array_f& out, const array_f&,
				     const float ) const
{
  out[0] = a;
}

SDE_AXpB::SDE_AXpB( const float drift, const float diff )
  : SDE( 1,1 )
  , a( drift )
  , b( diff )
{}

void SDE_AXpB::computeDrift( array_f& out, const array_f& in, const float )
  const
{
  out[0] = a*in[0];
}

void SDE_AXpB::computeDiffusion( matrix_f& out, const array_f&, const float )
  const
{
  out[0][0] = b;
}

SDE_AXpBX::SDE_AXpBX( const float drift, const float diff )
  : SDE( 1,1 )
  , a( drift )
  , b( diff )
{}

void SDE_AXpBX::computeDrift( array_f& out, const array_f& in, const float )
  const
{
  out[0] = a*in[0];
}

void SDE_AXpBX::computeDiffusion( matrix_f& out, const array_f& in,
				  const float ) const
{
  out[0][0] = b*in[0];
}

MultiplicativeConstantNoise::MultiplicativeConstantNoise( const float drift,
							  const float diffusion )
  : SDE( 1,1 )
  , a( drift )
  , b( diffusion )
{
  //empty
}

void MultiplicativeConstantNoise::computeDrift( array_f& out,
						const array_f& in,
						const float ) const
{
  out[0] = a*in[0];
}

void MultiplicativeConstantNoise::computeDiffusion( matrix_f& out,
						    const array_f& in,
						    const float ) const
{
  out[0][0] = b*in[0];
}

void MultiplicativeConstantNoise::computeDiffusionDerivatives( array3_f& out,
							       const array_f&,
							       const float )
  const
{
  out[0][0][0] = b;
}
