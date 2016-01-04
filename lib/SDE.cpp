// SDE.cpp -- implementation of stochastic differential equations as a
// Langevin Equation object.
// 
// Oliver Laslett <O.Laslett@soton.ac.uk>
// 2015
// 

#include <SDE.hpp>

OH::OH( float t, float m, float s )
  : LangevinEquation( 1 )
  , theta( t )
  , mu( m )
  , sigma( s )
{
  // empty
}
float OH::getTheta() const { return theta; }
float OH::getMu() const { return mu; }
float OH::getSigma() const { return sigma; }

void OH::computeDrift( array_f& out, array_f& in, float )
{
  out[0] = theta*( mu-in[0] );
}

void OH::computeDiffusion( matrix_f& out, array_f&, float )
{
  out[0][0] = sigma;
}

Wiener::Wiener()
  : LangevinEquation( 1 )
{
  // empty
}

void Wiener::computeDrift( array_f& out, array_f&, float )
{
  out[0] = 0;
}

void Wiener::computeDiffusion( matrix_f& out, array_f&, float )
{
  out[0][0] = 1;
}

ODEConstantDrift::ODEConstantDrift( float drift )
  : LangevinEquation( 1 )
  , a( drift )
{
  //empty
}

SDE_AXpB::SDE_AXpB( float drift, float diff )
  : LangevinEquation( 1 )
  , a( drift )
  , b( diff )
{}

void SDE_AXpB::computeDrift( array_f& out, array_f& in, float )
{
  out[0] = a*in[0];
}

void SDE_AXpB::computeDiffusion( matrix_f& out, array_f&, float )
{
  out[0][0] = b;
}

SDE_AXpBX::SDE_AXpBX( float drift, float diff )
  : LangevinEquation( 1 )
  , a( drift )
  , b( diff )
{}

void SDE_AXpBX::computeDrift( array_f& out, array_f& in, float )
{
  out[0] = a*in[0];
}

void SDE_AXpBX::computeDiffusion( matrix_f& out, array_f& in, float )
{
  out[0][0] = b*in[0];
}

void ODEConstantDrift::computeDrift( array_f& out, array_f&, float )
{
  out[0] = a;
}

MultiplicativeConstantNoise::MultiplicativeConstantNoise( float drift,
							  float diffusion )
  : LangevinEquation( 1 )
  , a( drift )
  , b( diffusion )
{
  //empty
}

void MultiplicativeConstantNoise::computeDrift( array_f& out, array_f& in,
						float )
{
  out[0] = a*in[0];
}

void MultiplicativeConstantNoise::computeDiffusion( matrix_f& out,
						    array_f& in, float )
{
  out[0][0] = b*in[0];
}

void MultiplicativeConstantNoise::computeDiffusionDerivatives( array3_f& out,
							       array_f&,
							       float )
{
  out[0][0][0] = b;
}
