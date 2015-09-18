// LangevinEquation.cpp
// Implementation for stochastic LLG Langevin Equation
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include <iostream>
using std::cout;
using std::endl;

#include <StocLLG.hpp>

// Constructor
StocLLG::StocLLG( float s, float a, float hx, float hy, float hz )
  : LangevinEquation( 3 )
  , h( boost::extents[3] )
{
  setSigma( s );
  setAlpha( a );
  setReducedHeff( hx, hy, hz );
}

// set the reduced effective field
void StocLLG::setReducedHeff( float hx, float hy, float hz )
{
  h[0] = hx;
  h[1] = hy;
  h[2] = hz;
}
// get the reduced effective field
array_f& StocLLG::getReducedHeff() { return h; }

// set and set sigma
void StocLLG::setSigma( float s ) { sigma = s; }
float StocLLG::getSigma() { return sigma; }

// set and set alpha
void StocLLG::setAlpha( float s ) { alpha = s; }
float StocLLG::getAlpha() { return alpha; }

// compute the drift term of the stochastic LLG equation
void StocLLG::computeDrift( array_f& out, array_f& state )
{
  out[0] = state[2]*h[1] - state[1]*h[2]
    + alpha*(h[0]*(state[1]*state[1] + state[2]*state[2])
             - state[0]*(state[1]*h[1]
                         + state[2]*h[2]));
  out[1] = state[0]*h[2] - state[2]*h[0]
    + alpha*(h[1]*(state[0]*state[0] + state[2]*state[2])
             - state[1]*(state[0]*h[0]
                         + state[2]*h[2]));
  out[2] = state[1]*h[0] - state[0]*h[1]
    + alpha*(h[2]*(state[0]*state[0] + state[1]*state[1])
             - state[2]*(state[0]*h[0]
                         + state[1]*h[1]));
}

// compute the diffusion matrix for the stochastic LLG equation
void StocLLG::computeDiffusion( matrix_f& out, array_f& state)
{
  out[0][0] = alpha*sigma*(state[1]*state[1]+state[2]*state[2]);
  out[0][1] = sigma*(state[2]-alpha*state[0]*state[1]);
  out[0][2] = -sigma*(state[1]+alpha*state[0]*state[2]);
  out[1][0] = -sigma*(state[2]+alpha*state[0]*state[1]);
  out[1][1] = alpha*sigma*(state[0]*state[0]+state[2]*state[2]);
  out[1][2] = sigma*(state[0]-alpha*state[1]*state[2]);
  out[2][0] = sigma*(state[1]-alpha*state[0]*state[2]);
  out[2][1] = -sigma*(state[0]+alpha*state[1]*state[2]);
  out[2][2] = alpha*sigma*(state[0]*state[0]+state[1]*state[1]);
}
