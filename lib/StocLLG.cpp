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
StocLLG::StocLLG( const float s, const float a,
		  const float hx, const float hy, const float hz )
  : SDE( 3,3 )
  , h( boost::extents[3] )
{
  setSigma( s );
  setAlpha( a );
  setReducedHeff( hx, hy, hz );
}

// set the reduced effective field
void StocLLG::setReducedHeff( const float hx, const float hy, const float hz )
{
  h[0] = hx;
  h[1] = hy;
  h[2] = hz;
}
// get the reduced effective field
array_f StocLLG::getReducedHeff() const { return h; }

// set and set sigma
void StocLLG::setSigma( const float s ) { sigma = s; }
float StocLLG::getSigma() const { return sigma; }

// set and set alpha
void StocLLG::setAlpha( const float s ) { alpha = s; }
float StocLLG::getAlpha() const { return alpha; }

// compute the drift term of the stochastic LLG equation
void StocLLG::computeDrift( array_f& out, const array_f& state,
			    const float ) const
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
void StocLLG::computeDiffusion( matrix_f& out, const array_f& state,
				const float ) const
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

// compute the Taylor derivative sums L
void StocLLG::computeDiffusionDerivatives( array3_f &out,
					   const array_f &state,
					   const float ) const
{
  out[0][0][0] = 0;
  out[0][0][1] = 2*alpha*sigma*state[1];
  out[0][0][2] = 2*alpha*sigma*state[2];

  out[0][1][0] = -alpha*sigma*state[1];
  out[0][1][1] = -alpha*sigma*state[2];
  out[0][1][2] = sigma;

  out[0][2][0] = -alpha*sigma*state[2];
  out[0][2][1] = -sigma;
  out[0][2][2] = -alpha*sigma*state[0];

  out[1][0][0] = -alpha*sigma*state[1];
  out[1][0][1] = -alpha*sigma*state[0];
  out[1][0][2] = -sigma;

  out[1][1][0] = 2*alpha*sigma*state[0];
  out[1][1][1] = 0;
  out[1][1][2] = 2*alpha*sigma*state[2];

  out[1][2][0] = sigma;
  out[1][2][1] = -alpha*sigma*state[2];
  out[1][2][2] = -alpha*sigma*state[1];

  out[2][0][0] = -alpha*sigma*state[2];
  out[2][0][1] = sigma;
  out[2][0][2] = -alpha*sigma*state[0];

  out[2][1][0] = -sigma;
  out[2][1][1] = -alpha*sigma*state[2];
  out[2][1][2] = -alpha*sigma*state[1];

  out[2][2][0] = 2*alpha*sigma*state[0];
  out[2][2][1] = 2*alpha*sigma*state[2];
  out[2][2][2] = 0;
}

// constructor for the stoc llg
StocLLGIto::StocLLGIto( const float s, const float a, const float hx,
                        const float hy, const float hz )
    : StocLLG( s, a, hx, hy, hz )
{}

// compute the drift in ito form
void StocLLGIto::computeDrift(array_f& out, const array_f& in,
                              const float t) const
{
    const int d{ getDim() };
    const int m{ getWDim() };
    matrix_f b( boost::extents[d][m] );
    cube_f bDash( boost::extents[d][m][d]);

    StocLLG::computeDrift(out, in, t);
    computeDiffusion(b, in, t );
    computeDiffusionDerivatives(bDash, in, t );

    for ( array_f::index i=0; i!=d; ++i )
        for( array_f::index j=0; j!=m; ++j )
            for( array_f::index k=0; j!=d; ++k )
                out[i] += 0.5 * b[k][j] * bDash[i][j][k];
}
