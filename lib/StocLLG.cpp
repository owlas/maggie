// LangevinEquation.cpp
// Implementation for stochastic LLG Langevin Equation
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include "../include/StocLLG.hpp"

using namespace maggie;

// Constructor
StocLLG::StocLLG( const stability s, const damping a, const field h )
    : SDE<3,3>()
{
  setSigma( s );
  setAlpha( a );
  setReducedHeff( h );
}

// set the reduced effective field
void StocLLG::setReducedHeff( const field _h )
{
    h = _h;
}
// get the reduced effective field
field StocLLG::getReducedHeff() const { return h; }

// set and set sigma
void StocLLG::setSigma( const stability s ) { sigma = s; }
stability StocLLG::getSigma() const { return sigma; }

// set and set alpha
void StocLLG::setAlpha( const damping a ) { alpha = a; }
damping StocLLG::getAlpha() const { return alpha; }

// get a modifiable reference to the reduced field
field& StocLLG::getReducedFieldRef() { return h; }

// compute the drift term of the stochastic LLG equation
void StocLLG::computeDrift( array& out, const array& state,
			    const double ) const
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
void StocLLG::computeDiffusion( matrix& out, const array& state,
				const double ) const
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
void StocLLG::computeDiffusionDerivatives( array3 &out,
                                           const array& state,
                                           const double ) const
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

StocLLGIto::StocLLGIto( const stability s, const damping a, const field h )
    : StocLLG( s, a, h )
{}

// compute the drift in ito form
void StocLLGIto::computeDrift( array& out, const array& in,
                               const double t) const
{
    const int d{ StocLLGIto::dim };
    const int m{ StocLLGIto::wdim };
    StocLLGIto::matrix b;
    StocLLGIto::array3 bDash;

    StocLLG::computeDrift(out, in, t);
    this->computeDiffusion(b, in, t );
    this->computeDiffusionDerivatives(bDash, in, t );

    for ( unsigned int i=0; i!=d; ++i )
        for( unsigned int j=0; j!=m; ++j )
            for( unsigned int k=0; k!=d; ++k )
                out[i] += 0.5 * b[k][j] * bDash[i][j][k];
}
