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
template <typename T>
StocLLG<T>::StocLLG( const T s, const T a,
		  const T hx, const T hy, const T hz )
    : SDE<T>( 3,3 )
  , h( boost::extents[3] )
{
  setSigma( s );
  setAlpha( a );
  setReducedHeff( hx, hy, hz );
}

// set the reduced effective field
template <typename T>
void StocLLG<T>::setReducedHeff( const T hx, T hy, const T hz )
{
  h[0] = hx;
  h[1] = hy;
  h[2] = hz;
}
// get the reduced effective field
template <typename T>
array<T> StocLLG<T>::getReducedHeff() const { return h; }

// set and set sigma
template <typename T>
void StocLLG<T>::setSigma( const T s ) { sigma = s; }
template <typename T>
T StocLLG<T>::getSigma() const { return sigma; }

// set and set alpha
template <typename T>
void StocLLG<T>::setAlpha( const T s ) { alpha = s; }
template <typename T>
T StocLLG<T>::getAlpha() const { return alpha; }

// get a modifiable reference to the reduced field
template <typename T>
array<T>& StocLLG<T>::getReducedFieldRef() { return h; }

// compute the drift term of the stochastic LLG equation
template <typename T>
void StocLLG<T>::computeDrift( array<T>& out, const array<T>& state,
			    const T ) const
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
template <typename T>
void StocLLG<T>::computeDiffusion( matrix<T>& out, const array<T>& state,
				const T ) const
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
template <typename T>
void StocLLG<T>::computeDiffusionDerivatives( array3<T> &out,
                                              const array<T> &state,
                                              const T ) const
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
template <typename T>
StocLLGIto<T>::StocLLGIto( const T s, const T a, const T hx,
                           const T hy, const T hz )
    : StocLLG<T>( s, a, hx, hy, hz )
{}

// compute the drift in ito form
template <typename T>
void StocLLGIto<T>::computeDrift(array<T>& out, const array<T>& in,
                                 const T t) const
{
    const int d{ this->getDim() };
    const int m{ this->getWDim() };
    matrix<T> b( boost::extents[d][m] );
    array3<T> bDash( boost::extents[d][m][d]);

    StocLLG<T>::computeDrift(out, in, t);
    this->computeDiffusion(b, in, t );
    this->computeDiffusionDerivatives(bDash, in, t );

    for ( bidx i=0; i!=d; ++i )
        for( bidx j=0; j!=m; ++j )
            for( bidx k=0; k!=d; ++k )
                out[i] += 0.5 * b[k][j] * bDash[i][j][k];
}

// Explicit template class instatiation
template class StocLLG<float>;
template class StocLLG<double>;
template class StocLLGIto<float>;
template class StocLLGIto<double>;
