// Heun.cpp
// Implementation for the Heun integration scheme
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<Heun.hpp>


// Constructor
Heun::Heun( const SDE &le, const array_f& init_state,
	    const float time, const float dt, mt19937& rng )
  : Integrator<SDE>( le, init_state, time)
  , h( dt )
  , dim( le.getDim() )
  , dw( boost::extents[dim])
  , xPred( boost::extents[dim])
  , tmp1( boost::extents[dim])
  , tmp1Up( boost::extents[dim])
  , tmp2( boost::extents[dim][dim] )
  , tmp2Up( boost::extents[dim][dim] )
  , dist(0,1)
  , gen( rng, dist )
{
  // empty
}

void Heun::step()
{
  // Generate 1D array of Wiener increments
  if( !manualWiener )
    for( auto &i : dw )
      i = sqrt(h) * gen();

  // PREDICTION
  getLE().computeDrift( tmp1, state, getTime() );
  getLE().computeDiffusion( tmp2, state, getTime() );
  for( array_f::index i=0; i!=dim; i++ )
    {
      xPred[i] = state[i] + tmp1[i]*h;
      for ( array_f::index j=0; j!=dim; j++ )
				xPred[i] += tmp2[i][j]*dw[j];
    }

  getLE().computeDrift( tmp1Up, xPred, getTime()+h );
  getLE().computeDiffusion( tmp2Up, xPred, getTime()+h );

  for( int i=0; i<dim; i++ )
    {
      xPred[i] = state[i] + 0.5*h*( tmp1Up[i] + tmp1[i] );
      for( int j=0; j<dim; j++ )
	xPred[i] += 0.5*dw[j]*( tmp2Up[i][j] + tmp2[i][j] );
    }

  setState( xPred );
  setTime( getTime()+h );
}

// Set the manual wiener process mode
void Heun::setManualWienerMode( const bool s ) { manualWiener=s; }
void Heun::setWienerIncrements( const array_f a ) { dw=a; }
