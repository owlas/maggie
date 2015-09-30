// Heun.cpp
// Implementation for the Heun integration scheme
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<Heun.hpp>


// Constructor
Heun::Heun( LangevinEquation &le, array_f& init_state,
	    float time, float dt, mt19937& rng )
  : Integrator( le, init_state, time)
  , h( dt )
  , dim( le.getDim() )
  , dw( boost::extents[dim])
  , xPred( boost::extents[dim])
  , tmp1( boost::extents[dim])
  , tmp1Up( boost::extents[dim])
  , wienerProducts( boost::extents[dim])
  , tmp2( boost::extents[dim][dim] )
  , tmp2Up( boost::extents[dim][dim] )
  , dist(0,1)
  , gen( rng, dist )
{
  // empty
}

void Heun::step()
{
  for( int i=0; i<dim; i++ )
    dw[i] = sqrt(h) * gen();
  
  getLE().computeDrift( tmp1, state, getTime() );
  getLE().computeDiffusion( tmp2, state, getTime() );
  
  for( int i=0; i<dim; i++ )
    {
      wienerProducts[i] = 0;
      for( int j=0; j<dim; j++ )
	wienerProducts[i] += tmp2[i][j]*dw[j];
      xPred[i] = state[i] + tmp1[i]*h + wienerProducts[i];
    }

  getLE().computeDrift( tmp1Up, xPred, getTime() );
  getLE().computeDiffusion( tmp2Up, xPred, getTime() );
  
  for( int i=0; i<dim; i++ )
    {
      for( int j=0; j<dim; j++ )
	wienerProducts[i] += tmp2Up[i][j]*dw[j];
      wienerProducts[i] *= 0.5;
      xPred[i] = state[i] + 0.5*h*( tmp1[i] + tmp1Up[i] )
	+ wienerProducts[i];
    }
  
  setState( xPred );
  setTime( getTime()+h );
}
