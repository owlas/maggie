// Milstein.cpp 
// Implementation of the Milstein Taylor expansion numerical
// solver for stochastic differential equations. Implemented here for the
// stochastic Landau-Lifshitz-Gilbert equation in its reduced form
// 
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
// 
#include <Milstein.hpp>

// Constructor
Milstein::Milstein( LangevinEquation &le, array_f& init_state,
		    float time, float dt, mt19937& rng_1, mt19937& rng_2 )
  : Integrator( le, init_state, time )
  , h( dt )
  , dim( le.getDim() )
  , dist_1( 0,sqrt( h ) )
  , dist_2( 0,1 )
  , gen_1( rng_1, dist_1 )
  , gen_2( rng_2, dist_2 )
  , next_state( boost::extents[dim] )
  , dw( boost::extents[dim] )
  , dw2( boost::extents[dim] )
  , tmp1( boost::extents[dim] )
  , tmp2( boost::extents[dim][dim] )
  , tmp22( boost::extents[dim][dim] )
  , tmp3( boost::extents[dim][dim][dim] )
  , mu( boost::extents[dim] )
  , eta( boost::extents[dim][p] )
  , zeta( boost::extents[dim][p] )
{
  // empty
}

void Milstein::step()
{
  // Compute the drift and diffusion at the current time step
  getLE().computeDrift( tmp1, state, getTime() );
  getLE().computeDiffusion( tmp2, state, getTime() );
  getLE().computeDiffusionDerivatives( tmp3, state, getTime() );

  // Compute a Wiener process increment (unless manual mode)
  if( !manualWiener )
    for( auto &i : dw )
      i = gen_1();

  // Scaled Wiener process for computation of multiple integrals below
  for( array_f::index i=0; i!=dim; i++ )
    dw2[i] = dw[i]/sqrt( h );

  // rh0 for the calculation below
  float rho=0.0;
  for( unsigned int r=0; r!=p; r++ )
    rho += 1/pow( r, 2 );
  rho *= 1/( 2*pow( M_PI, 2 ) );
  rho += 1.0/12.;

  // Independent Gaussian RVs for integrals
  for( unsigned int i=0; i!=mu.num_elements(); i++ )
    *( mu.data()+i ) = gen_2();
  for( unsigned int i=0; i!=eta.num_elements(); i++ )
    *( eta.data()+i ) = gen_2();
  for( unsigned int i=0; i!=zeta.num_elements(); i++ )
    *( zeta.data()+i ) = gen_2();

  // Compute the multiple Stratonovich integrals
  for( matrix_f::index j1=0; j1!= dim; j1++ )
    for( matrix_f::index j2=0; j2!=dim; j2++ )
      {
	if( j1==j2 )
	  tmp22[j1][j2] = 0.5*pow( dw[j1],2 );
	else
	  {
	    tmp22[j1][j2] = 0;
	    for( array_f::index r=0; r!=p; r++ )
	      tmp22[j1][j2] += ( zeta[j1][r]*( sqrt( 2 )*dw2[j2]
					       + eta[j2][r] )
				 - zeta[j2][r]*( sqrt( 2 )*dw2[j1]
						 + eta[j1][r] ) ) / r;
	    tmp22[j1][j2] *= h/( 2*M_PI );
	    tmp22[j1][j2] += h*( 0.5*dw2[j1]*dw2[j2] + sqrt( rho ) 
				 *( mu[j1]*dw2[j2] - mu[j2]*dw2[j1] ) );
	  }
      }

  // Compute the next state for each vector
  for( array_f::index k=0; k!=dim; k++ )
    {
      next_state[k] = state[k] + tmp1[k]*h;

      for( array_f::index j=0; j!=dim; j++ )
	  next_state[k] += tmp2[k][j]*dw[j];

      for( array_f::index j1=0; j1!=dim; j1++ )
	for( array_f::index j2=0; j2!=dim; j2++ )
	  for( array_f::index j3=0; j3!=dim; j3++ )
	    next_state[k] += tmp2[j3][j1]*tmp3[k][j2][j1]*tmp22[j1][j2];
    }

  setState( next_state );
  setTime( getTime()+h );
}

// Getters and setters for the truncation number
void Milstein::setP( const int pset ) { p=pset; }
int Milstein::getP() const { return p; }

// Set the manual wiener process mode
void Milstein::setManualWienerMode( const bool s ) { manualWiener=s; }
void Milstein::setWienerIncrements( const array_f a ) { dw=a; }
