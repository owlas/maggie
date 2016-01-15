// llg_solver_convergence
// Produces convergence data for Heun and Milstein schemes on the llg equation
// 
// Oliver W. Laslett
// O.Laslett@soton.ac.uk
// 
#include <maggie.hpp>
#include <boost/multi_array.hpp>
using array_f = boost::multi_array<float,1>;
using matrix_f =  boost::multi_array<float,2>;
typedef boost::multi_array_types::index bIndex;
#include <boost/random.hpp>
using mt19937=boost::random::mt19937;
using normal_f=boost::random::normal_distribution<float>;
#include <cmath>
#include <iostream>
using std::cout; using std::endl;

int main() {

  cout << "Running convergence tests for integrators" << endl;

  // Run two numerical solvers for
  int N_dt{ 5 }; // different time steps and
  int N_samples{ 100 }; // different runs

  // Simulation parameters
  const float K{ 1 }, D{ 5e-9 }, V{ 4.0/3.0*M_PI*pow( D/2, 3 ) }
  , KB{ 1.38064852e-23 }, T{ 300 }, alpha{ 0.1 }
  , MU0{ 1.25663706e-6 }, Ms{ 1400e3 }, Hk{ 2*K/( MU0*Ms ) }, hz{ 100e3/Hk }
  , gamma{ 1.7609e11 }, tfactor{ gamma*MU0*Hk/( 1+pow( alpha,2 ) ) };
  const double s{ std::sqrt( alpha*KB*T/( K*V*( 1+pow( alpha,2 ) ) ) ) };
  const float sigma = float( s );

  // Set up the LLG
  const StocLLG sde( sigma, alpha, 0.0, 0.0, hz );

  // Simulation length
  const float sim_length{ 2e-10 };

  // Store the strong convergence errors
  array_f e_strong_heu( boost::extents[N_dt] );
  for( auto& x : e_strong_heu ) x=0.0;
  array_f e_strong_mil{ e_strong_heu };
  array_f dt_values{ e_strong_mil };
  
  for( int i=0; i!=N_dt; i++ )
    {
      // compute the time interval and reduced time interval
      const float dt{ float( pow(10,-12.8)*pow( 10, -0.3*i ) ) };
      const float dtau{ dt*tfactor };

      cout << "dt = " << dt << endl;

      // Create random number generators
      mt19937 rng( 99 );
      mt19937 rng_mil( 188888 ); // needed for milstein

      // Variate generator for the Wiener process
      normal_f dist( 0, 1 );
      boost::random::variate_generator< mt19937, normal_f >
	vg( rng, dist );

      // Initial condition for the system
      array_f init( boost::extents[3] );
      init[0] = 1.0; init[1] = 0.0; init[2] = 0.0;

      // Set up the integrators
      Heun inte_heu( sde, init, 0.0, dtau, rng );
      Milstein inte_mil( sde, init, 0.0, dtau, rng, rng_mil );
 
      // Set up manual mode for the integrators
      inte_heu.setManualWienerMode( true );
      inte_mil.setManualWienerMode( true );
      array_f dw( boost::extents[3] );
      for( auto& i : dw )
	i = 0.0;

      // Store the errors
      float err_heu{ 0.0 }, err_mil{ 0.0 };

      // Run each integrator multiple times
      for( int j=0; j!=N_samples; j++ )
	{
	  inte_heu.reset();
	  inte_mil.reset();
	  
	  // Run the integrators for simulation length and compute
	  // the Wiener process manually
	  int N_steps{ int( sim_length / dt ) };
	  float W{ 0 };
	  for( int n=0; n!=N_steps; n++ )
	    {
	      // get Gaussian RV and scale for both solutions
	      float val = vg();
	      dw[2] = val*std::sqrt( dtau ); // reduced time
	      W += val*std::sqrt( dt ); // real time

	      // save the wiener increment vector to the integrators
	      inte_heu.setWienerIncrements( dw );
	      inte_mil.setWienerIncrements( dw );
	      
	      // step the integrators
	      inte_heu.step();
	      inte_mil.step();
	    }
	  
	  // Analytic solution from Hannay
	  array_f anSol( boost::extents[3] );
	  float kb{ 1.38064852e-16 }, gamma{ 1.7609e7 }, alpha{ 0.1 }
	  , H{ 400*M_PI }, V{ 4.0/3.0*M_PI*pow( 2.5e-7,3 ) }, T{ 300 }
	  , Ms{ 1400 }, sigma{ std::sqrt( 2*alpha*kb*T/( gamma*Ms*V ) ) }
	  , t{ dt*N_steps };
	  anSol[0] =
	    1/std::cosh( alpha*gamma*( H*t+sigma*W ) / ( 1+pow( alpha, 2 ) ) )
	    * std::cos( gamma*( H*t+sigma*W )/( 1+pow( alpha,2 ) ) );

	  anSol[1] =
	    1/std::cosh( alpha*gamma*( H*t+sigma*W ) / ( 1+pow( alpha,2 ) ) )
	    * std::sin( gamma*( H*t+sigma*W )/( 1+pow( alpha,2 ) ) );

	  anSol[2] =
	    std::tanh( alpha*gamma*( H*t+sigma*W ) / ( 1+pow( alpha,2 ) ) );

	  // Add the errors
	  for ( bIndex i=0; i!=3; ++i )
	    {
	      err_heu += std::abs( inte_heu.getState()[i] - anSol[i] );
	      err_mil += std::abs( inte_mil.getState()[i] - anSol[i] );
	    }
	  // DEBUG
	  cout << "heun: " << inte_heu.getState()[2] << endl;

	} // end for multiple samples

      // Compute the expected absolute error
      e_strong_heu[i] += err_heu/N_samples;
      e_strong_mil[i] += err_mil/N_samples;

      // store time step
      dt_values[i] = dt;

    } // end for each dt

  // Save the convergence values
  boostToFile( dt_values, "convergence.dt" );
  boostToFile( e_strong_heu, "convergence.heun" );
  boostToFile( e_strong_mil, "convergence.mils" );
}
