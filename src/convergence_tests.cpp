// convergence_tests.cpp
// Produces convergence data for Heun and Milstein schemes
//
// Oliver W. Laslett
// O.Laslett@soton.ac.uk
//
#include <maggie.hpp>
#include <boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;
typedef boost::multi_array<float,2> matrix_f;
#include <boost/random.hpp>
using mt19937=boost::random::mt19937;
using normal_f=boost::random::normal_distribution<float>;
#include <cmath>

int main() {

  // Run two numerical solvers for
  const int N_dt{ 5 }; // different time steps and
  const int N_samples{ 2000 }; // different runs

  const float sim_length{ 1.0 }; // time of simualtion

  // Stochastic differential equation is scalar constant drift
  const float a{ 1.5 };
  const float b{ 0.1 };
  SDE_AXpBX sde( a, b );

  // Store the strong convergence errors
  array_f e_strong_heu( boost::extents[N_dt] );
  for( auto& x : e_strong_heu ) x=0.0;
  array_f e_strong_mil{ e_strong_heu };
  array_f dt_values{ e_strong_mil };

  for( int i=0; i!=N_dt; i++ )
    {
      // compute the next time step
      const float dt{ float( pow(2,-3)*pow( 2, -1*i ) ) };

      // Create random number generators
      mt19937 rng_heu( 99 );
      mt19937 rng_mil( 99 );
      mt19937 rng_ana( 99 );
      mt19937 rng_mil2( 188888 ); // needed for milstein

      // Variate generator for the analytical solution
      normal_f dist( 0, std::sqrt( dt ) );
      boost::random::variate_generator< mt19937, normal_f >
	vg( rng_ana, dist );

      // Initial condition for the system
      array_f init( boost::extents[1] ); init[0] = 1.0;

      // Set up the integrators
      Heun<float> inte_heu( sde, init, 0.0, dt, rng_heu );
      Milstein<float> inte_mil( sde, init, 0.0, dt, rng_mil, rng_mil2 );

      // Store the errors
      matrix_f errors( boost::extents[2][N_samples] );
      // Run each integrator multiple times
      for( int j=0; j!=N_samples; j++ )
	{
	  inte_heu.reset();
	  inte_mil.reset();

	  // Run the integrators for simulation length
	  const int N_steps{ int( sim_length / dt ) };
	  for( int n=0; n!=N_steps; n++ )
	    {
	      inte_heu.step();
	      inte_mil.step();
	    }

	  // compute the analytic solution
	  // see: Kloen "Numerical solution of SDEs"
	  float W{ 0 };
	  for( int n=0; n!=N_steps; n++ )
	    W += vg();
	  float analytic =
	    init[0]*std::exp( (a-0.5*std::pow( b,2 ) )*sim_length + b*W );

	  // store the absolute errors
	  errors[0][j] = std::abs( inte_heu.getState()[0] - analytic );
	  errors[1][j] = std::abs( inte_mil.getState()[0] - analytic );

	  //
	}

      // Compute the expected absolute error
      for( auto& x : errors[0] )
	e_strong_heu[i] += x;
      e_strong_heu[i] /= N_samples;
      for( auto& x : errors[1] )
	e_strong_mil[i] += x;
      e_strong_mil[i] /= N_samples;

      // store time step
      dt_values[i] = dt;
    }

  // Save the convergence values
  boostToFile( dt_values, "convergence.dt" );
  boostToFile( e_strong_heu, "convergence.heun" );
  boostToFile( e_strong_mil, "convergence.mils" );
}
