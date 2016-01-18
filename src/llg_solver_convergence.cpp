// llg_solver_convergence
// Produces convergence data for Heun and Milstein schemes on the llg equation
//
// Oliver W. Laslett
// O.Laslett@soton.ac.uk
//
#include <maggie.hpp>
#include <boost/multi_array.hpp>
using array_d = boost::multi_array<double,1>;
using matrix_d =  boost::multi_array<double,2>;
typedef boost::multi_array_types::index bIndex;
#include <boost/random.hpp>
using mt19937=boost::random::mt19937;
using normal_f=boost::random::normal_distribution<double>;
#include <cmath>
#include <iostream>
using std::cout; using std::endl;

int main() {

  cout << "Running convergence tests for integrators" << endl;

  // Run two numerical solvers for
  int N_dt{ 5 }; // different time steps and
  int N_samples{ 100 }; // different runs

  // Simulation parameters
  const double K{ 1 }, D{ 5e-9 }, V{ 4.0/3.0*M_PI*pow( D/2, 3 ) }
  , KB{ 1.38064852e-23 }, T{ 300 }, alpha{ 0.1 }
  , MU0{ 1.25663706e-6 }, Ms{ 1400e3 }, Hk{ 2*K/( MU0*Ms ) }, hz{ 100e3/Hk }
  , gamma{ 1.7609e11 }, tfactor{ gamma*MU0*Hk/( 1+pow( alpha,2 ) ) };
  const double sigma{ std::sqrt( alpha*KB*T/( K*V*( 1+pow( alpha,2 ) ) ) ) };

  // Set up the LLG
  const StocLLG<double> sde( sigma, alpha, 0.0, 0.0, hz );

  // Initial condition for the system
  array_d init( boost::extents[3] );
  init[0] = 1.0; init[1] = 0.0; init[2] = 0.0;

  // Simulation length
  const double sim_length{ 2e-10 };

  // Store the strong convergence errors
  array_d e_strong_heu( boost::extents[N_dt] );
  for( auto& x : e_strong_heu ) x=0.0;
  array_d e_strong_mil{ e_strong_heu };
  array_d dt_values{ e_strong_mil };

  // Compute a very accurate solution with the Euler scheme
  cout << "computing accurate solution" << endl;
  mt19937 rng_e( 99 );
  double euler_dt{ 1e-19 };
  double euler_dtau{ euler_dt * tfactor};
  int euler_steps{ int( sim_length / euler_dt ) };
  Euler<double> inte_euler( sde, init, 0.0, euler_dtau, rng_e );

  // store the solutions here
  matrix_d sols( boost::extents[3][N_samples] );

  for( int i=0; i!=N_samples; ++i )
  {
      inte_euler.reset();
      for( int n=0; n!=euler_steps; ++n )
          inte_euler.step();

      if( !( i%( N_samples/10 ) ) )
          cout << ".";


      for( array_d::index k=0; k!=3; ++k )
          sols[k][i] = inte_euler.getState()[k];
  }
  cout << endl;

  // Now run the other integrators for various time steps
  for( int i=0; i!=N_dt; i++ )
    {
      // compute the time interval and reduced time interval
      const double dt{ double( pow(10,-12.8)*pow( 10, -0.3*i ) ) };
      const double dtau{ dt*tfactor };

      cout << "dt = " << dt << endl;

      // Create random number generators
      mt19937 rng( 99 );
      mt19937 rng_mil( 188888 ); // needed for milstein
      mt19937 rng_h( 99 );

      // Set up the integrators
      Heun<double> inte_heu( sde, init, 0.0, dtau, rng );
      Milstein<double> inte_mil( sde, init, 0.0, dtau, rng, rng_mil );

      // Store the errors
      double err_heu{ 0.0 }, err_mil{ 0.0 };

      // Run each integrator multiple times
      for( int j=0; j!=N_samples; j++ )
      {
          inte_heu.reset();
          inte_mil.reset();

          // Run the integrators for simulation length and compute
          // the Wiener process manually
          int N_steps{ int( sim_length / dt ) };
          for( int n=0; n!=N_steps; n++ )
          {

              // step the integrators
              inte_heu.step();
              inte_mil.step();
          }

          // Compute a very accurate solution using the Euler scheme
          for( int n=0; n!=euler_steps; ++n )
              inte_euler.step();

          // Add the errors
          for ( bIndex k=0; k!=3; ++k )
          {
              err_heu += std::abs( inte_heu.getState()[k] -
                                   sols[k][i]);
              err_mil += std::abs( inte_mil.getState()[k] -
                                   sols[k][i] );
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
  boostToFile<double>( dt_values, "convergence.dt" );
  boostToFile<double>( e_strong_heu, "convergence.heun" );
  boostToFile<double>( e_strong_mil, "convergence.mils" );
}
