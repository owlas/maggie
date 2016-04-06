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
using range = boost::multi_array_types::index_range;
#include <boost/random.hpp>
using mt19937=boost::random::mt19937;
using normal_d=boost::random::normal_distribution<double>;
#include <cmath>
#include <iostream>
using std::cout; using std::endl;
#include <string>

int main() {

  cout << "Running convergence tests for integrators" << endl;

  // Run two numerical solvers for
  int N_dt{ 6 }; // different time steps and
  int N_samples{ 1000 }; // different runs

  // Simulation parameters
  const double K{ 1 }, D{ 5e-9 }, V{ 4.0/3.0*M_PI*pow( D/2, 3 ) }
  , KB{ 1.38064852e-23 }, T{ 300 }, alpha{ 0.1 }
  , MU0{ 1.25663706e-6 }, Ms{ 1400e3 }, Hk{ 2*K/( MU0*Ms ) }, hz{ 100e3/Hk }
  , gamma{ 1.7609e11 }, tfactor{ gamma*MU0*Hk/( 1+pow( alpha,2 ) ) };
  const double sigma{ std::sqrt( alpha*KB*T/( K*V*( 1+pow( alpha,2 ) ) ) ) };

  // Set up the LLG
  const StocLLG<double> sde( sigma, alpha, 0.0, 0.0, hz );
  // and an ito version
  const StocLLGIto<double> sde_ito( sigma, alpha, 0.0, 0.0, hz );

  // Initial condition for the system
  array_d init( boost::extents[3] );
  init[0] = 1.0; init[1] = 0.0; init[2] = 0.0;

  // Simulation length
  const double sim_length{ 2e-10 };

  // most accurate simulation step
  const double ref_dt{ 1e-14 };
  const int ref_steps{ int( sim_length / ref_dt ) };

  // define the rng
  mt19937 rng( 99 );
  normal_d dist( 0,sqrt( ref_dt ) );
  boost::variate_generator<mt19937, normal_d> vg( rng, dist );

  // use this to store 3d wiener vectors
  array_d ww( boost::extents[3] );

  // store the mz values for each run and time step
  matrix_d sols_heun( boost::extents[N_samples][N_dt] );
  matrix_d sols_eule( boost::extents[N_samples][N_dt] );


  // Now run the solver for Nsamples different wiener paths
  for( int i=0; i!=N_samples; i++ )
    {

        // compute a wiener path
        matrix_d w( boost::extents[3][ref_steps] );
        for( bIndex n=0; n!=3; ++n )
            w[n][0] = 0.0;
        for( bIndex n=1; n!=ref_steps; ++n )
        {
            w[0][n] = w[0][n-1] + vg();
            w[1][n] = w[1][n-1] + vg();
            w[2][n] = w[2][n-1] + vg();
        }

        // some debugging
        cout << "sample " << i << " of " << N_samples << endl;

      // Run an integrator with different step sizes
      for( int j=0; j!=N_dt; j++ )
      {

          // Factor of step size multiplier
          int N_mul = pow( 2, j );

          // compute a solution with the Heun scheme
          const int steps{ ref_steps / N_mul };
          const double dt{ ref_dt * N_mul };
          const double dtau{ dt*tfactor };

          // construct numerical solvers
          // and set the wiener increments manually
          auto heun = Heun<double>( sde, init, 0.0, dtau, rng );
          heun.setManualWienerMode( true );
          auto euler = Euler<double>( sde_ito, init, 0.0, dtau, rng );
          euler.setManualWienerMode( true );

          // step the integrator
          for( int n=0; n!=steps-1; n++ )
          {
              // compute the wiener step
              int nthis = n*N_mul;
              int nnext = nthis + N_mul;
              ww[0] = w[0][nnext] - w[0][nthis];
              ww[1] = w[1][nnext] - w[1][nthis];
              ww[2] = w[2][nnext] - w[2][nthis];

              // set the increments and step
              heun.setWienerIncrements( ww );
              euler.setWienerIncrements( ww );
              heun.step();
              euler.step();
          }

          // store the final mz state
          sols_heun[i][j] = heun.getState()[0];
          sols_eule[i][j] = euler.getState()[0];

      } // end for each time step

    } // end for each sample

  // save the data
  boostToFile(sols_heun, "llg_heun.mz" );
  boostToFile(sols_eule, "llg_euler.mz" );

} // end main
