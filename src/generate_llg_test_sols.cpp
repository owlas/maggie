// generate_llg_test_sols.cpp
// Very simply uses the euler scheme to solve the LLG with high accuracy
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
#include <chrono>

int main() {

  cout << "Generating solutions" << endl;

  int N_samples{ 100 }; // different runs

  // Simulation parameters
  const double K{ 1 }, D{ 5e-9 }, V{ 4.0/3.0*M_PI*pow( D/2, 3 ) }
  , KB{ 1.38064852e-23 }, T{ 300 }, alpha{ 0.1 }
  , MU0{ 1.25663706e-6 }, Ms{ 1400e3 }, Hk{ 2*K/( MU0*Ms ) }, hz{ 100e3/Hk }
  , gamma{ 1.7609e11 }, tfactor{ gamma*MU0*Hk/( 1+pow( alpha,2 ) ) };
  const double s{ std::sqrt( alpha*KB*T/( K*V*( 1+pow( alpha,2 ) ) ) ) };
  const double sigma = double( s );

  // Set up the LLG
  const StocLLGIto<double> sde( sigma, alpha, 0.0, 0.0, hz );

  // Initial condition for the system
  array_d init( boost::extents[3] );
  init[0] = 1.0; init[1] = 0.0; init[2] = 0.0;

  // Simulation length
  const double sim_length{ 1e-9 };

  // Compute a very accurate solution with the Euler scheme
  mt19937 rng( 99 );
  double dt{ 1e-15 };
  double dtau{ dt * tfactor};
  int steps{ int( sim_length / dt ) };
  auto inte = Euler<double>( sde, init, 0.0, dtau, rng );

  // store the solutions here
  matrix_d sols( boost::extents[3][N_samples] );

  // start timer
  auto start = std::chrono::steady_clock::now();

  // generate samples
  int N_savepoints{ 1000 };
  int stride{ steps/N_savepoints };
  matrix_d res( boost::extents[3][N_savepoints]);
  for( int i=0; i!=N_samples; ++i )
  {
      inte.reset();
      for( int n=0; n!=steps; ++n )
      {
          inte.step();
          if( !( n%stride ) )
              for( array_d::index k=0; k!=3; ++k )
                  res[k][n/stride] = inte.getState()[k];
      }

      for( array_d::index k=0; k!=3; ++k )
          sols[k][i] = inte.getState()[k];
  }

  // end timer print diff
  auto end = std::chrono::steady_clock::now();

  // Save the results
  boostToFile( sols, "llg_sols.out" );
  boostToFile( res, "test.out");

  // print time taken
  auto diff = end - start;
  cout << std::chrono::duration<double, std::milli> (diff).count() << " ms" << endl;
}
