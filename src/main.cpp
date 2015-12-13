// main.cpp
// Main function for testing Maggie code
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<iostream>
using std::cout;
using std::endl;

#include<cmath>
using std::cos; using std::sin;
using std::sqrt; using std::exp;

#include<maggie.hpp>
#include<gtest/gtest.h>
#include<gmock/gmock.h>
#include<boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;
typedef boost::multi_array<float,2> matrix_f;
typedef boost::multi_array<float,3> cube_f;
#include<boost/random.hpp>
using mt19937=boost::random::mt19937;
using normal_f=boost::random::normal_distribution<float>;

#include<vector>

// Test construction of StocLLG
TEST(StochasticLLG, ContructAndDim)
{
  StocLLG llg( 1, 2, 3, 4, 5 );
  EXPECT_EQ( 3, llg.getDim() );
}

// Test deterministic part of StocLLG
TEST(StochasticLLG, Drift)
{
  // declare in and out arrays
  array_f out( boost::extents[3] );
  array_f in( boost::extents[3] );
  in[0] = 2;
  in[1] = 3;
  in[2] = 4;

  // Set up StocLLG and compute drift
  StocLLG llg( 2, 3, 1, 2, 3 );
  llg.computeDrift( out, in, 0 );

  // Are the solutions correct?
  EXPECT_EQ( -34, out[0] );
  EXPECT_EQ( -4, out[1] );
  EXPECT_EQ( 20, out[2] );
}

// Test the stochastic part of the StocLLG
TEST(StochasticLLG, Diffusion)
{
  // declare in and out arrays
  array_f in( boost::extents[3] );
  matrix_f out( boost::extents[3][3] );
  in[0] = 2;
  in[1] = 3;
  in[2] = 4;

  // set up the StocLLG and compute diffusion
  StocLLG llg( 2, 3, 1, 2, 3 );
  llg.computeDiffusion( out, in, 0 );
  
  // are the solutions correct?
  EXPECT_EQ( 150, out[0][0] );
  EXPECT_EQ( -28, out[0][1] );
  EXPECT_EQ( -54, out[0][2] );
  EXPECT_EQ( -44, out[1][0] );
  EXPECT_EQ( 120, out[1][1] );
  EXPECT_EQ( -68, out[1][2] );
  EXPECT_EQ( -42, out[2][0] );
  EXPECT_EQ( -76, out[2][1] );
  EXPECT_EQ( 78, out[2][2] );
}

// Test the RK4 algorithm
TEST(RK4, BasicCircle)
{
  // declare array for initial state
  array_f state( boost::extents[2] );

  float t=0; // initial time

  // fill the array state
  state[0] = 1.0;
  state[1] = 0.0;

  // Basic differential equation
  class ode : public LangevinEquation
  {
  public:
    ode() : LangevinEquation( 2 ) {}; // constructor

    // Differential equation
    virtual void computeDrift( array_f& out, array_f& in, float )
    {
      out[0] = in[1];
      out[1] = -in[0];
    }
  } testOde; 
  

  // Create an instance of the RK4 integrator
  RK4 inte( testOde, state, t, 0.000001 );

  // Run the integrator for 2000 steps
  for( int i=0; i<500; i++ )
    inte.step();

  // Get the state and time
  state = inte.getState();
  t = inte.getTime();

  // Check the solution
  ASSERT_FLOAT_EQ( cos( t ), state[0] );
  ASSERT_FLOAT_EQ( -sin( t ), state[1] );
}

// Test for the Particles
TEST(Particle, SetSize)
{
  array_f uea( boost::extents[3] );
  uea[0] = 1;
  uea[1] = 0;
  uea[2] = 0;
  Particle p( 1, 2, 3, 4, 5, uea );
  ASSERT_FLOAT_EQ( 33.51032, p.getV() );
}

// Check that uea is a unit vector
TEST(Particle, UnitVector)
{
  array_f uea( boost::extents[3] );
  // Not a unit vector
  uea[0] = 0.3;
  uea[1] = 0.5;
  uea[2] = 0.7;
  ASSERT_THROW( Particle p( 1, 2, 3, 4, 5, uea ),
		std::invalid_argument );
  // is a unit vector
  uea[0] = 0.5773503;
  uea[1] = 0.5773503;
  uea[2] = 0.5773503;
  ASSERT_NO_THROW( Particle p( 1, 2, 3, 4, 5, uea ) );
}

// Check computation of the energy barriers
TEST(Particle, EnergyBarriers)
{
  array_f uea( boost::extents[3] );
  uea[0] = 0;
  uea[0] = 0;
  uea[0] = 1;
  Particle p( Constants::GYROMAG, 0.1, 3.5e5, 4e-9, 1e4, uea );

  // Cannot have reduced fields > 1 (i.e. H>H_k)
  ASSERT_THROW( p.alignedEnergyBarriers( 1.5 ),
                std::invalid_argument );

  // Check computation
  array_f ans( boost::extents[2] );
  ans = p.alignedEnergyBarriers( 0.2 );

  ASSERT_FLOAT_EQ( 2.14466058e-22, ans[0] );
  ASSERT_FLOAT_EQ( 4.82548631e-22, ans[1] );
}

// Check the computaiton of the transition rates
TEST(Particle, TransitionRates)
{
  // test here
}

// Test the particle cluster, are the distances correct?
TEST(ParticleCluster, Distances)
{
  // create list of 3 identical particles
  array_f uea( boost::extents[3] );
  uea[0] = 1.0;
  uea[1] = 0.0;
  uea[2] = 0.0;  
  Particle p1( 1,2,3,4,5,uea );
  std::vector<Particle> plist = {p1, p1, p1};

  // Set the particle locations
  matrix_f locs( boost::extents[3][3] );
  
  // particle 1 at location (0,0,0)
  for( int i=0; i<3; i++ )
    locs[0][i] = 0;
  // particle 2 at location (1,0,0)
  locs[1][0] = 1;
  locs[1][1] = 0;
  locs[1][2] = 0;
  // particle 3 at location (0,1,1)
  locs[2][0] = 0;
  locs[2][1] = 1;
  locs[2][2] = 1;

  // Create the cluster and get the distances
  ParticleCluster const pclust( plist, locs );
  cube_f dists( boost::extents[3][3][3] );
  dists = pclust.getDistances();
  
  // assert the distances are correct
  // p1 - p1
  EXPECT_EQ( 0, dists[0][0][0] );
  EXPECT_EQ( 0, dists[0][0][1] );
  EXPECT_EQ( 0, dists[0][0][2] );
  // p1 - p2
  EXPECT_EQ( 1, dists[1][0][0] );
  EXPECT_EQ( 0, dists[1][0][1] );
  EXPECT_EQ( 0, dists[1][0][2] );
  // p1 - p3
  EXPECT_EQ( 0, dists[2][0][0] );
  EXPECT_EQ( 1, dists[2][0][1] );
  EXPECT_EQ( 1, dists[2][0][2] );
  // p2 - p1
  EXPECT_EQ( 1, dists[0][1][0] );
  EXPECT_EQ( 0, dists[0][1][1] );
  EXPECT_EQ( 0, dists[0][1][2] );
  // p2 - p2
  EXPECT_EQ( 0, dists[1][1][0] );
  EXPECT_EQ( 0, dists[1][1][1] );
  EXPECT_EQ( 0, dists[1][1][2] );
  // p2 - p3
  EXPECT_EQ( 1, dists[2][1][0] );
  EXPECT_EQ( 1, dists[2][1][1] );
  EXPECT_EQ( 1, dists[2][1][2] );
  // p3 - p1
  EXPECT_EQ( 0, dists[0][2][0] );
  EXPECT_EQ( 1, dists[0][2][1] );
  EXPECT_EQ( 1, dists[0][2][2] );
  // p3 - p2
  EXPECT_EQ( 1, dists[1][2][0] );
  EXPECT_EQ( 1, dists[1][2][1] );
  EXPECT_EQ( 1, dists[1][2][2] );
  // p3 - p3
  EXPECT_EQ( 0, dists[2][2][0] );
  EXPECT_EQ( 0, dists[2][2][1] );
  EXPECT_EQ( 0, dists[2][2][2] );
}

// Test the stochastic integrator on a simple Wiener process
TEST( Heun, HeunWiener )
{
  // Time step
  const float dt = 1e-10;

  // Langevin equation representing a Wiener process
  class Wiener : public LangevinEquation
  {
  public:
    Wiener() : LangevinEquation( 1 ) {};
    virtual void computeDrift( array_f& out, array_f&, float )
    {
      out[0] = 0;
    }
    virtual void computeDiffusion( matrix_f& out, array_f&, float )
    {
      out[0][0] = 1;
    }
  } testSDE;

  // Create random number generators with same seed
  mt19937 rng(1234);
  mt19937 rng_an(1234);

  // Create a variate generator for the analytic solution
  normal_f dist(0, sqrt( dt ) );
  boost::random::variate_generator<mt19937, normal_f>
    vg( rng_an, dist );

  // Set up numerical integration
  array_f init( boost::extents[1] ); init[0]=0;
  auto inte = Heun( testSDE, init, 0.0, dt, rng );

  // Set the initial condition for the analytic solution
  float analyticSol = init[0];

  // Integrate for 5000 steps and compare against analytical
  for( array_f::index i=0; i!=5000; i++ )
  {
    inte.step();
    float numericalSol = inte.getState()[0];
    analyticSol += vg();

    ASSERT_LE( std::abs( analyticSol - numericalSol ), 1e-5 )
                  << "Steps completed =" << i;
  }
}


// Test the integration scheme
TEST( Heun, HeunOrnsteinUhlenbeck )
{
  // Time step
  float dt = 1e-10;

  // Class for the Ornstein-Uhlenbeck process
  class OH : public LangevinEquation
  {
  public:
    OH() : LangevinEquation( 1 ) {};
    virtual void computeDrift( array_f& out, array_f& in, float )
    {
      out[0] = theta*( mu-in[0] );
    }
    virtual void computeDiffusion( matrix_f& out, array_f&, float )
    {
      out[0][0] = sigma;
    }

    // Parameters for OU process
    float theta=10, mu=-1.0, sigma=0.8;
  } testSDE;

  // create the RNG
  mt19937 rng(555);
  mt19937 rng_an(555);

  // Create a variate generator for the analytic solution
  normal_f dist(0, 1 );
  boost::random::variate_generator<mt19937, normal_f>
    vg( rng_an, dist );

  // Create the integrator
  array_f init( boost::extents[1] ); init[0]=-3;
  auto inte = Heun( testSDE, init, 0.0, dt, rng );

  // Track the solutions
  float numericalSol = init[0];
  float analyticSol = init[0];

  float theta = testSDE.theta;
  float mu = testSDE.mu;
  float sigma = testSDE.sigma;

  for( int i=0; i<300; i++ )
    {
      // Numerical solution
      inte.step();
      numericalSol = inte.getState()[0];

      // Analytic solution
      // FROM http://henley.ac.uk/web/FILES/REP/Monte_Carlo_Simulation_of_Stochastic_Processes.pdf

      analyticSol = analyticSol*exp( -theta*dt ) + mu*( 1-exp( -theta*dt ) )
	+ sigma*sqrt( ( 1-exp( -2*theta*dt ) )/( 2*theta ) ) * vg();

      ASSERT_LE( std::abs( analyticSol - numericalSol ), 1e-4 )
        << "Steps completed: " << i << std::endl;
    }
  }

// Test run of the Heun and LLG algorithms
TEST( StocLLG, HeunIntegrationSolution )
{
  // write test for heun integration here.
}

// Test the quadrature rule
// First test is simple, and checks against a straight line
// integration.
TEST( Quadrature, LinearFunction )
{
  // integrate between 0 and 7 with default tolerance
  // on a simple linear function
  Quad::Quadrature q( []( float x ){ return 2*x+3; }, 0.0, 7.0 );
  float res = q.qTrap();

  // check that the solution is correct
  ASSERT_FLOAT_EQ( 70.0, res );
}

// Does it work for an odd function?
TEST( Quadrature, OddFunctionOver0 )
{
  // integrate a linear function odd function
  // from -2 to 2 to  get 0 area
  Quad::Quadrature q( []( float x){ return 1.2*x; }, -2.0, 2.0 );
  float res = q.qTrap();

  ASSERT_FLOAT_EQ( 0.0, res );
}

// Test the vector function
TEST( Quadrature, LinearFunctionFromVector )
{
  // generate a vector of data
  const int N=10000;
  const float h = 7.0/N;
  array_f vec( boost::extents[N] );
  for( int i=0; i<N; i++ ) vec[i] = 2*i*h+3;
  
  // get the result
  float res = Quad::trapVec( vec, h );
  
  ASSERT_FLOAT_EQ( 70.0, res );
}


// Test the two state master equation
// Is the sum of the probabilities always unity?
TEST( TwoStateMasterEquation, ConservationOfProbability )
{
  using func1DFloat=std::function<float( float ) >;

  cout.precision(10);

  // two simple functions for the transition rates
  func1DFloat w1 = []( float t ){ return t*0.01 + 0.3; };
  func1DFloat w2 = []( float t ){ return std::pow( t,2 )*0.1 + 0.5; };

  // create the two state master eqution
  TwoStateMasterEquation eq( w1, w2 );
  
  // Integrate from an initial condition
  array_f init( boost::extents[2] );
  init[0] = 0.5;
  init[1] = 0.5;

  RK4 integrator( eq, init, 0, 1e-7 );

  // step the integrator and check that the probability is conserved
  // throughout.
  array_f state( boost::extents[2] );

  state = integrator.getState();
  for( int i=0; i<200; i++ )
    {
      integrator.step();
      state = integrator.getState();
      ASSERT_LE( std::abs( 1.0-state[0]-state[1] ), 1e-5 );
    }
}


// Write test here for input output behaviour
TEST( IOTests, ReadAndWrite )
{
        
}

// Test the results of the Kramers calculations
TEST( Kramers, MaxLocation )
{
  float h=0.3, psi=0.7;  
  ASSERT_LE( std::abs( 1.862366031 -
		       KramersTrig::theta_max( h, psi ) )
	     , pow( h,4 ) );
}

TEST( Kramers, MinLocationOne )
{
  float h=0.3, psi=0.7;  
  ASSERT_LE( std::abs( 0.157479602 -
		       KramersTrig::theta_min1( h, psi ) )
	     , pow( h,4 ) );
}

TEST( Kramers, MinLocationTwo )
{
  float h=0.3, psi=0.7;
  ASSERT_LE( std::abs( 2.885440350 -
		       KramersTrig::theta_min2( h, psi ) )
	     , pow( h,4 ) );
}

// Check the energy barrier calculations in reduced energy
// i.e. e=E/(2VK)
TEST( Kramers, EBarOne )
{
  float h=0.3, psi=0.7, s=3.0;
  float betaV = KramersTrig::k_ebar_1( s, h, psi );
  float res = betaV/s/2.0;
  ASSERT_LE( std::abs( 0.584159039 - res )
	     , pow( h,4 ) );
}
TEST( Kramers, EBarTwo )
{
  float h=0.3, psi=0.7, s=3.0;
  float betaV = KramersTrig::k_ebar_2( s, h, psi );
  float res = betaV/s/2.0;
  ASSERT_LE( std::abs( 0.134437709 - res )
	     , pow( h,4 ) );
}
  

int main( int argc, char **argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
