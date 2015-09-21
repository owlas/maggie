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
using std::cos;
using std::sin;

#include<maggie.hpp>
#include<gtest/gtest.h>
#include<gmock/gmock.h>
#include<boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;
typedef boost::multi_array<float,2> matrix_f;
typedef boost::multi_array<float,3> cube_f;

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
  llg.computeDrift( out, in );

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
  llg.computeDiffusion( out, in );
  
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
  state[0] = 1;
  state[1] = 0;

  // Basic differential equation
  class ode : public LangevinEquation
  {
  public:
    ode() : LangevinEquation( 2 ) {}; // constructor

    // Differential equation
    virtual void computeDrift( array_f& out, array_f& in)
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
  
  

int main( int argc, char **argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
