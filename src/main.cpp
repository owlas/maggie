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
using std::abs;
using std::cos;
using std::sin;

#include<maggie.hpp>
#include<gtest/gtest.h>
#include<boost/multi_array.hpp>
typedef boost::multi_array<float,1> array_f;
typedef boost::multi_array<float,2> matrix_f;

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
  for( int i=0; i<1000; i++ )
    inte.step();

  // Get the state and time
  state = inte.getState();
  t = inte.getTime();

  // Check the solution
  EXPECT_LE( abs(  cos( t ) - state[0]), 1e-6 );
  EXPECT_LE( abs( -sin( t ) - state[1]), 1e-6 );
}

  

int main( int argc, char **argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
