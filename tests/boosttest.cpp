#include <iostream>
using std::cout;
using std::endl;

#include <boost/multi_array.hpp>

typedef boost::multi_array<float,1> array_f;

// class definition
class Test
{
public:
  Test( float, float, float );
  void printarr();
  void setarr( float, float, float );
  array_f& getarr();

private:
  array_f arrmem;
};

// Define some function
float vectorsum( array_f& );

// implement function
float vectorsum( array_f& h )
{
  float sum = 0;
  for ( int i=0; i<3; i++)
    sum += h[i];
  return sum;
}

// constructor
Test::Test( float x, float y, float z )
  : arrmem( boost::extents[3] )
{
  setarr( x, y, z );
}

// Set arr function
void Test::setarr( float x, float y, float z )
{
  arrmem[0] = x;
  arrmem[1] = y;
  arrmem[2] = z;
}

// get arr function
array_f& Test::getarr()
{
  return arrmem;
}

// print function
void Test::printarr()
{
  cout << arrmem[0] << arrmem[1] << arrmem[2] << endl;
}

int main( void )
{
  Test t( 1.0, 2.0, 3.0 );
  cout << "Created test class instance with 1.0 2.0 3.0" << endl;
  t.printarr();

  t.setarr( 6.0, 5.0, 4.0 );
  cout << "Set array to 4.0 5.0 6.0" << endl;
  t.printarr();

  array_f arr( boost::extents[3] );
  cout << "Created boost array" << endl;

  arr[0] = 1.2;
  arr[1] = 2.2;
  arr[2] = 3.2;

  cout << "Sum of vector [1.2,2.2,3.2] is..." << endl;
  cout << vectorsum( arr ) << endl;

  array_f x( boost::extents[2] );
  x[0] = 1;
  x[1] = 0;
  array_f& xref = x;

  array_f y( boost::extents[2] );
  y = xref;
  y[1] = 2;
  cout << "y is [" << y[0] << "," << y[1] << "]" << endl;  
  cout << "x is [" << x[0] << "," << x[1] << "]" << endl;

  // What is the index thing for?
  typedef boost::multi_array<double,3> cube_d;
  cube_d c( boost::extents[2][3][4] );
  
  // Set array to zero
  for( int i=0; i<2; i++ )
    for( int j=0; j<3; j++ )
      for( int k=0; k<4; k++ )
	c[i][j][k] = 0;
  cout << "3d array initialised" << endl;

  // Now access out of bounds
  typedef cube_d::index idx;
  cout << "expect out of bounds error" << endl;
  for( idx i=0; i<2; i++ )
    for( idx j=0; j<4; j++ )
      for( idx k=0; k<4; k++ )
	cout << "index cube at " << i << " " << j << " " << k << endl;
  return 1;
}
