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
  
  return 1;
}
