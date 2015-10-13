// Utility.cpp - implementation of utility functions
//
// Oliver W. laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<Utility.hpp>

int boostToFile( boost::multi_array<float,1> &array, std::string name )
{
  // Open the file
  std::ofstream fid;
  fid.open( name );

  for( int i=0; i<array.shape()[0]; i++ )
    fid << array[i] << " ";

  fid.close();
  return 1;

}

int boostToFile( boost::multi_array<float,2> &array, std::string name )
{
  // Open the file
  std::ofstream fid;
  fid.open( name );

  for( int unsigned i=0; i<array.shape()[0]; i++ )
    {
      for( int unsigned j=0; j<array.shape()[1]; j++ )
        fid << array[i][j] << " ";
      fid << "\n";
    }

  fid.close();
  return 1;
}
