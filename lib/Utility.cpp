// Utility.cpp - implementation of utility functions
//
// Oliver W. laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<Utility.hpp>

template<typename T>
int boostToFile( T &array, std::string name )
{
  // Get the dimensions
  int nDims = array.Dimensionality();

  // Open the file
  std::ofstream fid;
  fid.open( name );

  if( nDims == 1 )
    {
      for( int i=0; i<array.shape()[0]; i++ )
        fid << array[i] << " ";
    }
  else if( nDims == 2 )
    {
      for( int i=0; i<array()[0]; i++ )
        {
          for( int j=0; j<array()[1]; j++ )
            fid << array[i][j] << " ";
          fid << "\n";
        }
    }
  else
    throw std::invalid_argument( "Can only save boost arrays with less"
                                 " than two dimensions." );

  fid.close();
  return 1;

}
