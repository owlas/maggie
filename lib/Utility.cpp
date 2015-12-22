// Utility.cpp - implementation of utility functions
//
// Oliver W. laslett (2015)
// O.Laslett@soton.ac.uk
//
#include<Utility.hpp>

int boostToFile( array_f &array, std::string name )
{
  // Open the file
  std::ofstream fid;
  fid.open( name );

  for( const auto &i : array )
    fid << i << " ";

  fid.close();
  return 1;

}

int boostToFile( matrix_f &array, std::string name )
{
  // Open the file
  std::ofstream fid;
  fid.open( name );

  for( const auto &i : array )
    {
      for( const auto &j : i )
        fid << j << " ";
      fid << "\n";
    }

  fid.close();
  return 1;
}
