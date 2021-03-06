// Utility.hpp - header file for utility functions used by the maggie
// library
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//

#ifndef UTIL_H
#define UTIL_H

#include<iostream>
using std::cout; using std::endl;
using std::flush;
#include<fstream>
#include<boost/multi_array.hpp>
using bidx = boost::multi_array_types::index;
template <typename T> using array = boost::multi_array<T,1>;
template <typename T> using matrix = boost::multi_array<T,2>;
#include <stdexcept>

// saves a boost multiarray to a file with the supplied name
template <typename T>
int boostToFile( array<T> a, std::string name )
{
    // Open the file
    std::ofstream fid;
    fid.open( name );
    for( const auto& i : a )
        fid << i << " ";

    fid.close();
    return 1;
}

template <typename T>
int boostToFile( matrix<T> m, std::string name )
{
    // Open the file
    std::ofstream fid;
    fid.open( name );
    for( const auto &i : m )
    {
        for( const auto& j : i )
            fid << j << " ";
        fid << "\n";
    }

    fid.close();
    return 1;
}

template <class Container>
void containerToFile( const Container &in, const std::string name, const size_t Nwrite )
{
    // set to raise errors for fail read/write
    std::ofstream fid;
    fid.exceptions( std::ofstream::failbit | std::ofstream::badbit );

    // open the file
    fid.open( name );

    // write each element
    for( unsigned int i=0; i!=Nwrite; ++i )
        fid << in[i] << " ";
    fid.close();
}

template <class Container>
void containerToFile( const Container &in, const std::string name )
{
    containerToFile( in, name, in.size() );
}


// Prints a progress bar
void progress_bar( std::string, double );
#endif
