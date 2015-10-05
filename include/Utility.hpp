// Utility.hpp - header file for utility functions used by the maggie
// library
//
// Oliver W. Laslett (2015)
// O.Laslett@soton.ac.uk
//

#ifndef UTIL_H
#define UTIL_H

#include<iostream>
#include<fstream>
#include<stdexcept>
#include<boost/multi_array.hpp>


// saves a boost array to a file with the supplied name
template<typename T>
int boostToFile( T &array, std::string name );
#endif
