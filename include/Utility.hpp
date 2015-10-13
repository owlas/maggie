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
#include<boost/multi_array.hpp>


// saves a boost array to a file with the supplied name
int boostToFile( boost::multi_array<float,1> &array, std::string name );
int boostToFile( boost::multi_array<float,2> &array, std::string name );
#endif
