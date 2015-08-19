// main.cpp
// Main function for testing Maggie code
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<iostream>
using std::cout;
using std::endl;

#include<LangevinEquation.h>

int main(void)
{
  cout << "Attempting to create LangevinEquation object..." << endl;
  LangevinEquation myLE(2);
  return 0;
}
