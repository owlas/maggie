// main.cpp
// Main function for testing Maggie code
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#include<iostream>
using std::cout;
using std::endl;

#include<LangevinEquation.hpp>

int main(void)
{
  cout << "Attempting to create StocLLG object..." << endl;
  StocLLG llg(1,2,3,4,5,6);
  return 0;
}
