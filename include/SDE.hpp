// SDE.h
// Header for abstract base class
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef SDE_H
#define SDE_H

#include <ODE.hpp>

template <size_t DIM>
using darray = std::array<double,DIM>;
template <size_t DIM1, size_t DIM2>
using dmatrix = std::array<std::array<double,DIM2>,DIM1>;
template <size_t DIM1, size_t DIM2, size_t DIM3>
using darray3 = std::array<std::array<std::array<double,DIM3>,DIM2>,DIM1>;

template <size_t DIM, size_t WDIM>
class SDE : public ODE<DIM>
{
public:
  SDE() {};

  // Useful aliases
  using matrix = dmatrix<DIM,WDIM>;
  using array3 = darray3<DIM,WDIM,WDIM>;
  static const size_t wdim = DIM;

  // Return Langevin components from a given state vector
  virtual void computeDiffusion( matrix& out, const typename SDE<DIM,WDIM>::array& in,
                                 const double t ) const = 0;

  // Langevin equations can specify derivative terms for the diffusxion matrix
  // BB[i][j][k] = PD of B[i][j] w.r.t. state[k]
  //
  virtual void computeDiffusionDerivatives( array3& out,
                                            const typename SDE<DIM,WDIM>::array& in,
                                            const double t ) const;
};

#include <tpp/SDE.cpp>
#endif
