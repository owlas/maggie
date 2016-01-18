// StocLLG.h
// Header for the stochastic LLG Langevin equation
//
// Oliver W. Laslett 2015
// O.Laslett@soton.ac.uk
//
#ifndef STOCLLG_H
#define STOCLLG_H

#include <SDE.hpp>

#include <boost/multi_array.hpp>
template <typename T> using array = boost::multi_array<T,1>;
template <typename T> using matrix = boost::multi_array<T,2>;
template <typename T> using array3 = boost::multi_array<T,3>;
using bidx = boost::multi_array_types::index;

template <typename T>
class StocLLG : public SDE<T>
{
 public:
        StocLLG( const T s, const T a,
		 const T hx, const T hy, const T hz );

  // returns the drift vector
    virtual void computeDrift( array<T>&, const array<T>&, const T ) const;
  // returns the diffusion marix
    virtual void computeDiffusion( matrix<T>&, const array<T>&,
				 const T ) const;
  // returns a vector of derivative sums for Taylor
    virtual void computeDiffusionDerivatives( array3<T> &out, const array<T> &in,
					    const T ) const;

  void setReducedHeff( const T, const T, const T );
    array<T> getReducedHeff() const;

  void setSigma( const T );
  T getSigma() const;

  void setAlpha( const T );
  T getAlpha() const;

 private:
    array<T> h;
    float sigma;
    float alpha;
};

// This is the Ito version of the Stochastic LLG
template <typename T>
class StocLLGIto : public StocLLG<T>
{
public:
    StocLLGIto( const T s, const T a,
                const T hx, const T hy, const T hz );

    // returns the drift vector
    void computeDrift( array<T>&, const array<T>&, const T ) const;

};

#include <tpp/StocLLG.cpp>
#endif
