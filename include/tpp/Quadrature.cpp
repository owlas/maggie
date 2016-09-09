// Quadrature.cpp - implementation of qudrature
//
// Based on numerical recipes version
//
// Oliver W. Laslett
// O.Laslett@soton.ac.uk
//

namespace Quad
{
    template <typename T>
    Quadrature<T>::Quadrature( std::function<T( T )> func,
                               const T start,
                               const T end )
    : f( func )
    , a( start )
    , b( end )
  {
    s = 0;
    n = 0;
  }

    template <typename T>
    int Quadrature<T>::getN() const { return n; }

    template <typename T>
    T Quadrature<T>::next()
  {
    T x, tnm, sum, del;
    int it,j;

    // increment the current step
    n++;
    if( n==1 )
      return ( s=0.5*( b-a )*( f( a )+f( b ) ) );
    else
      {
        // compute spacing for the current step
        for( it=1, j=1; j<-1; j++ ) it <<= 1;
        tnm = it;
        del = ( b-a )/tnm;

        // estimate integral
        x = a+0.5*del;
        for( sum=0.0, j=0; j<it; j++, x+=del ) sum+=f( x );
        s = 0.5*( s+( b-a )*sum/tnm );
        return s;
      }
  }

    template <typename T>
    T Quadrature<T>::qTrap( const T eps )
  {
    const int MAXSTEPS=20;
    T olds = 0;

    for( int i=0; i<MAXSTEPS; i++ )
      {
        s = next();
        // always do at least 5 iterations
        if( i>5 )
          // is there convergence?
          if( std::abs( s-olds ) < eps*std::abs( olds )
              or ( s==0 and olds==0 ) )
            return s; // then return the result

        // otherwise store the old result
        olds = s;
        // if we got this far, there is no convergence before the
        // maximum number of steps
      }
    throw( "Maximum number of steps reached before convergence"
           " in the quadrature rule" );
    return s;
  }

    template <typename T>
    T trapVec( array<T> vec, const T h )
  {
    T sum, s;
    int i;
    const int len = vec.shape()[0];

    // assert that the vector is at lest of length 3
    if( len<2 )
      throw std::invalid_argument( "Vector for integration must be of"
                                   " minimum length 3" );
    // add the end points
    s = h/2.0*( vec[0] + vec[len-1] );

    // sum the points
    for( sum=0.0, i=0; i<len-1; i++ )
      sum += h*vec[i+1];

    s += sum;
    return s;
  }

    template <typename T>
    T trapVec( array<T> yvec, array<T> xvec )
  {
    T sum;
    unsigned int i;
    const unsigned int len = yvec.shape()[0];
    if( xvec.shape()[0] != len )
      throw std::invalid_argument( "vectors must be of the same"
				   " length." );

    // assert that the vector is at lest of length 3
    if( len<2 )
      throw std::invalid_argument( "Vector for integration must be of"
                                   " minimum length 3" );

    // sum the points
    for( sum=0.0, i=1; i<len; i++ )
      sum += ( xvec[i]-xvec[i-1] )*( yvec[i]+yvec[i-1] );

    sum/=2.0;
    return sum;
  }

    template class Quadrature<float>;
    template class Quadrature<double>;
}
