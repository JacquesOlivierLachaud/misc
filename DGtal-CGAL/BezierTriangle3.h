/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file BezierTriangle3.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module BezierTriangle3.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(BezierTriangle3_RECURSES)
#error Recursive header files inclusion detected in BezierTriangle3.h
#else // defined(BezierTriangle3_RECURSES)
/** Prevents recursive inclusion of headers. */
#define BezierTriangle3_RECURSES

#if !defined BezierTriangle3_h
/** Prevents repeated inclusion of headers. */
#define BezierTriangle3_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/math/MPolynomial.h>
#include <DGtal/math/linalg/SimpleMatrix.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class BezierTriangle3
  /**
     Description of class 'BezierTriangle3' <p> \brief Aim: This class
     represents a Bezier triangle of degree 3, as a
     polynomial R^2 -> R^M. 

     @tparam TSpace the digital space for images (choose Z2i::Space).

     @tparam M the number of scalar per value, generally 1 for
     gray-level images and 3 for color images.
  */
  template <typename TSpace, int M>
  class BezierTriangle3
  {
  public:
    typedef TSpace                            Space;
    typedef BezierTriangle3< TSpace, M> Self;
    BOOST_CONCEPT_ASSERT(( concepts::CSpace< Space > ));
    BOOST_STATIC_ASSERT (( Space::dimension == 2 ));
    BOOST_STATIC_ASSERT (( M >= 1 ));

    typedef typename Space::Point                 Point;
    typedef typename Space::Vector                Vector;
    typedef typename Space::RealPoint             RealPoint;
    typedef typename Space::RealVector            RealVector;
    typedef typename Space::Size                  Size;
    typedef typename RealVector::Component        Scalar;
    typedef PointVector< M, Scalar >              Value;
    typedef PointVector< 3, Scalar >              RealPoint3;
    struct VectorValue { Value x; Value y; };

    /// A Bezier triangle is a polynomial in three variables.
    typedef MPolynomial<3, Scalar >               BPolynomial;

    /// The Bezier polynomials (one for each dimension of the values).
    BPolynomial _P[ M ];
    
    /// There are 10 bezier points in a Bezier triangle of degree 3.
    /// Their order is inverse lexicographic according to variables
    /// ijk, i.e. b300, b210, b201, b120, b111, b102, b030, b021,
    /// b012, b003
    std::array<RealPoint,10>   _bpoints;

    /// Useful table to get the linear index from the three indices of
    /// a Bezier point.
    inline
    static constexpr int _index[ 4 ][ 4 ][ 4 ]
    = { { { -1, -1, -1, 9 /* b003 */ }, { -1, -1, 8 /* b012 */, -1 },
	  { -1, 7 /* b021 */, -1, -1 }, { 6 /* b030 */, -1, -1, -1 } },
	{ { -1, -1, 5 /* b102 */, -1 }, { -1, 4 /* b111 */, -1, -1 },
	  { 3 /* b120 */, -1, -1, -1 }, { -1, -1, -1, -1 } },
	{ { -1, 2 /* b201 */, -1, -1 }, { 1 /* b210 */, -1, -1, -1 },
	  { -1, -1, -1, -1 }, { -1, -1, -1, -1 } },
	{ { 0 /* b300 */, -1, -1, -1 }, { -1, -1, -1, -1 },
	  { -1, -1, -1, -1 }, { -1, -1, -1, -1 } } };

    /// The indices of the vertices of the Bezier triangle.
    inline
    static constexpr int _vtx[ 3 ] = { 0 /* b300 */, 6 /* b030 */, 9 /* b003 */ };
    /// The values at the ten bezier points, in order b300, b210,
    /// b201, b120, b111, b102, b030, b021, b012, b003
    std::array<Value,10>       _bvalues;

    /// Given i,j,k with i+j+k = 3, return the corresponding index.
    static inline
    int index( int i, int j, int k )
    {
      ASSERT( ( i+j+k == 3 ) && ( i >= 0 ) && ( j >= 0 ) && ( k >= 0 ) );
      return _index[ i ][ j ][ k ];
    }

    static inline
    Scalar det( const RealVector& v, const RealVector& w )
    {
      return v[ 0 ] * w[ 1 ] - v[ 1 ] * w[ 0 ];
    }

    static inline
    Scalar bound( Scalar s, Scalar lo = 0.0, Scalar up = 255.0 )
    { return std::min( up, std::max( lo, s ) ); }
    
    static inline
    Value bound( const Value& s,
		 const Value& lo = Value { 0.0, 0.0, 0.0 },
		 const Value& up = Value { 255.0, 255.0, 255.0 } )
    {
      Value r;
      for ( int m = 0; m < M; ++m )
	r[ m ] = bound( s[ m ], lo[ m ], up[ m ] );
      return r;
    }
    
    // ----------------------- Standard services ------------------------------
  public:
  
    /// Destructor.
    ~BezierTriangle3() = default;
    /// Default constructor. The object is invalid.
    BezierTriangle3() = default;
    /// Default copy constructor. 
    BezierTriangle3( const BezierTriangle3& ) = default;
    /// Default move constructor. 
    BezierTriangle3( BezierTriangle3&& ) = default;
    /// Default assignment.
    BezierTriangle3& operator=( BezierTriangle3&& ) = default;

    /// Constructeur from vertex points and all values in order b200,
    /// b020, b002, b110, b011, b101
    BezierTriangle3( const std::array<RealPoint,3>& vtcs,
		     const std::array<Value,10>&    bvalues )
      : _bvalues( bvalues )
    {
      getBezierPoints( _bpoints, vtcs );
      makePolynomials();
    }

    /// Initializes all Bezier points from the three vertex points of
    /// the Bezier triangle.
    static void getBezierPoints
    ( std::array<RealPoint,10>& bp, const std::array<RealPoint,3>& vtcs )
    {
      for ( int i = 0; i < 4; ++i ) {
	for ( int j = 0; j < 4; ++j ) {
	  if ( ( i + j ) > 3 ) continue;
	  int k = 3 - ( i + j ); // i+j+k == 3
	  int n = index( i, j, k );
	  bp[ n ] = ( (Scalar) i * vtcs[ 0 ]
		      + (Scalar) j * vtcs[ 1 ]
		      + (Scalar) k * vtcs[ 2 ] )
	    / 3.0 ;
	  std::cout << "bp[" << n << "]=" << bp[ n ] << std::endl;
	}
      }
    }

    // Definition of a (local) gradient operator that assigns vectors to triangles
    static
    VectorValue grad( const std::array<RealPoint,3>& pos,
		      const std::array<Value,3>&     val )
    {
      // [ yj-yk yk-yi yi-yk ] * [ ui ]
      // [ xk-xj xi-xk xj-xi ]   [ uj ]
      //                         [ uk ]
      const RealPoint& pi = pos[ 0 ];
      const RealPoint& pj = pos[ 1 ];
      const RealPoint& pk = pos[ 2 ];
      const Value&     ui = val[ 0 ];
      const Value&     uj = val[ 1 ];
      const Value&     uk = val[ 2 ];
      VectorValue G;
      for ( int m = 0; m < M; ++m ) {
	G.x[ m ] = ui[ m ] * ( pj[ 1 ] - pk[ 1 ] )
	  +        uj[ m ] * ( pk[ 1 ] - pi[ 1 ] )
	  +        uk[ m ] * ( pi[ 1 ] - pj[ 1 ] );
	G.y[ m ] = ui[ m ] * ( pk[ 0 ] - pj[ 0 ] )
	  +        uj[ m ] * ( pi[ 0 ] - pk[ 0 ] )
	  +        uk[ m ] * ( pj[ 0 ] - pi[ 0 ] );
      }
      Scalar area = ( pj - pi ).crossProduct( pk - pi ).norm();
      G.x /= area;
      G.y /= area;
      return G;
    }

    // gradient operator that assigns vectors to triangles
    static
    RealVector grad( const std::array<RealPoint,3>& pos,
		     const std::array<Scalar,3>&    val )
    {
      // [ yj-yk yk-yi yi-yk ] * [ ui ]
      // [ xk-xj xi-xk xj-xi ]   [ uj ]
      //                         [ uk ]
      const RealPoint& pi = pos[ 0 ];
      const RealPoint& pj = pos[ 1 ];
      const RealPoint& pk = pos[ 2 ];
      const Scalar&    ui = val[ 0 ];
      const Scalar&    uj = val[ 1 ];
      const Scalar&    uk = val[ 2 ];
      RealVector G = { ui * ( pj[ 1 ] - pk[ 1 ] )
		       + uj * ( pk[ 1 ] - pi[ 1 ] )
		       + uk * ( pi[ 1 ] - pj[ 1 ] ),
		       ui * ( pk[ 0 ] - pj[ 0 ] )
		       + uj * ( pi[ 0 ] - pk[ 0 ] )
		       + uk * ( pj[ 0 ] - pi[ 0 ] ) };
      Scalar area = ( pj - pi ).crossProduct( pk - pi ).norm();
      return G / area;
    }
    /// From a gradient, positions and one value at pos[ 2 ], find the
    /// values at pos[ 0 ] and pos[ 1 ].
    static
    void getVab( Scalar& v0, Scalar& v1,
		 const Scalar                    v2,
		 const std::array<RealPoint,3>& pos,
		 const RealVector&                G )
    {
      v0 = v2 + G[ 0 ] * ( pos[ 0 ][ 0 ] - pos[ 2 ][ 0 ] )
	+ G[ 1 ] * ( pos[ 0 ][ 1 ] - pos[ 2 ][ 1 ] );
      v1 = v2 + G[ 0 ] * ( pos[ 1 ][ 0 ] - pos[ 2 ][ 0 ] )
	+ G[ 1 ] * ( pos[ 1 ][ 1 ] - pos[ 2 ][ 1 ] );
    }
    /// From a gradient, positions and one value at pos[ 2 ], find the
    /// values at pos[ 0 ] and pos[ 1 ].
    static
    void getVab( Value& v0, Value& v1,
		 const Value                     v2,
		 const std::array<RealPoint,3>& pos,
		 const RealVector&                G )
    {
      for ( int m = 0; m < M; ++m ) {
	getVab( v0[ m ], v1[ m ], v2[ m ], pos, G );
      }
    }

    static
    Value midValue( const std::array<Value,3>      & tv )
    {
      Value v;
      for ( int m = 0; m < M; ++m )
	v[ m ] = 0.5 *
	  ( std::max( tv[ 0 ][ m ], std::max( tv[ 1 ][ m ], tv[ 2 ][ m ] ) )
	    + std::min( tv[ 0 ][ m ], std::min( tv[ 1 ][ m ], tv[ 2 ][ m ] ) ) );
      return v;
    }
    /// Constructeur from Bezier points b200, b020, b002, and
    /// orthogonality conditions at b110, b011, b101 and values at
    /// b200, b020, b002
    ///
    /// @note The orthogonality conditions is that at b110 (resp. b011
    /// and b101), the gradient of the Bezier polynomial should be
    /// orthogonal to orthvecs[ 0 ] (resp. orthvecs[ 1 ] and orthvecs[
    /// 2 ]).
    BezierTriangle3( const std::array<RealPoint,3>  & vtcs,
		     const std::array<Value,3>      & vtcs_values,
		     const std::array<RealVector,3> & orthvecs,
		     const Value                    & c,
		     const std::array<RealPoint,3>  & dpoints,
		     const std::array<Value,3>      & dvalues )
    {
      getBezierPoints( _bpoints, vtcs );
      _bvalues[ index( 3, 0, 0 ) ] = vtcs_values[ 0 ];
      _bvalues[ index( 0, 3, 0 ) ] = vtcs_values[ 1 ];
      _bvalues[ index( 0, 0, 3 ) ] = vtcs_values[ 2 ];
      _bvalues[ index( 1, 1, 1 ) ] = c;
      // We must determine the six remaining values, that lie along
      // edges b300 -- b030 -- b003 -- b300.
      typedef SimpleMatrix<Scalar,2,2>      Matrix;
      typedef typename Matrix::ColumnVector ColumnVector;
      const int        ie[ 3 ] = { index(3,0,0), index(0,3,0), index(0,0,3) };
      const int        ia[ 3 ] = { index(2,1,0), index(0,2,1), index(1,0,2) };
      const int        ib[ 3 ] = { index(1,2,0), index(0,1,2), index(2,0,1) };
      const int        ic = index(1,1,1);
      const RealPoint& Xc = _bpoints[ ic ];
      for ( int n = 0; n < 3; ++n ) {
	const RealPoint& Xa = _bpoints[ ia[ n ] ];
	const RealPoint& Xb = _bpoints[ ib[ n ] ];
	const RealPoint& Xd = dpoints[ n ];
	const RealVector& G = orthvecs[ n ];
	const Scalar  alpha = G[ 0 ];
	const Scalar   beta = G[ 1 ];
	std::cout << "Xa=" << Xa << " Xb=" << Xb << " Xc=" << Xc
		  << " Xd=" << Xd << " alpha=" << alpha << " beta=" << beta
		  << std::endl;
	Value va_c, vb_c, va_d, vb_d;
	getVab( va_c, vb_c, c, { Xa, Xb, Xc }, G ); 
	getVab( vb_d, va_d, dvalues[ n ], { Xb, Xa, Xd }, G ); 
	std::cout << "va_c=" << va_c << " va_d=" << va_d << std::endl;
	std::cout << "vb_c=" << vb_c << " vb_d=" << vb_d << std::endl;
	_bvalues[ ia[ n ] ] = bound( 0.5 * ( va_c + va_d ) );
	_bvalues[ ib[ n ] ] = bound( 0.5 * ( vb_c + vb_d ) );
	// Scalar Aabc = (Xb-Xa).crossProduct(Xc-Xa).norm();
	// Scalar Abad = (Xa-Xb).crossProduct(Xd-Xb).norm();
	// Matrix Tabc =
	//   { ( Xb[ 1 ] - Xc[ 1 ] ), ( Xc[ 0 ] - Xb[ 0 ] ),
	//     ( Xc[ 1 ] - Xa[ 1 ] ), ( Xa[ 0 ] - Xc[ 0 ] ) };
	// Matrix Tbad =
	//   { ( Xa[ 1 ] - Xd[ 1 ] ), ( Xd[ 0 ] - Xa[ 0 ] ),
	//     ( Xd[ 1 ] - Xb[ 1 ] ), ( Xb[ 0 ] - Xd[ 0 ] ) };
	// Scalar dTabc = Tabc.determinant();
	// std::cout << "Tabc=" << Tabc << std::endl;
	// std::cout << "dTabc=" << dTabc << std::endl;
	// Scalar dTbad = Tbad.determinant();
	// std::cout << "Tbad=" << Tbad << std::endl;
	// std::cout << "dTbad=" << dTbad << std::endl;
	// Matrix Tabc_inv = Tabc.inverse();
	// Matrix Tbad_inv = Tbad.inverse();
	// for ( int m = 0; m < M; ++m ) {
	//   ColumnVector Bc = { Aabc * alpha - ( Xa[ 1 ] - Xb[ 1 ] ) * c[ m ],
	// 		      Aabc * beta  - ( Xb[ 0 ] - Xa[ 0 ] ) * c[ m ] }; 
	//   ColumnVector Bd = { Abad * alpha - ( Xb[ 1 ] - Xa[ 1 ] ) * dvalues[ n ][ m ],
	// 		      Abad * beta  - ( Xa[ 0 ] - Xb[ 0 ] ) * dvalues[ n ][ m ] };
	//   ColumnVector Rc = Tabc_inv * Bc;
	//   ColumnVector Rd = Tbad_inv * Bd;
	//   // NB: inversion a and b between abc and bad
	//   std::cout << "Rc=" << Rc << " Rd=" << Rd << std::endl;
	//   _bvalues[ ia[ n ] ][ m ] = bound( 0.5 * ( Rc[ 0 ] + Rd[ 1 ] ) );
	//   _bvalues[ ib[ n ] ][ m ] = bound( 0.5 * ( Rc[ 1 ] + Rd[ 0 ] ) );
	// }
	// if ( dT == 0 ) {
	//   // We have to find some optimal value according to least-square sense
	//   Scalar delta = ( T( 1, 0 ) != 0.0 )
	//     ? T( 1, 0 ) / T( 0, 0 ) : T( 1, 1 ) / T( 0, 1 );
	//   _bvalues[ ia[ n ] ] = _bvalues[ ib[ n ] ] =
	//     bound( ( c + dvalues[ n ] ) / ( 1.0 + delta ) );
	// } else {
	//   Matrix    Tinv = T.inverse();
	//   Scalar   delta = - ( alpha * ( Xa[ 1 ] - Xb[ 1 ] )
	// 		       + beta * ( Xb[ 0 ] - Xa[ 0 ] ) );
	//   for ( int m = 0; m < M; ++m ) {
	//     ColumnVector B = { delta * c[ m ], -delta * dvalues[ n ][ m ] };
	//     ColumnVector R = Tinv * B;
	//     _bvalues[ ia[ n ] ][ m ] = bound( R[ 0 ] );
	//     _bvalues[ ib[ n ] ][ m ] = bound( R[ 1 ] );
	//   }
	// }
      } // for ( int n = 0; n <3; ++n ) {
      for ( int n = 0; n < 10; ++n ) {
	std::cout << "bv[" << n << "]=" << _bvalues[ n ] << std::endl;
      }
      makePolynomials();
    }

    /// @return the i-th Bezier control point
    const RealPoint& b( int i ) const { return _bpoints[ i ]; }

    /// @return the i-th Bezier value point
    const Value& v( int i ) const { return _bvalues[ i ]; }

    /// @return the i-th triangle vertex for i in 0,1,2.
    const RealPoint& vertex( int i ) const { return _bpoints[ _vtx[ i ] ]; }

    /// @param[in] p the point where you want to evaluate the
    /// BezierTriangle.
    ///
    /// @return the value at this point
    Value operator()( const RealPoint& p ) const
    {
      return operator()( barycentric( p ) );
    }
    
    /// @param[in] b the barycentric coordinates where you want to
    /// evaluate the BezierTriangle: (1,0,0) is b200, (0,1,0) is b020,
    /// (0,0,1) is b002, (a,b,c) with a+b+c=1, 0<=a<=1, 0<=b<=1,
    /// 0<=c<=1 is anywhere between these points.
    ///
    /// @return the value at this point
    Value operator()( const RealPoint3& bc ) const
    {
      return operator()( bc[ 0 ], bc[ 1 ], bc[ 2 ] );
    }

    /// @param[in] r,s,t the barycentric coordinates where you want to
    /// evaluate the BezierTriangle: (1,0,0) is b300, (0,1,0) is b030,
    /// (0,0,1) is b003, (a,b,c) with a+b+c=1, 0<=a<=1, 0<=b<=1,
    /// 0<=c<=1 is anywhere between these points.
    ///
    /// @return the value at this point
    Value operator()( Scalar r, Scalar s, Scalar t ) const
    {
      return Value( _P[ 0 ]( r )( s )( t ),
		    _P[ 1 ]( r )( s )( t ),
		    _P[ 2 ]( r )( s )( t ) );
    }

    /// @return the barycentric coordinates of point \a p in this triangle.
    RealPoint3 barycentric( const RealPoint& p ) const
    {
      RealPoint3 bc = { det( b( _vtx[ 1 ] ) - p, b( _vtx[ 2 ] ) - p ),
			det( b( _vtx[ 2 ] ) - p, b( _vtx[ 0 ] ) - p ),
			det( b( _vtx[ 0 ] ) - p, b( _vtx[ 1 ] ) - p ) };
      return bc / ( bc[ 0 ] + bc[ 1 ] + bc[ 2 ] );
    }

    /// @param bc any barycentric coordinates (sum is 1).
    ///
    /// @return 'true' iff the given barycentric coordinates \a bc
    /// indicates a point inside or on the boundary of the triangle.
    bool isInTriangle( const RealPoint3& bc ) const
    {
      return ( 0 <= bc[ 0 ] ) && ( bc[ 0 ] <= 1 )
	&&   ( 0 <= bc[ 1 ] ) && ( bc[ 1 ] <= 1 )
	&&   ( 0 <= bc[ 2 ] ) && ( bc[ 2 ] <= 1 );
    }

      
    
    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:

    /// Builds the Bezier polynomials from the current values
    void makePolynomials()
    {
      for ( int m = 0; m < M; ++m ) {
	for ( int i = 0; i < 4; ++i ) {
	  for ( int j = 0; j < 4; ++j ) {
	    if ( ( i + j ) > 3 ) continue;
	    int k = 3 - ( i + j ); // i+j+k == 3
	    int n = index( i, j, k );
	    _P[ m ] += ( ( n == 4 ? 3.0 : 1.0 ) * _bvalues[ n ][ m ] )
	      * mmonomial<Scalar>( i, j, k );
	  }
	}
	std::cout << "_P[ " << m << " ] = " << _P[ m ] << std::endl;
      }
    }
    
    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class BezierTriangle3

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined BezierTriangle3_h

#undef BezierTriangle3_RECURSES
#endif // else defined(BezierTriangle3_RECURSES)
