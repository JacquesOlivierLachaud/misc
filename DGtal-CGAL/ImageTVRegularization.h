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
 * @file ImageTVRegularization.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module ImageTVRegularization.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ImageTVRegularization_RECURSES)
#error Recursive header files inclusion detected in ImageTVRegularization.h
#else // defined(ImageTVRegularization_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ImageTVRegularization_RECURSES

#if !defined ImageTVRegularization_h
/** Prevents repeated inclusion of headers. */
#define ImageTVRegularization_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/kernel/domains/Linearizer.h>
#include <DGtal/helpers/StdDefs.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class ImageTVRegularization
  /**
     Description of class 'ImageTVRegularization' <p> \brief Aim:
     This class regularizes an input image with respect to its total variation.
     @tparam TSurface any type of digital surface.
  */
  template <typename TSpace, int M>
  class ImageTVRegularization
  {
  public:
    typedef TSpace                            Space;
    typedef ImageTVRegularization< TSpace, M> Self;
    // BOOST_CONCEPT_ASSERT(( concepts::CSpace< Space > ));
    BOOST_STATIC_ASSERT (( Space::dimension <= 3 ));
    BOOST_STATIC_ASSERT (( M >= 1 ));

    typedef typename Space::Point                 Point;
    typedef typename Space::Vector                Vector;
    typedef typename Space::RealPoint             RealPoint;
    typedef typename Space::RealVector            RealVector;
    typedef typename Space::Size                  Size;
    typedef HyperRectDomain<Space>                Domain;
    typedef typename RealVector::Component        Scalar;
    typedef PointVector< M, Scalar >              Value;
    typedef std::array< Value, Space::dimension > VectorValue;
    ///static constants to store the dimension.
    static const Dimension N = Space::dimension;
    typedef std::vector<Scalar>                   ScalarForm;
    typedef std::vector<Value>                    ValueForm;
    typedef std::vector<VectorValue>              VectorValueForm;

    /// The image values at each vertex
    ValueForm            _I;
    /// The domain of the image and of the computations.
    Domain               _domain;
    /// The extent of the domain of the image and of the computations.
    Vector               _extent;
    /// The regularized values at each vertex
    ValueForm            _u;
    /// The TV-regularized vectors
    VectorValueForm      _p;
    
    // ----------------------- Standard services ------------------------------
  public:
  
    /**
     * Destructor.
     */
    ~ImageTVRegularization() {}

    /// Default constructor. The object is invalid.
    ImageTVRegularization() : _domain( Point(), Point() ) {}

    /// Functor used to feed the TV with a color image (M should be 3).
    struct Color2ValueFunctor {
      Scalar operator()( unsigned int color, unsigned int m ) const {
	switch( m ) {
	case 0: return (color >> 16) & 0xff;
	case 1: return (color >> 8) & 0xff;
	default: return color & 0xff;
	}
      }
    };
    /// Functor used to feed the TV with a gray-level image (M should be 1).
    struct GrayLevel2ValueFunctor {
      Scalar operator()( unsigned int color, unsigned int m ) const {
	return color & 0xff;
      }
    };
    
    /// @tparam Functor the type of function: Image Value x int --> Scalar,
    /// where int is the dimension in the Value.
    /// @see Color2ValueFunctor
    /// @see GrayLevel2ValueFunctor
    template <typename Image, typename Functor>
    void init( const Image& I, Functor f )
    {
      _I.reserve( I.size() );
      for ( auto val_I : I ) {
	Value v;
	for ( unsigned int m = 0; m < M; ++m )
	  v[ m ] = f( val_I, m );
	_I.push_back( v );
      }
      _u = _I;                  // u = image at initialization
      _p.resize( _I.size() );   // p = 0     at initialization
      _domain = I.domain();
      _extent = I.extent();
    }

    /// Bounds the scalar value in [0,255] and rounds it to the nearest integer.
    static unsigned int dig( Scalar v )
    {
      return (unsigned int ) round( std::min( 255.0, std::max( 0.0, v ) ) );
    }

    /// Functor for outputing the result into a color image.
    struct Value2ColorFunctor {
      unsigned int operator()( const Value& v ) const
      {
	if ( M <= 2 ) 
	  return ( dig( v[ 0 ] ) << 16 )
	    +    ( dig( v[ 0 ] ) << 8 )
	    +    ( dig( v[ 0 ] ) );
	else 
	  return ( dig( v[ 0 ] ) << 16 )
	    +    ( dig( v[ 1 ] ) << 8 )
	    +    ( dig( v[ 2 ] ) );
      }
    };
    
    /// Functor for outputing the result into a gray-level image.
    struct Value2GrayLevelFunctor {
      unsigned int operator()( const Value& v ) const
      {
	unsigned int gl = 0;
	for ( unsigned int m = 0; m < M; ++m )
	  gl += dig( v[ m ] );
	return gl / M;
      }
    };
    
    /// @see Value2ColorFunctor
    /// @see Value2GrayLevelFunctor
    template <typename Image, typename Functor>
    bool outputU( Image& J, Functor f ) const
    {
      Size i = 0;
      for ( auto & val_J : J ) {
	val_J = f( _u[ i ] );
	i++;
      }
      return i == _u.size();
    }

    
    /// Linearize the point to an index.
    Size index( const Point &aPoint ) const
    {
      return Linearizer<Domain, ColMajorStorage>::getIndex( aPoint, _extent );
    }
    /// Delinearize the index to a point.
    Point point( const Size index ) const
    {
      return Linearizer<Domain, ColMajorStorage>::getPoint( index, _extent );
    }

    // ---------------- Calculus norms and operators ------------------
    
    /// Computes the square of x.
    static Scalar square( Scalar x ) { return x*x; }

    /// The norm for Value (Euclidean)
    Scalar normX( const Value& v ) const
    {
      Scalar x = square( v[ 0 ] );
      for ( unsigned int m = 1; m < M; ++m )
	x += square( v[ m ] );
      return sqrt( x );
    }

    /// The norm for VectorValue (TV and VTV)
    Scalar normY( const VectorValue& vv ) const
    {
      Scalar xx = 0.0;
      for ( unsigned int n = 0; n < N; ++n )
	for ( unsigned int m = 0; m < M; ++m )
	  xx += square( vv[ n ][ m ] );
      return sqrt( xx );
    }

    /// @return the vector of the norms of each vector in p (one per vertex). 
    ScalarForm norm( const VectorValueForm& p ) const
    {
      ScalarForm S( p.size() );
      for ( Size i = 0; i < p.size(); i++ )
	S[ i ] = normY( p[ i ] );
      return S;
    }

    // Definition of a global gradient operator that assigns vectors to triangles
    VectorValueForm grad( const ValueForm& u ) const
    { // it suffices to traverse all (valid) triangles.
      VectorValueForm G( u.size() );
      Size i = 0;
      for ( auto p : _domain ) {
	G[ i ] = grad( u, p, i );
	i++;
      }
      return G;
    }

    inline
    VectorValue grad( const ValueForm& u, const Point& p ) {
      return grad( u, p, index( p ) );
    } 

    inline
    VectorValue grad( const ValueForm& u, Size i ) {
      return grad( u, point( i ), i );
    } 

    // Definition of a (local) gradient operator that assigns vectors to triangles
    VectorValue grad( const ValueForm& u, const Point& p,
		      const Size i ) const
    {
      VectorValue G;
      for ( unsigned int n = 0; n < N; ++n )
	if ( p[ n ] < _extent[ n ] - 1 )
	  for ( unsigned int m = 0; m < M; ++m )
	    G[ n ][ m ] = u[ index( p + Point::base( n ) ) ][ m ] - u[ i ][ m ];
      // otherwise zero.
      return G;
    }

    // Definition of a (global) divergence operator that assigns
    // scalars to vertices from a vector field.
    ValueForm div( const VectorValueForm& G ) const
    {
      ValueForm S( G.size() );
      for ( Size i = 0; i < G.size(); ++i )
	S[ i ] = div( G, point( i ), i );
      return S;
    }

    Value div( const VectorValueForm& G, const Point& p ) const {
      return div( G, p, index( p ) );
    }
    Value div( const VectorValueForm& G, Size i ) const {
      return div( G, point( i ), i );
    }

    // Definition of a (local) divergence operator that assigns
    // scalars to vertices from a vector field.
    Value div( const VectorValueForm& G, const Point& p,
	       const Size i ) const
    {
      Value v;
      for ( unsigned int n = 0; n < N; ++n )
	if ( p[ n ] < _extent[ n ] - 1 ) {
	  for ( unsigned int m = 0; m < M; ++m )
	    v[ m ] += G[ i ][ n ][ m ];
	} else if ( p[ n ] > 0 ) {
	  Size i_1 = index( p - Point::base( n ) );
	  for ( unsigned int m = 0; m < M; ++m )
	    v[ m ] -= G[ i_1 ][ n ][ m ];
	}
      return v;
    }

    /// @return the scalar form lambda.u
    ValueForm multiplication( Scalar lambda, const ValueForm& u ) const
    {
      ValueForm S = u;
      for ( Size i = 0; i < S.size(); ++i )
	S[ i ] *= lambda;
      return S;
    }
    
    /// @return the scalar form u - v
    ValueForm subtraction( const ValueForm& u, const ValueForm& v ) const
    {
      ValueForm S = u;
      for ( Size i = 0; i < S.size(); ++i )
	S[ i ] -= v[ i ];
      return S;
    }

    /// u -= v
    void subtract( ValueForm& u, const ValueForm& v ) const
    {
      for ( Size i = 0; i < u.size(); ++i )
	u[ i ] -= v[ i ];
    }

    /// @return the scalar form a.u + b.v
    ValueForm combination( const Scalar a, const ValueForm& u,
			   const Scalar b, const ValueForm& v ) const
    {
      ValueForm S( u.size() );
      for ( Size i = 0; i < S.size(); ++i )
	S[ i ] = a * u[ i ] + b * v[ i ];
      return S;
    }

    /// @return the Total Variation of _u (current approximation of image _I).
    Scalar energyTV() const
    {
      ScalarForm energies = norm( grad( _u ) );
      Scalar            E = 0.0;
      for ( Scalar e_vtx : energies ) E += e_vtx;
      return E;
    }
    
    /// Does one pass of TV regularization (u, p and I must have the
    /// meaning of the previous iteration).
    /// @note Chambolle, Pock primal-dual algorithm 1
    Scalar optimize( Scalar lambda,
		     Scalar dt = 0.248, Scalar tol = 0.01, int max_iter = 15 )
    {
      trace.info() << "TV( u ) = " << energyTV() << std::endl;
      //trace.info() << "lambda.f" << std::endl;
      VectorValueForm p( _p.size() ); // this is p^{n+1}
      ValueForm      lf = multiplication( lambda, _I );
      Scalar     diff_p = 0.0;
      int          iter = 0; // iteration number
      do {
	// trace.info() << "div( p ) - lambda.f" << std::endl;
	ValueForm        dp_lf = subtraction( div( _p ), lf );
	// trace.info() << "G := grad( div( p ) - lambda.f)" << std::endl;
	VectorValueForm gdp_lf = grad( dp_lf );
	// trace.info() << "N := | G |" << std::endl;
	ScalarForm     ngdp_lf = norm( gdp_lf );
	// trace.info() << "p^n+1 := ( p + dt * G ) / ( 1 + dt | G | )" << std::endl;
	diff_p = 0.0;
	for ( Size i = 0; i < p.size(); i++ ) {
	  Scalar alpha = 1.0 / ( 1.0 + dt * ngdp_lf[ i ] );
	  if ( alpha <= 0.0 )
	    trace.warning() << "Face " << i << " alpha=" << alpha << std::endl;
	  VectorValue delta;
	  for ( unsigned int n = 0; n < N; ++n ) {
	    p[ i ][ n ] = alpha * ( _p[ i ][ n ] + dt * gdp_lf[ i ][ n ] );
	    delta[ n ]  = p[ i ][ n ] - _p[ i ][ n ];
	  }
	  diff_p = std::max( diff_p, normY( delta ) );
	}
	trace.info() << "Iter n=" << (iter++) << " diff_p=" << diff_p
		     << " tol=" << tol << std::endl;
	std::swap( p, _p );
      } while ( ( diff_p > tol ) && ( iter < max_iter ) );
      _u = combination( 1.0, _I, -1.0/lambda, div( _p ) );
      trace.info() << "TV( u ) = " << energyTV() << std::endl;
      return diff_p;
    }

    
    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ImageTVRegularization

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ImageTVRegularization_h

#undef ImageTVRegularization_RECURSES
#endif // else defined(ImageTVRegularization_RECURSES)
