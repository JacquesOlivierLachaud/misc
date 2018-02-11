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
 * @file ImageConnecter.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module ImageConnecter.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ImageConnecter_RECURSES)
#error Recursive header files inclusion detected in ImageConnecter.h
#else // defined(ImageConnecter_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ImageConnecter_RECURSES

#if !defined ImageConnecter_h
/** Prevents repeated inclusion of headers. */
#define ImageConnecter_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/CImage.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class ImageConnecter
  /**
     Description of class 'ImageConnecter' <p> \brief Aim: This class 
     optimizes the initial connectedness (4 or 8) of a given image.

     In this class, colors are represented as a 3-vector with
     components in 0..255.

     @tparam TImage the type of 2D image.
  */
  template <typename TImage>
  class ImageConnecter {
  public:
    typedef TImage                            Image;
    typedef ImageConnecter< TImage >          Self;
    BOOST_CONCEPT_ASSERT(( concepts::CImage< Image > ));

    typedef typename Image::Domain                Domain;
    typedef typename Image::Value                 Value;
    typedef typename Domain::Space                Space;
    typedef typename Space::Point                 Point;
    typedef typename Space::Vector                Vector;
    typedef typename Space::RealVector            RealVector;
    typedef typename Space::Size                  Size;
    typedef typename Vector::Component            Integer;
    typedef typename RealVector::Component        Scalar;
    typedef std::function<Scalar( Value, Value )> Comparator; // 0 is equal

    BOOST_STATIC_ASSERT (( Space::dimension == 2 ));
    
    struct Quad {
      Point base;
      Size  min;
      Size  max;
      bool operator<( const Quad& other ) const
      {
	return ( ( max - min ) > ( other.max - other.min ) )
	    || ( ( ( max - min ) == ( other.max - other.min ) )
		 && ( min < other.min ) );
      }
      // { return ( min < other.min )
      // 	  || ( ( min == other.min )
      // 	       && ( max > other.max ) ); }
    };
    
    // Union-find data structure.
    struct Element {
      Point    point;
      Size     rank;
      Element* father;
      Size     nb;
      Element( Element* self = nullptr, Point p = Point( 0, 0 ) )
	: point( p ), rank( 0 ), father( self ), nb ( 1 ) {}
    };

    // return the root of the tree containing e
    static Element* find( Element* e )
    {
      if ( e != e->father )
	e->father = find( e->father );
      return e->father;
    }
    // link the two trees of roots x and y.
    static void link( Element* x, Element* y )
    {
      if ( x->rank > y->rank ) {
	y->father = x;
	x->nb    += y->nb; 
      }
      else {
	x->father = y;
	y->nb    += x->nb; 
	if ( x->rank == y->rank ) y->rank += 1;
      }
    }
    // Given two elements x and y in different trees, makes the union
    // of the two trees.
    static void merge( Element* x, Element* y )
    {
      link( find( x ), find( y ) );
    }
    enum Configuration { Default = 0, Diagonal00_11, Diagonal10_01 };

    std::vector< Configuration > connections;
    std::vector< Element >       labels;
    Integer                      width;
  public:
    /**
       Constructor. 
    */
    ImageConnecter()
    {}
    
    Configuration howConnected( Point p ) const
    {
      return connections[ p[ 1 ] * width + p[ 0 ] ];
    }
    void init( const Image& I, Comparator comp, Scalar same )
    {
      // Creates 4-connected components
      const Domain& domain = I.domain();
      const Domain rdomain ( domain.lowerBound(),
			     domain.upperBound() - Point::diagonal( 1 ) );
      const Vector& extent = I.extent();
      width  = extent[ 0 ];
      labels.resize( domain.size() );
      connections.resize( domain.size() );
      Size i = 0;
      for ( auto p : domain ) {
	Element& self = labels[ i++ ];
	self = Element( &self, p );
      }
      // Scan and merge
      Point up = domain.upperBound();
      for ( auto p : domain ) {
	Value v = I( p );
	if ( p[ 0 ] < up[ 0 ] ) {
	  Point q = p + Vector( 1, 0 );
	  if ( comp( v, I( q ) ) <= same ) {
	    Element* e1 = & labels[ p[ 1 ] * width + p[ 0 ] ];
	    Element* e2 = & labels[ q[ 1 ] * width + q[ 0 ] ];
	    if ( find( e1 ) != find( e2 ) ) {
	      merge( e1, e2 );
	      // std::cout << "merge " << e1->nb << " " << e2->nb << std::endl;
	    }
	  }
	}
	if ( p[ 1 ] < up[ 1 ] ) {
	  Point q = p + Vector( 0, 1 );
	  if ( comp( v, I( q ) ) <= same ) {
	    Element* e1 = & labels[ p[ 1 ] * width + p[ 0 ] ];
	    Element* e2 = & labels[ q[ 1 ] * width + q[ 0 ] ];
	    if ( find( e1 ) != find( e2 ) ) {
	      merge( e1, e2 );
	      // std::cout << "merge " << e1->nb << " " << e2->nb << std::endl;
	      //std::cout << "merge " << p << " " << q << std::endl;
	    }
	  }
	}
      } // end scan for 4-connectedness
      // Connect diagonals around big regions.
      std::vector<Quad> quads( rdomain.size() );
      for ( auto p : rdomain ) {
	const Point p00 = p;
	const Point p10 = p + Vector( 1, 0 );
	const Point p01 = p + Vector( 0, 1 );
	const Point p11 = p + Vector( 1, 1 );
	Element*    e00 = find( & labels[ p00[ 1 ] * width + p00[ 0 ] ] );
	Element*    e10 = find( & labels[ p10[ 1 ] * width + p10[ 0 ] ] );
	Element*    e01 = find( & labels[ p01[ 1 ] * width + p01[ 0 ] ] );
	Element*    e11 = find( & labels[ p11[ 1 ] * width + p11[ 0 ] ] );
	Quad          q = { p00,
			    std::min( std::min( e00->nb, e10->nb ),
				      std::min( e01->nb, e11->nb ) ),
			    std::max( std::max( e00->nb, e10->nb ),
				      std::max( e01->nb, e11->nb ) ) };
	quads.push_back( q );
      }
      std::sort( quads.begin(), quads.end() );
      
      // Connect diagonals around big regions.
      Size nb00_11 = 0;
      Size nb10_01 = 0;
      for ( auto q : quads ) {
	//	std::cout << q.base << " " << q.min << " " << q.max << std::endl;
	const Point p00 = q.base;
	const Point p10 = p00 + Vector( 1, 0 );
	const Point p01 = p00 + Vector( 0, 1 );
	const Point p11 = p00 + Vector( 1, 1 );
	Element*    e00 = find( & labels[ p00[ 1 ] * width + p00[ 0 ] ] );
	Element*    e10 = find( & labels[ p10[ 1 ] * width + p10[ 0 ] ] );
	Element*    e01 = find( & labels[ p01[ 1 ] * width + p01[ 0 ] ] );
	Element*    e11 = find( & labels[ p11[ 1 ] * width + p11[ 0 ] ] );
	const Value v00 = I( p00 );
	const Value v10 = I( p10 );
	const Value v01 = I( p01 );
	const Value v11 = I( p11 );
	Scalar   s00_11 = comp( v00, v11 );
	Scalar   s10_01 = comp( v10, v01 );
	Configuration c = Default;
	if ( ( s00_11 <= same ) && ( same < s10_01 ) )
	  c = Diagonal00_11;
	else if ( ( s10_01 <= same ) && ( same < s00_11 ) )
	  c = Diagonal10_01;
	else if ( ( std::min( e00->nb, e11->nb ) == q.min ) 
		  && ( s00_11 <= same ) )
	  c = Diagonal00_11;
	else if ( ( std::min( e10->nb, e01->nb ) == q.min ) 
		  && ( s10_01 <= same ) )
	  c = Diagonal10_01;

	if ( c == Diagonal00_11 ) nb00_11++;
	if ( c == Diagonal10_01 ) nb10_01++;
	connections[ p00[ 1 ] * width + p00[ 0 ] ] = c;
      }
      trace.info() << "nb00_11=" << nb00_11
		   << " nb10_01=" << nb10_01;
    }
    
    /// Destructor.
    ~ImageConnecter()
    {
    }
    
    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ImageConnecter

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ImageConnecter_h

#undef ImageConnecter_RECURSES
#endif // else defined(ImageConnecter_RECURSES)
