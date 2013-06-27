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
 * @file Triangulation3DHelper.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2013/06/27
 *
 * Header file for module Triangulation3DHelper.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(Triangulation3DHelper_RECURSES)
#error Recursive header files inclusion detected in Triangulation3DHelper.h
#else // defined(Triangulation3DHelper_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Triangulation3DHelper_RECURSES

#if !defined Triangulation3DHelper_h
/** Prevents repeated inclusion of headers. */
#define Triangulation3DHelper_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/base/ConstAlias.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class Triangulation3DHelper
  /**
     Description of template class 'Triangulation3DHelper' <p>
     \brief Aim: A class that gives a few methods for moving in a Triangulation3.
     
     @tparam TTriangulation3 any kind of CGAL 3D triangulation.
     @tparam TKernel3 any kind of CGAL 3D Kernel.
  */
  template <typename TTriangulation3, typename TKernel3>
  class Triangulation3DHelper
  {
    // ----------------------- public types ------------------------------
  public:
    typedef TTriangulation3                                  Triangulation3;
    typedef TKernel3                                          Kernel3;
    typedef Triangulation3DHelper<Triangulation3, Kernel3>       Self;
    typedef typename Triangulation3::Vertex_handle           VertexHandle;
    typedef typename Triangulation3::Edge                    Edge;
    typedef typename Triangulation3::Facet                   Facet;
    typedef typename Triangulation3::Facet_circulator        FacetCirculator;
    typedef typename Triangulation3::Cell_handle             CellHandle;
    typedef typename Triangulation3::Point                   Point;
    typedef typename Triangulation3::Triangle_3              Triangle;
    typedef typename Kernel3::Plane_3                        Plane;
    typedef typename Kernel3::Vector_3                       Vector;
    typedef typename Kernel3::FT                             Coordinate;
    typedef typename Kernel3::FT                             Component;


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~Triangulation3DHelper();

    /**
       @param T any 3D triangulation.
     */
    Triangulation3DHelper( ConstAlias<Triangulation3> T )
      : _T( T )
    {}

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    Triangulation3DHelper ( const Triangulation3DHelper & other )
      : _T( other._T )
    {}

    /// @return a reference to the triangulation.
    inline const Triangulation3 & T() const
    { return _T; }

    /**
       @return the cross product between \a a and \a b.
    */
    inline Vector cross( const Vector & a, const Vector & b ) const
    {
      return Vector( a.y() * b.z() - a.z() * b.y(),
		     a.z() * b.x() - a.x() * b.z(),
		     a.x() * b.y() - a.y() * b.x() );
    }

    /**
       @return the index l, such that {i,j,k,l} = {0,1,2,3}
     */
    int complementIn0123( int i, int j, int k ) const
    {
      ASSERT( (0 <= i) && (i <= 3) );
      ASSERT( (0 <= j) && (j <= 3) );
      ASSERT( (0 <= k) && (k <= 3) );
      return 6 - (i+j+k);
    }

    std::pair<int,int> complementIn0123( int i, int j ) const
    {
      switch ( i ) {
      case 0: switch( j ) {
	case 1: return std::make_pair( 2, 3 );
	case 2: return std::make_pair( 1, 3 );
	case 3: return std::make_pair( 1, 2 );
	default: return std::make_pair( -1, -1 );
	}
      case 1: switch( j ) {
	case 0: return std::make_pair( 2, 3 );
	case 2: return std::make_pair( 0, 3 );
	case 3: return std::make_pair( 0, 2 );
	default: return std::make_pair( -1, -1 );
	}
      case 2: switch( j ) {
	case 0: return std::make_pair( 1, 3 );
	case 1: return std::make_pair( 0, 3 );
	case 3: return std::make_pair( 0, 1 );
	default: return std::make_pair( -1, -1 );
	}
      case 3: switch( j ) {
	case 0: return std::make_pair( 1, 2 );
	case 1: return std::make_pair( 0, 2 );
	case 2: return std::make_pair( 0, 1 );
	default: return std::make_pair( -1, -1 );
	}
      }
      return std::make_pair( -1, -1 ); 
    }

    /**
       Given an edge \a e and a facet \f such that \a e is an edge of
       \a f, and \a e and \a f are defined within the same cell,
       returns the index in the cell of the only vertex of \a f that
       is not in \a e.

       @param e an edge (c,i,j)
       @param f a facet (c,k), k != i, k != j
       @return the index l, such that {i,j,k,l} = {0,1,2,3}
    */
    inline int indexThirdVertex( const Edge & e, const Facet & f ) const
    {
      ASSERT( e.first == f.first );
      return complementIn0123( e.second, e.third, f.second );
    }

    inline VertexHandle thirdVertex( const Edge & e, const Facet & f ) const
    {
      return f.first.vertex( indexThirdVertex( e, f ) );
    }
    
    /**
       Given an edge \a e and a facet \f such that \a e is an edge of
       \a f, and \a e and \a f are defined within the same cell,
       returns the next facet incident to \a e and opposite to the facet.

       @param[inout] e an edge (c,i,j)
       @param[inout] f a facet (c,k), k != i, k != j
    */
    void nextAroundEdge( Edge & e, Facet & f ) const
    {
      ASSERT( e.first == f.first );
      int l = indexThirdVertex( e, f );
      CellHandle c = e.first;
      CellHandle nc = c.neighbor( l ); // next cell
      f.first = nc;
      f.second = T().mirror_index( c, l );
      int i = nc.index( c.vertex( i ) );
      int j = nc.index( c.vertex( j ) );
      e.first = nc;
      e.second = i;
      e.third = j;
    }

    /**
       Given an edge \a e and a facet \f such that \a e is an edge of
       \a f, and \a e and \a f are defined within the same cell,
       returns the previous facet incident to \a e and opposite to the facet.

       @param[inout] e an edge (c,i,j)
       @param[inout] f a facet (c,k), k != i, k != j
    */
    void previousAroundEdge( Edge & e, Facet & f ) const
    {
      ASSERT( e.first == f.first );
      int l = f.second;
      CellHandle c = e.first;
      CellHandle pc = c.neighbor( l ); // previous cell
      f.first = pc;
      f.second = T().mirror_index( c, l );
      int i = pc.index( c.vertex( i ) );
      int j = pc.index( c.vertex( j ) );
      e.first = pc;
      e.second = i;
      e.third = j;
    }

    /**
       @return the inner angle at edge \a e=(c,i,j) within its cell c.
    */
    double innerAngle( const Edge & e ) const
    {
      std::pair<int,int> p = complementIn0123( e.second, e.third );
      Triangle t1 = T().triangle( Facet( e.first, p.first ) );
      Triangle t2 = T().triangle( Facet( e.first, p.second ) );
      Vector n1 = t1.supporting_plane().orthogonal_vector();
      Vector n2 = t2.supporting_plane().orthogonal_vector();
      Component s = ( n1 * n2 ) / sqrt( n1.squared_length() * n2.squared_length() );
      s = ( s < 0.0 ) ? -s : s;
      return acos( s );
    }

    /**
       @return 'true' iff the pentahedron formed by the two tetrahedra
       sharing facet \a f is convex.
    */
    bool isConvexFacet( const Facet & f ) const
    {
      ASSERT( ! T().is_infinite( f.first ) );
      ASSERT( ! T().is_infinite( T().mirror_facet( f ).first ) );
      Point e0 = f.first.vertex( f.second )->point();       // top
      Point f0 = f.first.vertex( (f.second+1)%4 )->point(); // facet p
      Point f1 = f.first.vertex( (f.second+2)%4 )->point(); // facet q
      Point f2 = f.first.vertex( (f.second+3)%4 )->point(); // facet r, (p,q,r) ccw from apex
      Point e1 = T().mirror_vertex( f )->point();           // bottom
      ASSERT( Plane( f0, f1, f2 ).has_on_positive_side( e0 ) );
      ASSERT( Plane( f0, f1, f2 ).has_on_negative_side( e1 ) );
      return Plane( f0, f1, e1 ).has_on_positive_side( e0 )
	&& Plane( f1, f2, e1 ).has_on_positive_side( e0 )
	&& Plane( f2, f0, e1 ).has_on_positive_side( e0 )
	&& Plane( f1, f0, e0 ).has_on_positive_side( e1 )
	&& Plane( f2, f1, e0 ).has_on_positive_side( e1 )
	&& Plane( f0, f2, e0 ).has_on_positive_side( e1 );
    } 

    /**
       \a e forms an edge of the triangulation that is surrounded
       by an umbrella of size 3. The vertices of this umbrella are
       v0, v1, v2.

       @return 'true' iff the pentahedron formed by the three
       tetrahedra sharing edge \a e is convex.
    */
    bool isConvexEdge3( const Edge & e, 
			VertexHandle v0, VertexHandle v1, VertexHandle v2 ) const
    {
      Plane P( v0->point(), v1->point(), v2->point() );
      Point e0 = e.first->vertex( e.second )->point();
      Point e1 = e.first->vertex( e.third )->point();
      return ( P.has_on_positive_side( e0 ) && P.has_on_negative_side( e1 ) )
	|| ( P.has_on_positive_side( e1 ) && P.has_on_negative_side( e0 ) );
    }

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
  protected:
    
    /// The referenced triangulation.
    const Triangulation3 & _T;

    // ------------------------- Private Datas --------------------------------
  private:

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    Triangulation3DHelper();

  private:

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    Triangulation3DHelper & operator= ( const Triangulation3DHelper & other );

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class Triangulation3DHelper


  /**
   * Overloads 'operator<<' for displaying objects of class 'Triangulation3DHelper'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'Triangulation3DHelper' to write.
   * @return the output stream after the writing.
   */
  template <typename TTriangulation3, typename TKernel3>
  std::ostream&
  operator<< ( std::ostream & out, 
	       const Triangulation3DHelper<TTriangulation3,TKernel3> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "Triangulation3DHelper.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined Triangulation3DHelper_h

#undef Triangulation3DHelper_RECURSES
#endif // else defined(Triangulation3DHelper_RECURSES)
