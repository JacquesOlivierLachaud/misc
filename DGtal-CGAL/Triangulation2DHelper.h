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
 * @file Triangulation2DHelper.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2013/06/27
 *
 * Header file for module Triangulation2DHelper.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(Triangulation2DHelper_RECURSES)
#error Recursive header files inclusion detected in Triangulation2DHelper.h
#else // defined(Triangulation2DHelper_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Triangulation2DHelper_RECURSES

#if !defined Triangulation2DHelper_h
/** Prevents repeated inclusion of headers. */
#define Triangulation2DHelper_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/base/ConstAlias.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class Triangulation2DHelper
  /**
     Description of template class 'Triangulation2DHelper' <p>
     \brief Aim: A class that gives a few methods for moving in a Triangulation3.
     
     @tparam TTriangulation2 any kind of CGAL 2D triangulation.
     @tparam TKernel2 any kind of CGAL 2D Kernel.
  */
  template <typename TTriangulation2, typename TKernel2>
  class Triangulation2DHelper
  {
    // ----------------------- public types ------------------------------
  public:
    typedef TTriangulation2                                 Triangulation2;
    typedef TKernel2                                        Kernel2;
    typedef Triangulation2DHelper<Triangulation2, Kernel2>  Self;
    typedef typename Triangulation2::Vertex_circulator      VertexCirculator;
    typedef typename Triangulation2::Vertex_handle          VertexHandle;
    typedef typename Triangulation2::Edge_iterator          EdgeIterator;
    typedef typename Triangulation2::Edge_circulator        EdgeCirculator;
    typedef typename Triangulation2::Edge                   Edge;
    typedef typename Triangulation2::Finite_faces_iterator  FiniteFacesIterator;
    typedef typename Triangulation2::Face_handle            FaceHandle;
    typedef typename Triangulation2::Point                  Point;
    typedef typename Kernel2::Vector_2                      Vector;
    typedef typename Kernel2::FT                            Coordinate;
    typedef typename Kernel2::FT                            Component;
    

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~Triangulation2DHelper();

    /**
       @param T any 3D triangulation.
     */
    Triangulation2DHelper( ConstAlias<Triangulation2> T )
      : myT( T )
    {}

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    Triangulation2DHelper ( const Triangulation2DHelper & other )
      : myT( other.myT )
    {}

    /// @return a reference to the triangulation.
    inline const Triangulation2 & T() const
    { return myT; }

    /**
       Computes the wedge of \a u and \v, otherwise said their determinant.
       @param u any 2d vector.
       @param v any 2d vector.
       @return the scalar det(u,v) or u \wedge v
    */
    inline
    Component det( const Vector & u, const Vector & v ) const
    {
      return u.x() * v.y() - u.y() * v.x();
    }

    /**
       @param e an oriented edge [st] within a ccw face [stu].
       @return the (positive) angle of [st] and [su].
    */
    inline 
    Component innerAngle( const Edge & e ) const
    {
      Vector v1 = target( e )->point() - source( e )->point();
      Vector v2 = ( e.first->vertex( e.second ) )->point() - source( e )->point();
      return acos( v1 * v2 / sqrt( v1.squared_length() ) / sqrt( v2.squared_length() ) );
    }

    /**
       @param e any edge of a ccw face.
       @return the handle of the source vertex of the edge.
    */
    VertexHandle source( const Edge & e ) const
    {
      return e.first->vertex( T().ccw( e.second ) );
    }
    
    /**
       @param e any edge of a ccw face.
       @return the handle of the source target of the edge.
    */
    VertexHandle target( const Edge & e ) const
    {
      return e.first->vertex( T().cw( e.second ) );
    }
    
    /**
       @param e any edge, say [st], of a ccw face [stu].
       @return the next ccw edge within this face, i.e. [tu].
    */
    Edge nextCCWAroundFace( const Edge & e ) const
    {
      return Edge( e.first, T().ccw( e.second ) );
    }
  
    /**
       @param e any edge, say [st], of a ccw face [stu].
       @return the next cw edge within this face, i.e. [us].
    */
    Edge nextCWAroundFace( const Edge & e ) const
    {
      return Edge( e.first, T().cw( e.second ) );
    }
    
    /**
       @param e any edge.

       @return the next edge around the source vertex, such that this
       edge has the same source and is CCW around it.
    */
    Edge nextCCWAroundSourceVertex( const Edge & e ) const
    {
      Edge e1 = nextCWAroundFace( e );
      Edge n = T().mirror_edge( e1 );
      ASSERT( ( source( n ) == source( e ) )
              && "[Triangulation2DHelper::nextCCWAroundSourceVertex] sources are not consistent." );
      ASSERT( ( target( n ) != target( e ) )
              && "[Triangulation2DHelper::nextCCWAroundSourceVertex] targets are equal." );
      return n;
    }
  
    /**
       @param e any edge.

       @return the next edge around the source vertex, such that this
       edge has the same source and is CW around it.
    */
    Edge nextCWAroundSourceVertex( const Edge & e ) const
    {
      return nextCCWAroundFace( T().mirror_edge( e ) );
    }

    /**
       The vertices are ordered ccw around some fictive point from v1 to v4.
       @param v1 a vertex handle
       @param v2 a vertex handle
       @param v3 a vertex handle
       @param v4 a vertex handle
       @return 'true' if the simple polygon [v1,v2,v3,v4] is a convex quadrilateral.
    */
    bool isConvexQuadrilateral( const VertexHandle & v1, 
                                const VertexHandle & v2,
                                const VertexHandle & v3,
                                const VertexHandle & v4 ) const
    {
      Point a( v1->point() );
      Point b( v2->point() );
      Point c( v3->point() );
      Point d( v4->point() );
      return ( det( ( b-a ), ( c-b ) ) > 0 )
        && ( det( ( c-b ), ( d-c ) ) > 0 )
        && ( det( ( d-c ), ( a-d ) ) > 0 )
        && ( det( ( a-d ), ( b-a ) ) > 0 );
    }

    /**
       @param e any edge
       @param v1 any vertex handle
       @param v2 any vertex handle
       @return 'true' iff the edge e is stabbed by the segment [v1v2],
       assuming the following point locations.
       .....x......
       .....|e.....
       .v1..|......
       .x...|...v2.
       .....|...x..
       .....|......
       .....v......
    */
    bool isStabbedCCW( const Edge & e, VertexHandle v1, VertexHandle v2 ) const
    {
      return isConvexQuadrilateral( source( e ), v1, target( e ), v2 );
    }
    
    /**
       @param e any edge
       @param v1 any vertex handle
       @param v2 any vertex handle
       @return 'true' iff the edge e is stabbed by the segment [v1v2],
       without assumption on point locations.
       .....x......
       .....|e.....
       .v1..|......
       .x...|...v2.
       .....|...x..
       .....|......
       .....x......
    */
    bool isStabbed( const Edge & e, VertexHandle v1, VertexHandle v2 ) const
    {
      return isConvexQuadrilateral( source( e ), v1, target( e ), v2 )
        || isConvexQuadrilateral( source( e ), v2, target( e ), v1 );
    }
    
    /**
       @param v1 any vertex handle.
       @param v2 any vertex handle. 
       @return 'true' iff [v1,v2] is an edge of the triangulation. In
       this case, \a e is this edge as an out parameter.
       @pre v1 and v2 are vertex of this triangulation.
    */
    bool findEdge( Edge & e, VertexHandle v1, VertexHandle v2 ) const
    {
      EdgeCirculator ci = T().incident_edges( v1 );
      if ( ci != 0 ) {
        e = *ci;
        if ( source( e ) != v1 ) e = T().mirror_edge( e );
        ASSERT( ( source( e ) == v1 )
                && "[TriangulationHelper::isEdge] Invalid incident edge." );
        Edge start = e;
        bool found = false;
        do {
          if ( target( e ) == v2 ) return true;
          e = nextCCWAroundSourceVertex( e );
        } while ( e != start );
        return false;
      }
      else return false;
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
    
    /// The referenced 2D triangulation.
    const Triangulation2 & myT;

    // ------------------------- Private Datas --------------------------------
  private:

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    Triangulation2DHelper();

  private:

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    Triangulation2DHelper & operator= ( const Triangulation2DHelper & other );

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class Triangulation2DHelper


  /**
   * Overloads 'operator<<' for displaying objects of class 'Triangulation2DHelper'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'Triangulation2DHelper' to write.
   * @return the output stream after the writing.
   */
  template <typename TTriangulation2, typename TKernel2>
  std::ostream&
  operator<< ( std::ostream & out, 
	       const Triangulation2DHelper<TTriangulation2,TKernel2> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "Triangulation2DHelper.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined Triangulation2DHelper_h

#undef Triangulation2DHelper_RECURSES
#endif // else defined(Triangulation2DHelper_RECURSES)
