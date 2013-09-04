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
 * @file UmbrellaPart2D.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2013/08/24
 *
 * Header file for module UmbrellaPart2D.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(UmbrellaPart2D_RECURSES)
#error Recursive header files inclusion detected in UmbrellaPart2D.h
#else // defined(UmbrellaPart2D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define UmbrellaPart2D_RECURSES

#if !defined UmbrellaPart2D_h
/** Prevents repeated inclusion of headers. */
#define UmbrellaPart2D_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <deque>
#include "DGtal/base/Common.h"
#include "Triangulation2DHelper.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class UmbrellaPart2D
  /**
     Description of template class 'UmbrellaPart2D' <p>
     \brief Aim:
     
     This class represents a sequence of edge-adjacent faces around a
     given vertex \a v: it is thus a part of the umbrella around \a
     v. There are several ways for creating an umbrella part, the most
     common one is to give a starting edge and a stopping predicate.
     
     @code
     Triangulation T;
     ...
     Edge e;
     SimplicialStrip<Triangulation,K> strip( T );
     strip.init( P, e ); // P is a predicate VertexHandle -> bool.
     std::cout << strip.angle() << std::endl; // angle is the sum of the inner angle of each face.
     @endcode
     
     @tparam TTriangulation2 any kind of CGAL 2D triangulation.
     @tparam TKernel2 any kind of CGAL 2D Kernel.
  */
  template <typename TTriangulation2, typename TKernel2>
  class UmbrellaPart2D
  {
    // ----------------------- Public types ------------------------------
  public:
    typedef TTriangulation2                                  Triangulation2;
    typedef TKernel2                                         Kernel2;
    typedef UmbrellaPart2D<Triangulation2, Kernel2>          Self;
    typedef typename Triangulation2::Vertex_circulator       VertexCirculator;
    typedef typename Triangulation2::Vertex_handle           VertexHandle;
    typedef typename Triangulation2::Edge_iterator           EdgeIterator;
    typedef typename Triangulation2::Edge                    Edge;
    typedef typename Triangulation2::Finite_faces_iterator   FiniteFacesIterator;
    typedef typename Triangulation2::Face_handle             FaceHandle;
    typedef typename Triangulation2::Point                   Point;
    typedef typename Kernel2::Vector_2                       Vector;
    typedef typename Kernel2::FT                             Coordinate;
    typedef typename Kernel2::FT                             Component;

    static const double EPSILON;

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~UmbrellaPart2D();

    /**
       Constructor from triangulation..
       @param T any 2D triangulation.
    */
    UmbrellaPart2D( const Triangulation2 & T ) : TH( T ) {};

    /**
       Define a part of the umbrella around the source edge of \a edge, starting
       around \a edge.  The \a predicate should return true for all
       vertices inside the strip.  As a result, when ! myLoop, we have
       ! predicate( target( myUmbrella.front() ) ) and ! predicate(
       target( myUmbrella.back() ) ), while the predicate is true for
       all edges in-between.

       @tparam VertexHandlePredicate the type that defines a predicate
       VertexHandle -> bool.

       @param predicate the concrete predicate instance.

       @param edge the edge whose source defines the pivot of the
       umbrella part, and whose target defines the part of the
       umbrella around the pivot.
    */
    template <typename VertexHandlePredicate>
    void init( const VertexHandlePredicate & predicate, 
               const Edge & edge )
    {
      myUmbrella.clear();
      VertexHandle pivot = TH.source( edge );
      myLoop = false;
      Edge e = edge;
      while ( (! myLoop) && predicate( TH.target( e ) ) )
        {
          myUmbrella.push_back( e );
          e = TH.nextCCWAroundSourceVertex( e );
          if ( e == edge )
            myLoop = true; 
          ASSERT( TH.source( e ) == pivot );
        }
      if ( ! myLoop )
        {
          myUmbrella.push_back( e ); // last
          Edge e = TH.nextCWAroundSourceVertex( edge );
          while ( predicate( TH.target( e ) ) )
            {
              myUmbrella.push_front( e );
              e = TH.nextCWAroundSourceVertex( e );
              ASSERT( TH.source( e ) == pivot );
            }
          myUmbrella.push_front( e );
        }
      checkInfinite();
    }

    /**
       Defines a part of the umbrella around the source edge of \a
       border, starting around \a edge.  The two predicates should
       return true for all edges and faces inside the umbrella part.

       @tparam FaceHandlePredicate the type that defines a predicate
       FaceHandle -> bool.
       @tparam EdgePredicate the type that defines a predicate
       Edge -> bool.

       @param facePredicate the concrete face predicate instance.
       @param edgePredicate the concrete edge predicate instance.

       @param edge the edge whose source defines the pivot of the
       umbrella part, and whose target defines the part of the
       umbrella around the pivot.
     */
    template <typename FaceHandlePredicate, typename EdgePredicate>
    void init( const FaceHandlePredicate & facePredicate,
               const EdgePredicate & edgePredicate, 
               const Edge & edge )
    {
      myUmbrella.clear();
      VertexHandle pivot = TH.source( edge );
      myLoop = false;
      Edge e = edge;
      while ( (! myLoop) && edgePredicate( e )
              && ( facePredicate( TH.T().mirror_edge( e ).first ) ) )
        {
          e = TH.nextCWAroundSourceVertex( e );
          if ( e == edge )
            myLoop = true; 
          ASSERT( TH.source( e ) == pivot );
        }
      if ( myLoop ) return; // umbrella is the whole neighborhood.
      while ( facePredicate( e.first ) ) // face has already been extended.
        {
          myUmbrella.push_back( e );
          e = TH.nextCCWAroundSourceVertex( e );
          ASSERT( TH.source( e ) == pivot );
          if ( ! edgePredicate( e ) ) break;
        }
      myUmbrella.push_back( e );
      // std::cout << "[Umbrella] v=" << pivot->point() << " s=" << umbrella.size()<< std::endl;
      ASSERT( edgePredicate( edge ) );
      checkInfinite();
    }

    /// @return 'true' iff the strip contains all incident edges to the
    /// pivot (ie the whole umbrella).
    inline bool isLoop() const 
    { return myLoop; }

    /// @return 'true' iff the strip contains at least one infinite
    /// edge/face. In this case, method angle() has no meaning.
    inline bool isInfinite() const 
    { return myInfinite; }

    /// @return 'true' iff the umbrella is reduced to the initial edge.
    inline bool isTrivial() const
    {
      ASSERT( isValid() );
      return myUmbrella.size() <= 1;
    }

    /// @return 'true' iff the strip was initialized. It must have an
    /// umbrella of size at least 2.
    inline VertexHandle pivot() const
    { 
      ASSERT( isValid() );
      return TH.source( myUmbrella.front() );
    }

    /**
       @param[in] an index between 0 and myUmbrella.size()-1.
       @return the \a i-th vertex of myUmbrella.
    */
    inline VertexHandle v( unsigned int i ) const
    {
      ASSERT( isValid() && ( i < myUmbrella.size() ) );
      return TH.target( myUmbrella.at( i ) );
    }

    /**
       @param[in] an index between 0 and myUmbrella.size()-1.
       @return the i-th edge of myUmbrella.
    */
    inline Edge e( unsigned int i ) const
    {
      ASSERT( isValid() && ( i < myUmbrella.size() ) );
      return myUmbrella.at( i );
    }

    /**
       @param[in] an index between 0 and myUmbrella.size()-1.
       @return the i-th face handle of myUmbrella.
    */
    inline FaceHandle f( unsigned int i ) const
    {
      ASSERT( isValid() && ( i < myUmbrella.size() ) );
      return e( i ).first;
    }

    /// first edge of umbrella, equivalent to e( 0 )
    inline Edge e0() const
    {
      ASSERT( isValid() );
      return myUmbrella.front();
    }
    /// last edge of umbrella, equivalent to e( size()-1 )
    inline Edge en() const
    {
      ASSERT( isValid() );
      return myUmbrella.back();
    }
    /// first vertex of umbrella, equivalent to v( 0 )
    inline VertexHandle v0() const
    {
      ASSERT( isValid() );
      return TH.target( myUmbrella.front() );
    }
    /// last vertex of umbrella, equivalent to v( size()-1 )
    inline VertexHandle vn() const
    {
      ASSERT( isValid() );
      return TH.target( myUmbrella.back() );
    }
    /// Number of edges/vertices in umbrella (at least 2, which are then equal).
    inline unsigned int size() const
    { 
      return myUmbrella.size();
    }
    /// A priority value for the strip, the smaller, the highest.
    inline double priority() const
    {
      return angle(); // ( vn()->point() - v0()->point() ).squared_length();
    }
    
    /// Angle of umbrella.
    inline Component angle() const
    {
      ASSERT( isValid() );
      ASSERT( ! isInfinite() );
      if ( isLoop() ) return 2.0*M_PI;
      if ( isTrivial() ) return 2.0*M_PI;
      Component totalAngle = 0.0;
      for ( typename std::deque<Edge>::const_iterator it = myUmbrella.begin(), ite = myUmbrella.end() - 1;
            it != ite; ++it )
        {
          totalAngle += TH.innerAngle( *it );
        }
      return totalAngle;
    }
    /// Concavity test.
    bool isConcaveAndFlippable() const
    {
      ASSERT( isValid() );
      if ( isTrivial() || isInfinite() ) return false;
      if ( angle() >= M_PI - EPSILON ) return false;
      bool quadrilaterals = true;
      for ( unsigned int i = 1; i < size() - 1; ++i )
        {
          if ( TH.T().is_constrained( e( i ) )
               || ( ! TH.isStabbedCCW( e( i ), v0(), vn() ) ) )
            {
              quadrilaterals = false;
              break;
            }
        }
      return quadrilaterals;
    }
    /// Concavity test.
    bool isConcave() const
    {
      ASSERT( isValid() );
      if ( isTrivial() || isInfinite() ) return false;
      return ( angle() < M_PI - EPSILON );
    }
    /// Convexity test.
    bool isConvex() const
    {
      ASSERT( isValid() );
      if ( isTrivial() ) return false;
      if ( isInfinite() ) return true;
      return ( angle() > M_PI + EPSILON );
    }
    /// flat test.
    bool isFlat() const
    {
      ASSERT( isValid() );
      if ( isTrivial() || isInfinite() ) return false;
      double a = angle();
      return ( a <= M_PI + EPSILON ) && ( a >= M_PI - EPSILON );
    }
    /// @return 'true' iff the umbrella part is flippable around its
    /// extremeties.
    bool isFlippable() const
    {
      ASSERT( isValid() );
      if ( isTrivial() || isInfinite() ) return false;
      bool quadrilaterals = true;
      for ( unsigned int i = 1; i < size() - 1; ++i )
        {
          if ( TH.T().is_constrained( e( i ) )
               || ( ! TH.isStabbedCCW( e( i ), v0(), vn() ) ) )
            {
              quadrilaterals = false;
              break;
            }
        }
      return quadrilaterals;
    }
    
    // @return the index of an edge that is flippable in the strip.
    unsigned int getFlippableEdgeIndex() const
    {
      for ( unsigned int i = 1; i < size() - 1; ++i )
        if ( TH.isStabbedCCW( e( i ), v( i-1 ), v( i+1 ) ) ) 
          return i;
      return size();
    } 
    
    // @return the index of the concavity in the strip.
    unsigned int getUnflippableEdgeIndex() const
    {
      ASSERT( size() > 2 );
      for ( unsigned int i = 1; i < size() - 1; ++i )
        if ( ! TH.isStabbedCCW( e( i ), 
                                v( i-1 ), 
                                v( i+1 ) ) )
          return i;
      // DGtal::trace.error() << "[SimplicialStrip::getUnflippableEdge] Unable to find a valid edge." << std::endl;
      return size();
    } 

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
       Checks the validity/consistency of the object. return 'true' iff
       the strip was initialized. It must have an umbrella of size at least
       1.  
       @return 'true' if the object is valid, 'false' otherwise.
     */
    inline bool isValid() const
    { return ( myUmbrella.size() >= 1 ) || isLoop(); }

    // ------------------------- Protected Datas ------------------------------
  protected:

    /**
       The sequence of edges (v,v_i). The face of the last edge is not
       included in the umbrella.
    */
    std::deque<Edge> myUmbrella;
    /// True iff the umbrella contains all incident edges to the pivot.
    bool myLoop;
    /// True iff the umbrella contains at least one infinite edge/face.
    bool myInfinite;
    /// Used for computations.
    Triangulation2DHelper<Triangulation2, Kernel2> TH;
    
    // ------------------------- Private Datas --------------------------------
  private:

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    UmbrellaPart2D();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    UmbrellaPart2D ( const UmbrellaPart2D & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    UmbrellaPart2D & operator= ( const UmbrellaPart2D & other );

    // ------------------------- Internals ------------------------------------
  private:

    void checkInfinite()
    {
      myInfinite = false;
      for ( unsigned int i = 0; i < size(); ++i )
        if ( TH.T().is_infinite( myUmbrella[ i ] ) )
          {
            myInfinite = true;
            break;
          }
    }


  }; // end of class UmbrellaPart2D


  /**
   * Overloads 'operator<<' for displaying objects of class 'UmbrellaPart2D'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'UmbrellaPart2D' to write.
   * @return the output stream after the writing.
   */
  template <typename TTriangulation2, typename TKernel2>
  std::ostream&
  operator<< ( std::ostream & out, const UmbrellaPart2D<TTriangulation2, TKernel2> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "UmbrellaPart2D.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined UmbrellaPart2D_h

#undef UmbrellaPart2D_RECURSES
#endif // else defined(UmbrellaPart2D_RECURSES)
