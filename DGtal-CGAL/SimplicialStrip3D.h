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
 * @file SimplicialStrip3D.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2013/06/27
 *
 * Header file for module SimplicialStrip3D.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(SimplicialStrip3D_RECURSES)
#error Recursive header files inclusion detected in SimplicialStrip3D.h
#else // defined(SimplicialStrip3D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define SimplicialStrip3D_RECURSES

#if !defined SimplicialStrip3D_h
/** Prevents repeated inclusion of headers. */
#define SimplicialStrip3D_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "Triangulation3DHelper.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class SimplicialStrip3D
  /**
     Description of template class 'SimplicialStrip3D' <p> \brief Aim:
     This class represents a sequence of face-adjacent tetrahedra
     around a given initial facet (\a s, \a t, \a u). The user must
     also provide a predicate \a P. Each tetrahedra of the strip must
     be incident to edge (\a s, \a t ). Furthermore, inside facets (\a
     s, \a t, \a u_i) of the strip satisfy \f$ P(u_i) \f$.

     @code
     Triangulation T;
     ...
     Edge e;
     SimplicialStrip<Triangulation,K> strip( T );
     strip.init( P, e ); // P is a predicate VertexHandle -> bool.
     std::cout << strip.angle() << std::endl; // angle is the sum of the inner angle of each face.
     @endcode

     @tparam TTriangulation3 any kind of CGAL 3D triangulation.
     @tparam TKernel3 any kind of CGAL 3D Kernel.
*/
  template <typename TTriangulation3, typename TKernel3>
  class SimplicialStrip3D
  {
    // ----------------------- Standard services ------------------------------
  public:
    typedef TTriangulation3                                  Triangulation3;
    typedef TKernel3                                         Kernel3;
    typedef SimplicialStrip3D<Triangulation3, Kernel3>       Self;
    typedef typename Triangulation3::Vertex_handle           VertexHandle;
    typedef typename Triangulation3::Edge                    Edge;
    typedef typename Triangulation3::Facet                   Facet;
    typedef typename Triangulation3::Facet_circulator        FacetCirculator;
    typedef typename Triangulation3::Cell_handle             CellHandle;
    typedef typename Triangulation3::Point                   Point;
    typedef typename Kernel3::Vector_3                       Vector;
    typedef typename Kernel3::FT                             Coordinate;
    typedef typename Kernel3::FT                             Component;

    static const double EPSILON;

  protected:
    /// The sequence of facets (s,t,u_i) (be careful, the index of the
    /// facet is the remaining index in \f$ \{0,1,2,3\} \setminus
    /// \{s,t,u_i\}.
    std::deque<Facet> _umbrella;
    std::deque<Edge> _edge_umbrella;
    /// True iff the umbrella contains all incident facets to the pivot.
    bool _loop;
    /// Used for computations.
    Triangulation3DHelper<Triangulation3, Kernel3> TH;

  public:
    /// Constructor from triangulation.
    SimplicialStrip3D( ConstAlias<Triangulation3> T ) : TH( T ) {};

    /// Define a strip around the source edge of \a border, starting
    /// around \a border.  The \a predicate should return true for all
    /// vertices inside the strip.  It also checks if faces have been
    /// extended in the Digital Affine Complex \a dac.
    template <typename DAC, typename VertexHandlePredicate>
    void init( const DAC & dac, 
	       const VertexHandlePredicate & predicate, 
	       const Edge & pivot,
	       const Facet & border )
    {
      _edge_umbrella.clear();
      _umbrella.clear();
      ASSERT( pivot.first == border.first );
      ASSERT( pivot.second != border.second );
      ASSERT( pivot.third != border.second );
      ASSERT( ( ! predicate( pivot.first->vertex( pivot.second ) ) )
	      || ( ! predicate( pivot.first->vertex( pivot.third ) ) ) );
      ASSERT( predicate( TH.thirdVertex( pivot, border ) ) );

      _loop = false;
      Edge e = pivot;
      Facet f = border;
      while ( (! _loop) && predicate( TH.thirdVertex( e, f ) )
	      && ( ! dac.isCellExtended( f.first->neighbor( f.second ) ) ) )
      {
	// DGtal::trace.info() << "[SimplicialStrip3D::init] L1 Facet"
	// 		    << " [(" << f.first->vertex( (f.second+1)%4 )->point() << ")"
	// 		    << ",(" << f.first->vertex( (f.second+2)%4 )->point() << ")"
	// 		    << ",(" << f.first->vertex( (f.second+3)%4 )->point() << ")]" 
	// 		    << std::endl;
	TH.previousAroundEdge( e, f );
	// DGtal::trace.info() << "[SimplicialStrip3D::init] L1 prevF"
	// 		    << " [(" << f.first->vertex( (f.second+1)%4 )->point() << ")"
	// 		    << ",(" << f.first->vertex( (f.second+2)%4 )->point() << ")"
	// 		    << ",(" << f.first->vertex( (f.second+3)%4 )->point() << ")]" 
	// 		    << std::endl;
	if ( f == border ) _loop = true; 
	ASSERT( e.first->vertex( e.second ) == pivot.first->vertex( pivot.second ) );
	ASSERT( e.first->vertex( e.third ) == pivot.first->vertex( pivot.third ) );
      }
      if ( _loop ) return; // umbrella is the whole neighborhood.
      while ( ! dac.isCellExtended( f.first ) ) // cell has already been extended.
	{
	  // DGtal::trace.info() << "[SimplicialStrip3D::init] L2 Facet"
	  // 		      << " [(" << f.first->vertex( (f.second+1)%4 )->point() << ")"
	  // 		      << ",(" << f.first->vertex( (f.second+2)%4 )->point() << ")"
	  // 		      << ",(" << f.first->vertex( (f.second+3)%4 )->point() << ")]" 
	  // 		      << std::endl;
	  _edge_umbrella.push_back( e );
	  _umbrella.push_back( f );
	  TH.nextAroundEdge( e, f );
	  ASSERT( e.first->vertex( e.second ) == pivot.first->vertex( pivot.second ) );
	  ASSERT( e.first->vertex( e.third ) == pivot.first->vertex( pivot.third ) );
	  if ( ! predicate( TH.thirdVertex( e, f ) ) ) break;
      }
      _edge_umbrella.push_back( e );
      _umbrella.push_back( f );
      // DGtal::trace.info() << *this;
    }
    
    /// @return 'true' iff the strip was initialized. It must have an
    /// umbrella of size at least 2.
    inline bool isValid() const
    { return ( _umbrella.size() >= 1 ) || isLoop(); }
    
    /// @return 'true' iff the strip contains all incident edges to the
    /// pivot (ie the whole umbrella).
    inline bool isLoop() const 
    { return _loop; }
    
    /// @return 'true' iff the umbrella is reduced to the initial edge.
    inline bool isTrivial() const
    {
      ASSERT( isValid() );
      return _umbrella.size() == 1;
    }
    
    /**
       @param[in] an index between 0 and umbrella.size()-1.
       @return the i-th facet of umbrella.
    */
    inline Facet f( unsigned int i ) const
    {
      ASSERT( isValid() && ( i < _umbrella.size() ) );
      return _umbrella.at( i );
    }

    /**
       @param[in] an index between 0 and size()-1.
       @return the i-th pivot of the umbrella (the same edge, but
       defined with the same cell as f(i)).
    */
    inline Edge pivot( unsigned int i ) const
    {
      ASSERT( isValid() && ( i < _umbrella.size() ) );
      return _edge_umbrella.at( i );
    }

    /**
       @param[in] an index between 0 and umbrella.size()-1.
       @return the \a i-th vertex of umbrella.
    */
    inline VertexHandle v( unsigned int i ) const
    {
      ASSERT( isValid() && ( i < _umbrella.size() ) );
      return TH.thirdVertex( pivot( i ), f( i ) );
    }


    /**
       @return the first vertex of the pivot.
    */
    inline VertexHandle pivotV0() const
    {
      ASSERT( isValid() );
      Edge e = pivot( 0 );
      return e.first->vertex( e.second );
    }

    /**
       @return the first vertex of the pivot.
    */
    inline VertexHandle pivotV1() const
    {
      ASSERT( isValid() );
      Edge e = pivot( 0 );
      return e.first->vertex( e.third );
    }

    inline int indexV( unsigned int i ) const
    {
      ASSERT( isValid() && ( i < _umbrella.size() ) );
      return TH.indexThirdVertex( pivot( i ), f( i ) );
    }

    /// first vertex of umbrella, equivalent to v( 0 )
    inline VertexHandle v0() const
    {
      ASSERT( isValid() );
      return v( 0 );
    }
    /// last vertex of umbrella, equivalent to v( size()-1 )
    inline VertexHandle vn() const
    {
      ASSERT( isValid() );
      return v( size() - 1 );
    }
    /// Number of edges/vertices in umbrella.
    inline unsigned int size() const
    { 
      return _umbrella.size();
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
      if ( isLoop() ) return 2.0*M_PI;
      if ( isTrivial() ) return 2.0*M_PI;
      // if ( ( size() == 2 ) && ( umbrella.front() == umbrella.back() ) )
      //   return 2.0*M_PI;
      Component totalAngle = 0.0;
      for ( typename std::deque<Edge>::const_iterator it = _edge_umbrella.begin(),
	      ite = _edge_umbrella.end() - 1;
	    it != ite; ++it )
	{
	  totalAngle += TH.innerAngle( *it );
	}
      return totalAngle;
    }
    /// Concavity test.
    bool isConcave() const
    {
      ASSERT( isValid() );
      if ( isTrivial() ) return false;
      return ( angle() < M_PI - EPSILON );
    }
    /// Concavity test.
    bool isConvex() const
    {
      ASSERT( isValid() );
      if ( isTrivial() ) return false;
      return ( angle() > M_PI + EPSILON );
    }
    /// Concavity test.
    bool isFlat() const
    {
      ASSERT( isValid() );
      if ( isTrivial() ) return false;
      double a = angle();
      return ( a >= M_PI - EPSILON ) && ( a <= M_PI + EPSILON );
    }

    /**
       @pre this->isConcave()
       @return 'true' iff this (concave) strip is a tetrahedron that
       has already all its vertices in the hull.
    */
    bool checkSimpleCase() const
    {
      ASSERT( isConcave() );
      return size() == 2;
    }
    /**
       @pre this->isConcave()
       @param[out] e returns when 'true' the edge that can be flipped.
       @return 'true' iff this is the special case of flip 3->2 (3 tetrahedra -> 2 tetrahedra).
    */
    bool checkFlippableEdge( Edge & e ) const
    {
      if ( size() == 3 )
	{
	  CellHandle c0 = pivot( 0 ).first->neighbor( pivot( 0 ).second );
	  CellHandle c1 = pivot( 1 ).first->neighbor( pivot( 1 ).second );
	  if ( c0 == c1 )
	    {
	      e = Edge( pivot( 1 ).first, pivot( 1 ).third, indexV( 1 ) );
	      return TH.isConvexEdge3( e, pivotV0(), v( 0 ), v( 2 ) );
	    }
	  c0 = pivot( 0 ).first->neighbor( pivot( 0 ).third );
	  c1 = pivot( 1 ).first->neighbor( pivot( 1 ).third );
	  if ( c0 == c1 )
	    {
	      e = Edge( pivot( 1 ).first, pivot( 1 ).second, indexV( 1 ) );
	      return TH.isConvexEdge3( e, pivotV1(), v( 2 ), v( 0 ) );
	    }
	}
      return false;
    }
    // @return the index of a facet that is 2->3 flippable in the
    // strip.(2 tetrahedra -> 3 tetrahedra).
    unsigned int getFlippableFacetIndex() const
    {
      for ( unsigned int i = 1; i < size() - 1; ++i )
	if ( TH.isConvexFacet( f( i ) ) )
	  return i;
      return size();
    } 

    /**
     * Destructor.
     */
    ~SimplicialStrip3D();

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    // ------------------------- Protected Datas ------------------------------
  private:
    // ------------------------- Private Datas --------------------------------
  private:

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    SimplicialStrip3D();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    SimplicialStrip3D ( const SimplicialStrip3D & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    SimplicialStrip3D & operator= ( const SimplicialStrip3D & other );

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class SimplicialStrip3D


  /**
   * Overloads 'operator<<' for displaying objects of class 'SimplicialStrip3D'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'SimplicialStrip3D' to write.
   * @return the output stream after the writing.
   */
  template <typename TTriangulation3, typename TKernel3>
  std::ostream&
  operator<< ( std::ostream & out, 
	       const SimplicialStrip3D<TTriangulation3, TKernel3> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "SimplicialStrip3D.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined SimplicialStrip3D_h

#undef SimplicialStrip3D_RECURSES
#endif // else defined(SimplicialStrip3D_RECURSES)
