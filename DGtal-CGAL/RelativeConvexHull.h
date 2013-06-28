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
 * @file RelativeConvexHull.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2013/06/27
 *
 * Header file for module RelativeConvexHull.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(RelativeConvexHull_RECURSES)
#error Recursive header files inclusion detected in RelativeConvexHull.h
#else // defined(RelativeConvexHull_RECURSES)
/** Prevents recursive inclusion of headers. */
#define RelativeConvexHull_RECURSES

#if !defined RelativeConvexHull_h
/** Prevents repeated inclusion of headers. */
#define RelativeConvexHull_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "Triangulation3DHelper.h"
#include "SimplicialStrip3D.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class RelativeConvexHull
  /**
     Description of template class 'RelativeConvexHull' <p> \brief
     Aim: This class represents the convex hull of a given set of a
     triangulation of inside points P relatively to a
     (non-intersecting) triangulation of outside points Q.

     @tparam TTriangulation3 any kind of CGAL 3D triangulation.
     @tparam TKernel3 any kind of CGAL 3D Kernel.
  */
  template <typename TTriangulation3, typename TKernel3>
  class RelativeConvexHull
  {
    // ----------------------- public types ------------------------------
  public:
    typedef TTriangulation3                                  Triangulation3;
    typedef TKernel3                                         Kernel3;
    typedef RelativeConvexHull<Triangulation3, Kernel3>      Self;
    typedef CGAL::Delaunay_triangulation_3<Kernel3>          DelaunayTriangulation3;
    typedef typename Triangulation3::Vertex_handle           VertexHandle;
    typedef typename Triangulation3::Edge                    Edge;
    typedef typename Triangulation3::Facet                   Facet;
    typedef typename Triangulation3::Cell_handle             CellHandle;
    typedef typename Triangulation3::Point                   Point;
    typedef typename Triangulation3::Finite_facets_iterator  FiniteFacetsIterator;
    typedef typename Triangulation3::Triangle                Triangle;
    typedef typename Kernel3::Vector_3                       Vector;
    typedef typename Kernel3::Line_3                         Line;
    typedef typename Kernel3::Plane_3                        Plane;
    typedef typename Kernel3::FT                             Coordinate;
    typedef typename Kernel3::FT                             Component;
    typedef int                                              Label;
    typedef typename std::map<VertexHandle,Label>      VertexLabeling;
    typedef Triangulation3DHelper<Triangulation3,Kernel3>    TriangulationHelper;
    typedef SimplicialStrip3D<Triangulation3,Kernel3>        Strip;

    /// @todo. Add extension and retraction mappings.
    struct Extension {
      CellHandle _c;
      Extension( const CellHandle & c )
        : _c( c ) {}
      Extension( const Extension & other )
        : _c( other._c ) {}
      Extension() {}
      Extension& operator=( const Extension & other )
      {
        _c = other._c;
        return *this;
      }
    };

    struct BorderFacet {
      std::pair<VertexHandle, VertexHandle> first; //< pivot edge (both of label l).
      VertexHandle second;                         //< external vertex (of label != l ).

      inline
      BorderFacet( const std::pair<VertexHandle, VertexHandle> & pEdge, 
                   VertexHandle extVertex )
      {
        if ( pEdge.first < pEdge.second ) first = pEdge;
        else first = std::make_pair( pEdge.second, pEdge.first );
        second = extVertex;
      }

      inline
      BorderFacet( VertexHandle pivotVertex0, VertexHandle pivotVertex1,
                   VertexHandle extVertex )
      {
        if ( pivotVertex0 < pivotVertex1 ) 
          first = std::make_pair( pivotVertex0, pivotVertex1 );
        else
          first = std::make_pair( pivotVertex1, pivotVertex0 );
        second = extVertex;
      }

      inline 
      bool operator<( const BorderFacet & other ) const
      {
        return ( first < other.first )
          || ( ( first == other.first ) && ( second < other.second ) );
      }

    };

    typedef typename std::map<CellHandle, Extension> CellMapping;
    typedef typename CellMapping::iterator CellMappingIterator;
    typedef typename CellMapping::const_iterator CellMappingConstIterator;

    /// A predicate that returns 'true' whenever the labeling is not the
    /// one given at instanciation.
    struct CheckVertexLabelingInequality {
      const VertexLabeling & _vl;
      Label _l;
      
      inline
      CheckVertexLabelingInequality( const VertexLabeling & vl, Label l )
        : _vl( vl ), _l( l ) 
      {}
      
      inline
      CheckVertexLabelingInequality( const CheckVertexLabelingInequality & other )
        : _vl( other._vl ), _l( other._l ) 
      {}
      
      inline
      bool operator()( const VertexHandle & v ) const
      {
        typename VertexLabeling::const_iterator it = _vl.find( v );
        return ( it == _vl.end() ) ? false
          : ( it->second != _l );
      }
    };


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~RelativeConvexHull();

    /// Constructor.
    RelativeConvexHull() : TH( _T ) {}

    /// The object is reseted as if just instanciated.
    void clear()
    {
      _T.clear();
      _DT.clear();
      _vLabeling.clear();
      _vMarked.clear();
      _cellsExt.clear();
    }

    /// @return a const reference to the current triangulation.
    inline const Triangulation3 & T() const
    { return _T; }

    /// @return a const reference to the labeling map of vertices.
    inline const VertexLabeling & labeling() const
    { return _vLabeling; }
    
    /**
       @param v any vertex.
       @return the label of this vertex (or INVALID if not labelled).
    */
    inline Label label( VertexHandle v ) const
    {
      typename VertexLabeling::const_iterator it = labeling().find( v );
      if ( it == labeling().end() ) return INVALID;
      return it->second;
    }

    /**
       @param v any vertex.
       @return the mark of this vertex  (or INVALID if not marked).
    */
    inline Label mark( VertexHandle v ) const
    {
      typename VertexLabeling::const_iterator it = _vMarked.find( v );
      if ( it == _vMarked.end() ) return INVALID;
      return it->second;
    }

    /**
       @param v any vertex.
       @return 'true' iff the mark of this vertex is not INVALID.
    */
    inline Label isMarked( VertexHandle v ) const
    {
      return mark( v ) != INVALID;
    }

    /**
       @param c any cell of the triangulation
       @return 'true' iff the cell has been extended during the relative
       hull process and is part of global retraction of the relative
       hull onto the canonic complex.
    */
    inline bool isCellExtended( CellHandle c ) const
    {
      CellMappingConstIterator it = _cellsExt.find( c );
      return it != _cellsExt.end();
    }
  
    /**
       @param c a cell of the triangulation that has been extended during the relative
       hull process.
       @return a const_reference to the associated extension.
       @pre 'isCellExtended( f ) == true'
    */
    const Extension & extension( CellHandle c ) const
    {
      ASSERT( isCellExtended( c ) );
      return *( _cellsExt.find( c ) );
    }
    
    inline void setInfiniteLabel( Label l )
    {
      _vLabeling[ _T.infinite_vertex() ] = l;
    }


    // ----------------------- Insertion services --------------------------------------
  public:
    
    /** 
        Adds points to the triangulation.
        
        @param k the digital topology: 0 is 4-adjacency, 1 is 8. @todo Ignored for now.
        @param l the label for the set of points (0 stands for P, 1 for others).
    */
    template < typename PointFunctor, 
               typename PointIterator >
    void add( PointFunctor f, PointIterator itb, PointIterator ite, int k, Label l )
    {
      for ( ; itb != ite; ++itb )
        {
          Point p = f( *itb );
          VertexHandle v = _DT.insert( p );
          _vLabeling[ v ] = l;
        }
      // addConstraints0( current );
      // if ( k == 1 ) addConstraints1( current );
    }

    /// Terminate the point insertion.
    void terminate() 
    {
      _T.swap( _DT );
    }

    /**
       @return 'true' iff the facet \a f is a boundary facet of label \a l.
    */
    bool isFacetOnBoundary( Facet f, Label l ) const
    {
      CellHandle c = f.first;
      int i = f.second;
      return ( ( label( c->vertex( i ) ) != l )
               && ( label( c->vertex( (i+1)%4 ) ) == l )
               && ( label( c->vertex( (i+2)%4 ) ) == l )
               && ( label( c->vertex( (i+3)%4 ) ) == l ) );
    }

    /**
       @return 'true' iff the facet \a f is a boundary facet of label \a l.
    */
    bool isFacetOnRelativeBoundary( Facet f, Label l ) const
    {
      CellHandle c = f.first;
      int i = f.second;
      return ( ( label( c->vertex( i ) ) != l )
               && ( ( label( c->vertex( (i+1)%4 ) ) == l ) || isMarked( c->vertex( (i+1)%4 ) ) )
               && ( ( label( c->vertex( (i+2)%4 ) ) == l ) || isMarked( c->vertex( (i+2)%4 ) ) )
               && ( ( label( c->vertex( (i+3)%4 ) ) == l ) || isMarked( c->vertex( (i+3)%4 ) ) ) );
    }

    /**
       @return 'true' iff the facet \a f is a boundary facet of label \a l.
    */
    bool isFacetSubdivided( Facet f, Label l ) const
    {
      CellHandle c = f.first;
      int i = f.second;
      return ( ( label( c->vertex( i ) ) != l )
               && ( ( label( c->vertex( (i+1)%4 ) ) == l )
		    || ( label( c->vertex( (i+1)%4 ) ) == 2 ) )
               && ( ( label( c->vertex( (i+2)%4 ) ) == l )
		    || ( label( c->vertex( (i+2)%4 ) ) == 2 ) )
               && ( ( label( c->vertex( (i+3)%4 ) ) == l )
		    || ( label( c->vertex( (i+3)%4 ) ) == 2 ) ) );
    }

    /**
       Used by relativeHull2 to insert candidate concavities into the Queue.
       Check all possibilities for a facet specified as 3 vertices.
    */
    void insertQueue( std::set< BorderFacet > & inQueue,
                      VertexHandle v[ 3 ], 
                      const CheckVertexLabelingInequality & predicate ) const
    {
      CellHandle c;
      int indices[ 3 ];
      if ( T().is_facet( v[ 0 ], v[ 1 ], v[ 2 ], 
                         c, indices[ 0 ], indices[ 1 ], indices[ 2 ] ) )
        {
          bool pred_v[ 3 ];
          bool mark_v[ 3 ];
          for ( unsigned int l = 0; l < 3; ++l ) 
            {
              pred_v[ l ] = predicate( v[ l ] );
              mark_v[ l ] = mark( v[ l ] ) != INVALID;
            }
          for ( unsigned int l = 0; l < 3; ++l ) 
            {
              if ( pred_v[ (l+2)%3 ] 
                   && ( ( ( ! pred_v[ l ] ) && ( ( ! pred_v[ (l+1)%3 ] )
                                                 || mark_v[ (l+1)%3 ] ) )
                        || ( ( ! pred_v[ (l+1)%3 ] ) && ( ( ! pred_v[ l ] )
                                                          || mark_v[ l ] ) ) ) )
                {
                  Edge pivot( c, indices[ l ], indices[ (l+1)%3 ] );
                  Facet border( c, 
                                TH.complementIn0123( indices[ 0 ], indices[ 1 ], indices[ 2 ] ) );
                  Strip strip( T() );
                  strip.init( *this, predicate, pivot, border );
                  if ( strip.isConcave() ) 
                    inQueue.insert( BorderFacet( v[ l ], v[ (l+1)%3 ], v[ (l+2)%3 ] ) );
                }
            }
        }
    }

    /**
       Used by relativeHull2 to insert candidate concavities into the Queue.
       Check all possibilities for a facet specified as 3 vertices.
    */
    void insertQueue( std::set< BorderFacet > & inQueue,
                      const Facet & f,
                      const CheckVertexLabelingInequality & predicate ) const
    {
      int indices[ 3 ];
      VertexHandle v[ 3 ];
      bool pred_v[ 3 ];
      bool mark_v[ 3 ];
      CellHandle c = f.first;
      for ( unsigned int l = 0; l < 3; ++l ) 
        {
          indices[ l ] = (f.second+l+1) % 4;
          v[ l ] = c->vertex( indices[ l ] );
          pred_v[ l ] = predicate( v[ l ] );
          mark_v[ l ] = mark( v[ l ] ) != INVALID;
        }
      for ( unsigned int l = 0; l < 3; ++l ) 
        {
          // if ( ( pred_v[ (l+2)%3 ] && ( ! pred_v[ l ] ) && ( ! pred_v[ (l+1)%3 ] ) ) )
          if ( pred_v[ (l+2)%3 ] 
               && ( ( ( ! pred_v[ l ] ) && ( ( ! pred_v[ (l+1)%3 ] )
                                             || mark_v[ (l+1)%3 ] ) )
                    || ( ( ! pred_v[ (l+1)%3 ] ) && ( ( ! pred_v[ l ] )
                                                      || mark_v[ l ] ) ) ) )
            {
              Edge pivot( c, indices[ l ], indices[ (l+1)%3 ] );
              Facet border( c, 
                            TH.complementIn0123( indices[ 0 ], indices[ 1 ], indices[ 2 ] ) );
              Strip strip( T() );
              strip.init( *this, predicate, pivot, border );
              // DGtal::trace.info() << "- strip s=" << strip.size() << " a=" << strip.angle()
              //                     << std::endl;
              if ( strip.isConcave() ) 
                inQueue.insert( BorderFacet( v[ l ], v[ (l+1)%3 ], v[ (l+2)%3 ] ) );
            }
        }
    }


    /**
       This procedure computes the relative hull with the same principle
       as removeConcavities. around vertices with label \a l. Complexity
       is O(f), where f is the number of faces.
       
       This is a kind of simplest algorithm, where we do not try to sort
       concavities, or push again edges around
       concavities. Surprisingly, it is hard to do better without
       complexifying a lot the way concavities (v0,v1,vext) are pushed into
       the queue.

       @param[in] any label, as given with method \ref add.
       @return 'true' if some concavity was flipped, 'false' when no concavity was found.
    */
    bool relativeHull( Label l ) 
    {
      bool changes = false;
      static const Label rmark = 1;
      std::set< BorderFacet > inQueue;
      CheckVertexLabelingInequality predNotL( labeling(), l );
      Strip strip( T() );
      // DGtal::trace.beginBlock( "Searching concavities" );
      for ( FiniteFacetsIterator it = T().finite_facets_begin(), itend = T().finite_facets_end();
            it != itend; ++it )
        {
          Facet f = *it;
          if ( T().is_infinite( T().mirror_facet( f ) ) ) continue; // both facets should be finite.
          insertQueue( inQueue, f, predNotL );
        }
      DGtal::trace.info() << "- Found " << inQueue.size() << " potential concavities." 
                          << std::endl;
      // DGtal::trace.endBlock();
      // DGtal::trace.beginBlock( "Flipping concavities" );
      unsigned int nb_checked = 0;
      unsigned int nb_flipped3_2 = 0;
      unsigned int nb_flipped2_3 = 0;
      unsigned int nb_simple = 0;
      unsigned int nb_concave = 0;
      unsigned int nb_new = 0;
      Edge pivot;
      Facet border;
      CellHandle c;
      int i, j, k;
      while ( ! inQueue.empty() )
        {
          BorderFacet bf = *( inQueue.begin() );
          VertexHandle v0 = bf.first.first;
          VertexHandle v1 = bf.first.second;
          VertexHandle ve = bf.second;
          ++nb_checked;
          if ( nb_checked % 1000 == 0 ) 
            DGtal::trace.info() << "- Queue=" << inQueue.size()
				<< ", flipped 2->3=" << nb_flipped2_3
				<< ", 3->2=" << nb_flipped3_2
				<< ", ccv=" << nb_concave
				<< ", simple=" << nb_simple
				<< ", new=" << nb_new
				<< " / " << nb_checked << " concave facets." << std::endl;
          if ( ! T().is_facet( v0, v1, ve, c, i, j, k ) )
            {
              inQueue.erase( inQueue.begin() );
              continue; // (v0,v1,ve) is not a facet anymore.
            }
          border = Facet( c, TH.complementIn0123( i, j, k ) );
          pivot = Edge( c, i, j );
          strip.init( *this, predNotL, pivot, border );
          ASSERT( strip.isValid() );
          if ( strip.isConcave() ) 
            {
	      // DGtal::trace.info() << "Ccv strip " << strip << std::endl;
              Edge e;
              unsigned int i;
              if ( strip.checkSimpleCase() )
                {
                  _cellsExt[ strip.f( 0 ).first ] = Extension( strip.f( 0 ).first );
                  inQueue.erase( inQueue.begin() );
		  changes = true; ++nb_simple;
                }
              else if ( strip.checkFlippableEdge( e ) )
                {
                  bool ok = _T.flip( e );
                  if ( ! ok ) 
                    DGtal::trace.error() << "[relativeHull()] Error flipping 3 -> 2."
                                         << std::endl;
                  // the created inside cell will be added at next iteration with checkSimpleCase.
                  inQueue.erase( inQueue.begin() );
		  changes = true; ++nb_flipped3_2;
                }
              else if ( ( i = strip.getFlippableFacetIndex() ) != strip.size() )
                { // at least one facet is flippable.
                  Facet facet = strip.f( i );
                  ASSERT( ! isCellExtended( facet.first ) );
                  Facet mirror_first = T().mirror_facet( strip.f( 0 ) );
                  bool ok = _T.flip( facet );
                  if ( ! ok ) 
                    DGtal::trace.error() << "[relativeHull()] Error flipping 2 -> 3."
                                         << std::endl;
                  if ( strip.size() == 3 ) 
                    { // last flip made a tetrahedron
                      Facet nfacet = T().mirror_facet( mirror_first );
                      _cellsExt[ nfacet.first ] = Extension( nfacet.first );
                      inQueue.erase( inQueue.begin() );
                    }
                  changes = true; ++nb_flipped2_3;
                }
              else // none are flippable. Check is their a flat face.
		{
		  Edge e;
		  if ( ( strip.size() == 3 )
		       && TH.isFlatConvexFacet( strip.f( 1 ), e ) )
		    {
		      Point v0 = strip.v( 2 )->point();
		      Point v1 = strip.v( 0 )->point();
		      Point mid = 
			TH.intersect( T().triangle( strip.f( 1 ) ).supporting_plane(),
				      Line( v0, v1 ) );
		      VertexHandle nv = _T.insert_in_edge( mid, e );
		      // _vMarked[ nv ] = rmark;
		      DGtal::trace.info() << "- Added point " << mid 
					  << " (" << nv->point() << ")" << std::endl;
		      _vLabeling[ nv ] = l;
		      changes = true; ++nb_new;
		      inQueue.erase( inQueue.begin() );
		    }
		  else
		    { // This is a concave piece.
		      for ( unsigned int i = 1; i < strip.size() - 1; ++i ) 
			{
			  if ( TH.isConcaveFacet( strip.pivot( i ), strip.f( i ),
						  strip.v0(), strip.vn() ) 
			       && TH.isConcaveFacet( strip.pivot( i ), strip.f( i ),
						     strip.v( i - 1 ), strip.v( i + 1 ) ) )
			    {
			      Facet facet = strip.f( i - 1 );
			      if ( ! isCellExtended( facet.first ) )
				_cellsExt[ facet.first ] = Extension( facet.first );
			      facet = strip.f( i );
			      if ( ! isCellExtended( facet.first ) )
				_cellsExt[ facet.first ] = Extension( facet.first );
			      if ( predNotL( strip.v( i ) ) )
				_vMarked[ strip.v( i ) ] = rmark;
			      changes = true; ++nb_concave;
			    }
			}
		      inQueue.erase( inQueue.begin() );
		    }
		}
            }
          else
            inQueue.erase( inQueue.begin() );
        }
      ASSERT( inQueue.empty() );
      DGtal::trace.info() << "- Flipped 2->3=" << nb_flipped2_3
			  << ", 3->2=" << nb_flipped3_2
			  << ", ccv=" << nb_concave
			  << ", simple=" << nb_simple
			  << ", new=" << nb_new
			  << " / " << nb_checked << " concave facets." << std::endl;
      //if ( nb_flipped == 1 && nb_checked == 2 ) changes = false;
      // DGtal::trace.endBlock();
      return changes;
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
  private:
    /// The initial Delaunauy triangulation.
    DelaunayTriangulation3 _DT;
    /// The current triangulation.
    Triangulation3 _T;
    /// A mapping VertexHandle -> Label that stores for each vertex to which set it belongs.
    VertexLabeling _vLabeling;
    /// A mapping VertexHandle -> Label that stores for each vertex if it is marked.
    VertexLabeling _vMarked;
    /// A mapping Face -> Info that stores faces that have been added to
    /// the relative hull for each vertex if it is marked.
    CellMapping _cellsExt;

    // ------------------------- Private Datas --------------------------------
  private:

  public:
    /// Provides useful methods on triangulation
    Triangulation3DHelper<Triangulation3,Kernel3> TH;
    static const Label INVALID = -1;

    // ------------------------- Hidden services ------------------------------
  protected:

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    RelativeConvexHull ( const RelativeConvexHull & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    RelativeConvexHull & operator= ( const RelativeConvexHull & other );

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class RelativeConvexHull


  /**
   * Overloads 'operator<<' for displaying objects of class 'RelativeConvexHull'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'RelativeConvexHull' to write.
   * @return the output stream after the writing.
   */
  template <typename TTriangulation3, typename TKernel3>
  std::ostream&
  operator<< ( std::ostream & out, 
               const RelativeConvexHull<TTriangulation3,TKernel3> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "RelativeConvexHull.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined RelativeConvexHull_h

#undef RelativeConvexHull_RECURSES
#endif // else defined(RelativeConvexHull_RECURSES)
