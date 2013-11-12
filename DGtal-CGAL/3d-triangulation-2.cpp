#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <QtGui/qapplication.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>


#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/Shapes.h>
#include <DGtal/shapes/ShapeFactory.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "Auxiliary.h"
#include "Triangulation3DHelper.h"
// #include "SimplicialStrip3D.h"
// #include "RelativeConvexHull.h"

typedef DGtal::SpaceND<3, DGtal::int64_t> Z3;
typedef Z3::Point PointZ3;
typedef Z3::Point RealPoint;
typedef DGtal::HyperRectDomain<Z3> Domain;
typedef Domain::ConstIterator DomainConstIterator;

template <typename TKernel3>
struct toCGALFunctor {
  typedef typename TKernel3::Point_3 Point3;
  typedef typename TKernel3::Vector_3 Vector3;
  inline Point3 operator()( const PointZ3 & p ) const
  {
    return Point3( p[0], p[1], p[2] );
  }
};

template <typename TKernel3>
struct toDGtalFunctor {
  typedef typename TKernel3::Point_3 Point3;
  typedef typename TKernel3::Vector_3 Vector3;
  inline PointZ3 operator()( const Point3 & p ) const
  {
    return PointZ3( p.x(), p.y(), p.z() );
  }
};

inline
PointZ3 operator^( const PointZ3 & a, const PointZ3 & b )
{
  BOOST_STATIC_ASSERT( PointZ3::dimension == 3 );
  return PointZ3( a[ 1 ] * b[ 2 ] - a[ 2 ] * b[ 1 ],
		  a[ 2 ] * b[ 0 ] - a[ 0 ] * b[ 2 ],
		  a[ 0 ] * b[ 1 ] - a[ 1 ] * b[ 0 ] );
}



template <typename Viewer, typename ToDGtal, typename Facet>
void displayFacet( Viewer & viewer, const ToDGtal & toDGtal,
		   const Facet & f, const DGtal::Color & col, 
                   double retract )
{
  typedef typename ToDGtal::Point3 Point;
  typedef typename ToDGtal::Vector3 Vector;
  Point a( f.first->vertex( (f.second+1)%4 )->point() );
  Point b( f.first->vertex( (f.second+2)%4 )->point() );
  Point c( f.first->vertex( (f.second+3)%4 )->point() );
  Point mid = ( Vector( a ) + Vector( b ) + Vector( c ) ) / 3.0;
  a = a + ( mid-a )*retract;
  b = b + ( mid-b )*retract;
  c = c + ( mid-c )*retract;
  viewer.setFillColor( col );
  viewer.addTriangle( RealPoint( a.x(), a.y(), a.z() ),
		      RealPoint( c.x(), c.y(), c.z() ),
		      RealPoint( b.x(), b.y(), b.z() ) );
}

template <typename Viewer, typename ToDGtal>
void displayTriangle( Viewer & viewer, const ToDGtal & toDGtal,
                      typename ToDGtal::Point3 a,
                      typename ToDGtal::Point3 b,
                      typename ToDGtal::Point3 c,
                      const DGtal::Color & col, 
                      double retract )
{
  typedef typename ToDGtal::Point3 Point;
  typedef typename ToDGtal::Vector3 Vector;
  Vector vmid = ( Vector( a.x(), a.y(), a.z() ) + Vector( b.x(), b.y(), b.z() ) + Vector( c.x(), c.y(), c.z() ) ) / 3.0;
  Point mid( vmid.x(), vmid.y(), vmid.z() );
  a = a + ( mid-a )*retract;
  b = b + ( mid-b )*retract;
  c = c + ( mid-c )*retract;
  viewer.setFillColor( col );
  viewer.addTriangle( RealPoint( a.x(), a.y(), a.z() ),
		      RealPoint( c.x(), c.y(), c.z() ),
		      RealPoint( b.x(), b.y(), b.z() ) );
}

template <typename Viewer, typename ToDGtal, typename Cell>
void displayCell( Viewer & viewer, const ToDGtal & toDGtal,
		  const Cell & c, const DGtal::Color & col,
                  double retract )
{
  for ( int i = 0; i < 4; ++i )
    displayFacet( viewer, toDGtal, std::make_pair( c, i ), col, retract );
}



// A modifier creating a triangle with the incremental builder.
template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> {
public:
    Build_triangle() {}
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( 3, 1, 6);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        B.add_vertex( Point( 0, 0, 0));
        B.add_vertex( Point( 1, 0, 0));
        B.add_vertex( Point( 0, 1, 0));
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.end_surface();
    }
};


DGtal::int64_t
countLatticePointsInTetrahedra( const PointZ3 & a, const PointZ3 & b, const PointZ3 & c, const PointZ3 & d )
{
  DGtal::int64_t nb = 0;
  PointZ3 ab = b - a;
  PointZ3 bc = c - b;
  PointZ3 cd = d - c;  
  PointZ3 da = a - d;
  PointZ3 abc = ab^bc;
  PointZ3 bcd = bc^cd;
  PointZ3 cda = cd^da;
  PointZ3 dab = da^ab;
  PointZ3::Component abc_shift = abc.dot( a );
  PointZ3::Component bcd_shift = bcd.dot( b );
  PointZ3::Component cda_shift = cda.dot( c );
  PointZ3::Component dab_shift = dab.dot( d );
  if ( abc.dot( d ) < abc_shift )  { abc = -abc; abc_shift = -abc_shift; }
  if ( bcd.dot( a ) < bcd_shift )  { bcd = -bcd; bcd_shift = -bcd_shift; }
  if ( cda.dot( b ) < cda_shift )  { cda = -cda; cda_shift = -cda_shift; }
  if ( dab.dot( c ) < dab_shift )  { dab = -dab; dab_shift = -dab_shift; }
  PointZ3 inf = a.inf( b ).inf( c ).inf( d );
  PointZ3 sup = a.sup( b ).sup( c ).sup( d );
  Domain domain( inf, sup );
  for ( DomainConstIterator it = domain.begin(), itE = domain.end();
	it != itE; ++it )
    {
      PointZ3 p = *it;
      if ( ( abc.dot( p ) >= abc_shift )
	   && ( bcd.dot( p ) >= bcd_shift )
	   && ( cda.dot( p ) >= cda_shift )
	   && ( dab.dot( p ) >= dab_shift ) )
	++nb;
    }
  return nb;
}

/**
   @return 4 iff there are 4 lattice points in tetrahedra, otherwise returns 5.
*/
DGtal::int64_t
countLatticePointsInTetrahedraIf4( const PointZ3 & a, const PointZ3 & b, const PointZ3 & c, const PointZ3 & d )
{
  DGtal::int64_t nb = 0;
  PointZ3 ab = b - a;
  PointZ3 bc = c - b;
  PointZ3 cd = d - c;  
  PointZ3 da = a - d;
  PointZ3 abc = ab^bc;
  PointZ3 bcd = bc^cd;
  PointZ3 cda = cd^da;
  PointZ3 dab = da^ab;
  PointZ3::Component abc_shift = abc.dot( a );
  PointZ3::Component bcd_shift = bcd.dot( b );
  PointZ3::Component cda_shift = cda.dot( c );
  PointZ3::Component dab_shift = dab.dot( d );
  if ( abc.dot( d ) < abc_shift )  { abc = -abc; abc_shift = -abc_shift; }
  if ( bcd.dot( a ) < bcd_shift )  { bcd = -bcd; bcd_shift = -bcd_shift; }
  if ( cda.dot( b ) < cda_shift )  { cda = -cda; cda_shift = -cda_shift; }
  if ( dab.dot( c ) < dab_shift )  { dab = -dab; dab_shift = -dab_shift; }
  PointZ3 inf = a.inf( b ).inf( c ).inf( d );
  PointZ3 sup = a.sup( b ).sup( c ).sup( d );
  Domain domain( inf, sup );
  for ( DomainConstIterator it = domain.begin(), itE = domain.end();
	it != itE; ++it )
    {
      PointZ3 p = *it;
      if ( ( abc.dot( p ) >= abc_shift )
	   && ( bcd.dot( p ) >= bcd_shift )
	   && ( cda.dot( p ) >= cda_shift )
	   && ( dab.dot( p ) >= dab_shift ) )
	{ 
          ++nb;
          if ( nb > 4 ) return 5;
        }
    }
  return nb;
}


/**
   @return the number of lattice points in tetrahedra (or 5 if more
   than four). The user may specify if integer points on planes are
   counted or not (useful for checking if an open tetrahedra contains
   integer points).
  */
DGtal::int64_t
countLatticePointsInTetrahedra( const PointZ3 & a, const PointZ3 & b, 
				const PointZ3 & c, const PointZ3 & d, 
				bool abc_closed, bool bcd_closed,
				bool cda_closed, bool dab_closed )
{
  DGtal::int64_t nb = 0;
  PointZ3 ab = b - a;
  PointZ3 bc = c - b;
  PointZ3 cd = d - c;  
  PointZ3 da = a - d;
  PointZ3 abc = ab^bc;
  PointZ3 bcd = bc^cd;
  PointZ3 cda = cd^da;
  PointZ3 dab = da^ab;
  PointZ3::Component abc_shift = abc.dot( a );
  PointZ3::Component bcd_shift = bcd.dot( b );
  PointZ3::Component cda_shift = cda.dot( c );
  PointZ3::Component dab_shift = dab.dot( d );
  if ( abc.dot( d ) < abc_shift )  { abc = -abc; abc_shift = -abc_shift; }
  if ( bcd.dot( a ) < bcd_shift )  { bcd = -bcd; bcd_shift = -bcd_shift; }
  if ( cda.dot( b ) < cda_shift )  { cda = -cda; cda_shift = -cda_shift; }
  if ( dab.dot( c ) < dab_shift )  { dab = -dab; dab_shift = -dab_shift; }
  // if ( ! abc_closed ) ++abc_shift;
  // if ( ! bcd_closed ) ++bcd_shift;
  // if ( ! cda_closed ) ++cda_shift;
  // if ( ! dab_closed ) ++dab_shift;
  PointZ3 inf = a.inf( b ).inf( c ).inf( d );
  PointZ3 sup = a.sup( b ).sup( c ).sup( d );
  Domain domain( inf, sup );
  for ( DomainConstIterator it = domain.begin(), itE = domain.end();
	it != itE; ++it )
    {
      PointZ3 p = *it;
      if ( ( abc.dot( p ) >= abc_shift )
	   && ( bcd.dot( p ) >= bcd_shift )
	   && ( cda.dot( p ) >= cda_shift )
	   && ( dab.dot( p ) >= dab_shift ) )
	{ 
          ++nb;
          if ( nb > 4 ) return 5;
        }
    }
  return nb;
}


namespace DGtal {
  template <typename Kernel, typename HDS>
  class DigitalCore : public CGAL::Modifier_base<HDS> {

  public:
    typedef typename CGAL::Delaunay_triangulation_3<Kernel> Delaunay;
    typedef typename Delaunay::Cell_iterator  CellIterator;
    typedef typename Delaunay::Cell_handle    CellHandle;
    typedef typename Delaunay::Point          Point;
    typedef typename Delaunay::Facet          Facet;
    typedef typename Delaunay::Edge           Edge;
    typedef typename Delaunay::Facet_iterator FacetIterator;
    typedef typename Delaunay::Finite_facets_iterator FiniteFacetsIterator;
    typedef typename Delaunay::Finite_vertices_iterator FiniteVerticesIterator;
    typedef typename Delaunay::Vertex_handle  VertexHandle;
    typedef DGtal::uint64_t                   Size;
    typedef typename std::map< CellHandle, Size >  MapCell2Size;
    typedef typename std::set< Facet >             FacetSet;
    typedef typename std::set< CellHandle >        CellSet;
    typedef typename Kernel::Vector_3         Vector;
    typedef typename Kernel::Plane_3          Plane;

  private:
    /// The Delaunay triangulation that carries the digital core.
    Delaunay* myDelaunay;
    /// Stores the number of lattice of points within each triangulation Cell.
    MapCell2Size myNbLatticePoints;
    /// Current boundary as a set of facets.
    FacetSet myBoundary;
    /// Current interior cells.
    CellSet myInterior;
    /// Helper classes for some operations.
    Triangulation3DHelper<Delaunay,Kernel> myTH;
    /// For converting CGAL points to DGtal points.
    toDGtalFunctor<Kernel> toDGtal;

    int myVertexIdx;
    std::map<VertexHandle,int> myIndices;

    //--------------------------------------------------------------------------
  public:

    inline 
    DigitalCore( Alias<Delaunay> t ) 
      : myDelaunay( &t ), myTH( *myDelaunay )
    {
      countLatticePoints();
      computeBasicFacets();
    }

    //--------------------------------------------------------------------------
  public:


    Vector normal( const Facet & f ) const
    {
      Plane plane = myDelaunay->triangle( f ).supporting_plane();
      if ( plane.has_on_positive_side( f.first->vertex( f.second )->point() ) )
        return plane.orthogonal_vector();
      else
        {
          ASSERT( plane.has_on_negative_side( f.first->vertex( f.second )->point() ) );
          return -plane.orthogonal_vector();
        }
    }

    void getCCWVertices( VertexHandle & v0, VertexHandle & v1, VertexHandle & v2, 
                         const Facet & f ) const
    {
      v0 = f.first->vertex( (f.second+1)%4 );
      v1 = f.first->vertex( (f.second+2)%4 );
      v2 = f.first->vertex( (f.second+3)%4 );
      if ( myDelaunay->is_infinite( f.first->vertex( f.second ) ) )
        {
          Facet f2 = myDelaunay->mirror_facet( f );
          getCCWVertices( v0, v1, v2, f2 );
          std::swap( v1, v2 );
          // std::cout << "(Infinite)" << std::flush;
        }
      else
        {
          Plane plane( v0->point(), v1->point(), v2->point() );
          if ( ! plane.has_on_positive_side( f.first->vertex( f.second )->point() ) )
            {
              ASSERT( plane.has_on_negative_side( f.first->vertex( f.second )->point() ) );
              // std::cout << "(Inverted)" << std::flush;
              std::swap( v1, v2 );
            }
          // else 
          //   std::cout << "(Positive)" << std::flush;
        }
    }

    
    inline
    const FacetSet & boundary() const
    { return myBoundary; }

    inline
    FacetSet & boundary()
    { return myBoundary; }

    inline
    const CellSet & interior() const
    { return myInterior; }

    inline
    bool isFacetBasic( const Facet & f ) const
    {
      if ( ! myDelaunay->is_infinite( f ) )
        {
          int i = f.second;
          PointZ3 a( toDGtal( f.first->vertex( (i+1)%4 )->point() ) );
          PointZ3 b( toDGtal( f.first->vertex( (i+2)%4 )->point() ) );
          PointZ3 c( toDGtal( f.first->vertex( (i+3)%4 )->point() ) );
          return ( ( (b-a).norm( PointZ3::L_infty ) == 1 )
                   && ( (c-b).norm( PointZ3::L_infty ) == 1 )
                   && ( (a-c).norm( PointZ3::L_infty ) == 1 ) );
        }
      else
        return false;
    }

    inline 
    Size nbLatticePoints( const CellHandle & cell ) const
    {
      typename MapCell2Size::const_iterator it = myNbLatticePoints.find( cell );
      ASSERT( it != myNbLatticePoints.end() );
      return it->second;
    }

    inline 
    bool isCellInterior( const CellHandle & c ) const
    {
      return myInterior.find( c ) != myInterior.end();
    }

    inline 
    bool isFacetBoundary( const Facet & f ) const
    {
      return myBoundary.find( f ) != myBoundary.end();
    }

    inline 
    bool isFacetInterior( const Facet & f ) const
    {
      return myInterior.find( f.first ) != myInterior.end();
    }

    inline 
    bool isFacetExterior( const Facet & f ) const
    {
      return ( myInterior.find( f.first ) == myInterior.end() )
        && ! isFacetBoundary( f );
    }

    inline
    bool isVertexInterior( const VertexHandle & vh ) const
    {
      std::vector<CellHandle> incident_cells;
      myDelaunay->incident_cells( vh, std::back_inserter( incident_cells ) );
      for ( typename std::vector<CellHandle>::const_iterator it = incident_cells.begin(),
              itend = incident_cells.end(); it != itend; ++it )
        {
          if ( ! isCellInterior( *it ) ) return false;
        }
      return true;
    }

    inline
    double area( const Facet & f ) const
    {
      return sqrt( myDelaunay->triangle( f ).squared_area() );
    }

    void extend()
    {
      ASSERT( myDelaunay != 0 );
      trace.beginBlock("[DigitalCore] Extend boundary.");

      // Queue for computing the digital core.
      CellSet priorityQ;
      
      // Prepare queue.
      for ( typename FacetSet::const_iterator it = myBoundary.begin(), 
	      itend = myBoundary.end();
            it != itend; ++it )
        priorityQ.insert( it->first );
      
      // Start extension
      bool isMarked[ 4 ];
      unsigned int nb = 0;
      while ( ! priorityQ.empty() )
        {
          typename CellSet::iterator it = priorityQ.begin();
          CellHandle cell_h = *it;
          priorityQ.erase( it );
          // infinite cells are not processed.
          if ( myDelaunay->is_infinite( cell_h ) ) continue;
          // cells containing integer points in the interior are not processed.
          if ( nbLatticePoints( cell_h ) != 4 )         continue;
          // checking cell faces for extension.
          if ( checkCellExtension( isMarked, cell_h ) )
            {
              extendCell( priorityQ, isMarked, cell_h );
              ++nb;
            }
        }
      trace.info() << "- cells flipped = " << nb << std::endl;
      trace.info() << "- boundary has " << myBoundary.size() << " facets." << std::endl;
      trace.endBlock();
    }


    /**
       Given an edge \a e of a boundary facet \a f, returns the
       boundary facet that is incident to the same edge (according to
       the exterior adjacency).
    */
    Facet adjacentBoundaryFacet( const Edge & e, const Facet & f ) const
    {
      ASSERT( isFacetBoundary( f ) );
      Edge pivot = e;
      Facet facet = f;
      Facet result;
      do {
	myTH.nextAroundEdge( pivot, facet );
	if ( facet == f )
	  {
	    DGtal::trace.warning() << "[DigitalCore::adjacentBoundaryFacet] Has looped around the edge." << std::endl;
	    result = myDelaunay->mirror_facet( facet );
	    ASSERT( isFacetBoundary( result ) );
	    return result;
	  }
	result = myDelaunay->mirror_facet( facet );
	ASSERT( ( ! isFacetBoundary( facet ) )
		|| isFacetBoundary( result ) );
      } while ( ! isFacetBoundary( result ) );
      return result;
    }

     /**
       Given a boundary facet \a f, returns the boundary facets that
       are incident to it (according to the exterior adjacency).
    */
    void getAdjacentBoundaryFacets( std::vector<Facet> & adj_facets, const Facet & f ) const
    {
      ASSERT( isFacetBoundary( f ) );
      for ( int k = 0; k < 3; ++k )
	{
	  int i = (f.second+1+k)%4;
	  int j = (i+1)%4; 
	  j = ( j == f.second ) ? (j+1)%4 : j;
	  adj_facets.push_back( adjacentBoundaryFacet( Edge( f.first, i, j ), f ) );
	}
    }

    /**
       Removes biface facets that are dangling.
     */
    bool pruneBoundary()
    {
      FacetSet to_remove;
      std::vector<Facet> adj_facets;
      bool changes = false;
      for ( typename FacetSet::const_iterator it = boundary().begin(), 
	      itend = boundary().end();
            it != itend; ++it )
	{
	  Facet f = *it;
	  getAdjacentBoundaryFacets(adj_facets, f );
	  for ( unsigned int i = 0; i < adj_facets.size(); ++i )
	    if ( myDelaunay->mirror_facet( f ) == adj_facets[ i ] )
	      {
		to_remove.insert( f );
		to_remove.insert( adj_facets[ i ] );
	      }
	  adj_facets.clear();
	}
      DGtal::trace.info() << "Dangling facets: " << to_remove.size() << std::endl;
      for ( typename FacetSet::const_iterator it = to_remove.begin(), 
	      itend = to_remove.end();
            it != itend; ++it )
	boundary().erase( *it );
      return ! to_remove.empty();
    }

    void updateFacet( CGAL::Polyhedron_incremental_builder_3<HDS> & B,
		      const Facet & f )
    {
      int idx[ 3 ];
      VertexHandle vertices[ 3 ];
      getCCWVertices( vertices[ 0 ], vertices[ 1 ], vertices[ 2 ], f );
      for ( int j = 0; j < 3; ++j )
	{
	  VertexHandle v = vertices[ j ];
	  typename std::map<VertexHandle,int>::const_iterator 
	    it = myIndices.find( v );
	  if ( it == myIndices.end() )
	    {
	      idx[ j ] = myIndices[ v ] = myVertexIdx++;
	      B.add_vertex( v->point() );
	    } 
	  else 
	    idx[ j ] = myIndices[ v ];
	}
      B.begin_facet();
      // std::cout << "F(" << idx[0] << "," << idx[1] << "," << idx[2] << ")" << std::flush;
      B.add_vertex_to_facet( idx[ 0 ] );
      B.add_vertex_to_facet( idx[ 1 ] );
      B.add_vertex_to_facet( idx[ 2 ] );
      B.end_facet();
    }

    void operator()( HDS& hds) {
      // Postcondition: hds is a valid polyhedral surface.
      CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
      B.begin_surface( myDelaunay->number_of_vertices(), boundary().size() );

      FacetSet bdry_save = boundary();
      FacetSet visitedFacets;
      std::vector<Facet> componentFacets;
      myVertexIdx = 0;
      // Find a start facet in a connected component of the boundary.
      for ( typename FacetSet::const_iterator it = bdry_save.begin(), 
	      itend = bdry_save.end();
            it != itend; ++it )
        {
	  if ( visitedFacets.find( *it ) != visitedFacets.end() ) continue;
	  // new component;
	  std::cout << "New component..." << std::endl;
	  myIndices.clear();
	  componentFacets.clear();
	  std::vector<Facet> adj_facets;
	  std::queue<Facet> Q;
	  Q.push( *it );
	  while ( ! Q.empty() )
	    {
	      Facet f = Q.front(); Q.pop();
	      if ( visitedFacets.find( f ) != visitedFacets.end() ) continue;
	      visitedFacets.insert( f );
	      componentFacets.push_back( f );
	      updateFacet( B, f );
	      getAdjacentBoundaryFacets( adj_facets, f );
	      for ( unsigned int i = 0; i < adj_facets.size(); ++i )
		Q.push( adj_facets[ i ] );
	      adj_facets.clear();
	    }
	  for ( unsigned int i = 0; i < componentFacets.size(); ++i )
	    boundary().erase( componentFacets[ i ] );
	}
      B.end_surface();
      boundary() = bdry_save;
    }

        
    //--------------------------------------------------------------------------
  protected:

    bool checkCellExtension( bool isBoundary[ 4 ], const CellHandle & cell ) const
    {
      unsigned int n = 0;
      for ( int i = 0; i < 4; ++i )
        {
          isBoundary[ i ] = isFacetBoundary( Facet( cell, i ) );
          if ( isBoundary[ i ] ) ++n;
        }
      // At least 2 are marked, by convexity, we can close the gap
      // to the further faces of the cell.
      bool propagate = n >= 2;
      return propagate;
    }
    
    void extendCell( CellSet & priorityQ, bool isBoundary[ 4 ], const CellHandle & cell )
    { // Switch cell to interior
      myInterior.insert( cell );
      for ( unsigned int i = 0; i < 4; ++i ) {
        if ( ! isBoundary[ i ] )
          {
            Facet nfacet = myDelaunay->mirror_facet( Facet( cell, i ) );
	    if ( myInterior.find( nfacet.first ) == myInterior.end() )
	      {
		priorityQ.insert( nfacet.first );
		myBoundary.insert( nfacet );
	      }
          }
        else
          {
            myBoundary.erase( Facet( cell, i ) );
          }
      }
    }


    void computeBasicFacets()
    {
      ASSERT( myDelaunay != 0 );
      trace.beginBlock("[DigitalCore] Compute basic facets.");
      DGtal::uint64_t nb = 0;
      for( FacetIterator it = myDelaunay->facets_begin(), itend = myDelaunay->facets_end();
           it != itend; ++it )
        {
          if ( isFacetBasic( *it ) ) // || checkPlanarVFacet( f ) )
            {
	      myBoundary.insert( *it );
	      myBoundary.insert( myDelaunay->mirror_facet( *it ) );
              ++nb;
            }
        }
      trace.info() << "- nb basic facets = " << nb << std::endl;
      trace.endBlock();
    }

    void countLatticePoints()
    {
      ASSERT( myDelaunay != 0 );
      trace.beginBlock("[DigitalCore] Counting lattice points.");
      double maxbar = myDelaunay->number_of_cells();
      double bar = 0.0;
      int i = 0;
      for( CellIterator it = myDelaunay->cells_begin(), 
             itend = myDelaunay->cells_end();
           it != itend; ++it, ++bar, ++i )
        {
          if ( i % 1000 == 0 ) trace.progressBar( bar, maxbar );
          if ( myDelaunay->is_infinite( it ) ) 
            myNbLatticePoints[ it ] = -1;
          else
            {
              PointZ3 a( toDGtal(it->vertex(0)->point())),
                b(toDGtal(it->vertex(1)->point())),
                c(toDGtal(it->vertex(2)->point())),
                d(toDGtal(it->vertex(3)->point()));
              myNbLatticePoints[ it ] = countLatticePointsInTetrahedraIf4( a, b, c, d );
            }
        }
      trace.endBlock();
    }
    
  };

}

template <typename TKernel>
class PolyhedronFiller
{
  typedef TKernel                              Kernel;
  typedef typename CGAL::Polyhedron_3<Kernel>  Polyhedron;
  typedef toCGALFunctor<Kernel>       ToCGAL;
  typedef toDGtalFunctor<Kernel>      ToDGtal;
  typedef typename Polyhedron::Vertex          Vertex;
  typedef typename Polyhedron::Halfedge        HalfEdge;
  typedef typename Polyhedron::Facet           Facet;
  typedef typename Polyhedron::Vertex_handle   VertexHandle;
  typedef typename Polyhedron::Halfedge_handle HalfEdgeHandle;
  typedef typename Polyhedron::Facet_handle    FacetHandle;
  typedef typename Polyhedron::Edge_iterator   EdgeIterator;
  typedef typename Polyhedron::Facet_iterator  FacetIterator;
  typedef typename Polyhedron::Facet_const_iterator FacetConstIterator;
  typedef typename Polyhedron::Halfedge_around_facet_const_circulator 
  HalfedgeFacetConstCirculator;
  typedef typename Kernel::Point_3             Point;
  typedef typename Kernel::Vector_3            Vector;
  typedef typename Kernel::Plane_3             Plane;
  typedef typename Kernel::FT                  Component;


private:
  Polyhedron myP;
  ToDGtal toDGtal;

public:
  inline PolyhedronFiller() {}
  inline Polyhedron & P() { return myP; }
  inline const Polyhedron & P() const { return myP; }

  inline Vector cross( const Vector & a, const Vector & b ) const
  {
    return Vector( a.y() * b.z() - a.z() * b.y(),
		   a.z() * b.x() - a.x() * b.z(),
		   a.x() * b.y() - a.y() * b.x() );
  }
  
  /**
     @return the angle at half-edge \a he (positive is convex, zero is
     flat, negative is concave).
  */
  double signAngle( PointZ3 & zs, PointZ3 & zt, PointZ3 & zu, PointZ3 & zv,
		    HalfEdgeHandle he ) const
  {
    VertexHandle t = he->vertex();
    VertexHandle s = he->opposite()->vertex();
    VertexHandle u = he->next()->vertex();
    VertexHandle v = he->opposite()->next()->vertex();
    Point ps = s->point(); zs = toDGtal( ps );
    Point pt = t->point(); zt = toDGtal( pt );
    Point pu = u->point(); zu = toDGtal( pu );
    Point pv = v->point(); zv = toDGtal( pv );
    // Compute dihedral angle at (s,t)
    Plane stu( ps, pt, pu );
    Plane tsv( pt, ps, pv );
    Vector n1 = stu.orthogonal_vector(); //n1 = n1 / n1.squared_length();
    Vector n2 = tsv.orthogonal_vector(); //n2 = n2 / n2.squared_length();
    Vector c = cross( n1, n2 );
    return c * (pt -ps); // positive is convex
  }

  double signAngle( HalfEdgeHandle he ) const
  {
    PointZ3 zs, zt, zu, zv;
    return signAngle( zs, zt, zu, zv, he );
  }

  double angleFlippedEdge( HalfEdgeHandle he ) const
  {
    VertexHandle t = he->vertex();
    VertexHandle s = he->opposite()->vertex();
    VertexHandle u = he->next()->vertex();
    VertexHandle v = he->opposite()->next()->vertex();
    Point ps = s->point(); 
    Point pt = t->point(); 
    Point pu = u->point(); 
    Point pv = v->point(); 
    // Compute dihedral angle at (s,t)
    Plane svu( ps, pv, pu );
    Plane tuv( pt, pu, pv );
    double sa = sinAngle( pu - pv, svu, tuv );
    return asin( sa );
  }


  /**
     e is an oriented edge (a vector) separating two faces of
     supporting plane p1 (left of vector) and p2 (right of vector).
  */
  double sinAngle( const Vector & e, const Plane & p1, const Plane & p2 ) const
  {
    Vector n1 = p1.orthogonal_vector(); 
    Vector n2 = p2.orthogonal_vector(); 
    Vector c = cross( n1, n2 );
    return c * e 
      / sqrt( n1.squared_length() * n2.squared_length() * e.squared_length() ); 
    // sin(angle) positive is convex
  }

  bool checkFlipGeometry( HalfEdgeHandle he ) const
  {
    VertexHandle t = he->vertex();
    VertexHandle s = he->opposite()->vertex();
    VertexHandle u = he->next()->vertex();
    VertexHandle v = he->opposite()->next()->vertex();
    Point ps = s->point(); 
    Point pt = t->point(); 
    Point pu = u->point(); 
    Point pv = v->point(); 
    // Check angle at (u,s) 
    Plane stu( ps, pt, pu );
    Plane osu( he->prev()->opposite()->next()->vertex()->point(), ps, pu );
    double us_sa1 = sinAngle( ps - pu, stu, osu );
    Plane vsu( pv, ps, pu );
    double us_sa2 = sinAngle( ps - pu, stu, vsu );
    if ( us_sa2 >= us_sa1 ) 
      return false; // facet (o,s,u) would be inside tetrahedron (s,t,u,v)
    // Check angle at (s,v)
    Plane tsv( pt, ps, pv );
    Plane ovs( he->opposite()->next()->opposite()->next()->vertex()->point(), pv, ps );
    double sv_sa1 = sinAngle( pv - ps, tsv, ovs );
    double sv_sa2 = sinAngle( pv - ps, tsv, vsu );
    if ( sv_sa2 >= sv_sa1 ) 
      return false; // facet (o',v,s) would be inside tetrahedron (s,t,u,v)
    return true;
  }


  typedef typename std::set<HalfEdgeHandle> MarkHE;

  bool expand()
  {
    unsigned int nb_ccv = 0;
    unsigned int nb_edges = 0;
    PointZ3 zs, zt, zu, zv;
    MarkHE toFlip;
    MarkHE toErase;
    for ( EdgeIterator it = myP.edges_begin(), ite = myP.edges_end();
	  it != ite; ++it, ++nb_edges )
      {
	HalfEdgeHandle he = it;
        double sa = signAngle( zs, zt, zu, zv, he );
	if ( sa >= 0.0 ) continue; // positive is convex
	// Checks it it a simple hole
	if ( he->next()->opposite()->next() == he->opposite()->prev()->opposite() )
	  { // hole at t
	    int64_t nb = countLatticePointsInTetrahedra
	      ( zs, zt, zu, zv, true, true, false, true );
	    // bool stu_closed, bool tuv_closed, bool uvs_closed, bool vst_closed
	    if ( nb <= 4 ) toErase.insert( he );
	  }
	else if ( he->prev()->opposite()->prev() == he->opposite()->next()->opposite() )
	  { // hole at s
	    int64_t nb = countLatticePointsInTetrahedra
	      ( zs, zt, zu, zv, true, false, true, true );
	    // bool stu_closed, bool tuv_closed, bool uvs_closed, bool vst_closed
	    if ( nb <= 4 ) toErase.insert( he->opposite() );
	  }
	// normal concavities.
	int64_t nb = countLatticePointsInTetrahedra
	  ( zs, zt, zu, zv, true, false, false, true );
	if ( nb <= 4 )
	  {
	    // if ( ( signAngle( he->next() ) >= 0.0 )
	    // 	 && ( signAngle( he->prev() ) >= 0.0 )
	    // 	 && ( signAngle( he->opposite()->next() ) >= 0.0 )
	    // 	 && ( signAngle( he->opposite()->prev() ) >= 0.0 ) )
	    if ( checkFlipGeometry( he ) && checkFlipGeometry( he->opposite() ) )
              {
                double a = angleFlippedEdge( he );
                if ( ( a <= 0.0 ) || ( a >= M_PI ) )
                  DGtal::trace.warning() << "Invalid flipped edge angle: " << a << std::endl;
                else // if ( a > 0.0 && a < (M_PI / 2.0 ) )
                  toFlip.insert( he );
              }
	  }
	++nb_ccv;
      }
    DGtal::trace.info() << "Found concavities: " << nb_ccv  << std::endl
			<< "         to erase: " << toErase.size() << std::endl
			<< "          to flip: " << toFlip.size() << std::endl
			<< "                 / " << nb_edges << std::endl;
    bool changes = false;
    // Process holes.
    DGtal::trace.info() << "Processing holes..." << std::endl;
    while ( ! toErase.empty() )
      {
	HalfEdgeHandle st = *( toErase.begin() );
	toErase.erase( toErase.begin() );
	removeFacet( toErase, st );
	removeFacet( toErase, st->opposite() );
	removeFacet( toErase, st->next()->opposite() );
	// removeFacet( toErase, st->opposite()->next()->opposite() );
	removeFacet( toFlip, st );
	removeFacet( toFlip, st->opposite() );
	removeFacet( toFlip, st->next()->opposite() );
	// removeFacet( toFlip, st->opposite()->next()->opposite() );
	HalfEdgeHandle ut = st->next()->opposite();
	HalfEdgeHandle vt = st->opposite()->next()->next();
	bool tetra_st = st->opposite()->next()->opposite()->next() 
	  == st->next()->next()->opposite();
	bool tetra_ut = ut->opposite()->next()->opposite()->next() 
	  == ut->next()->next()->opposite();
	bool tetra_vt = vt->opposite()->next()->opposite()->next() 
	  == vt->next()->next()->opposite();
	if ( P().is_tetrahedron( st ) )
	  {
	    std::cout << "tetra to erase (stuv): " << st->vertex()->point() << std::endl;
	    P().erase_connected_component( st );
	  }
	else if ( tetra_st )
	  {
	    ASSERT( (! tetra_ut) && (! tetra_vt) );
	    std::cout << "tetra to erase (st): " << st->vertex()->point() << std::endl;
	    HalfEdgeHandle uv = st->next()->opposite()->prev();
	    removeFacet( toErase, uv );
	    removeFacet( toFlip, uv );
	    HalfEdgeHandle nuv = P().join_facet( uv );
	    st = P().erase_center_vertex( st );
	    P().erase_center_vertex( st );
	  }
	else if ( tetra_ut )
	  {
	    ASSERT( (! tetra_vt) && (! tetra_st) );
	    std::cout << "tetra to erase (ut): " << ut->vertex()->point() << std::endl;
	    HalfEdgeHandle vs = ut->next()->opposite()->prev();
	    removeFacet( toErase, vs );
	    removeFacet( toFlip, vs );
	    HalfEdgeHandle nvs = P().join_facet( vs );
	    ut = P().erase_center_vertex( ut );
	    P().erase_center_vertex( ut );
	  }
	else if ( tetra_vt )
	  {
	    ASSERT( (! tetra_st) && (! tetra_ut) );
	    std::cout << "tetra to erase (vt): " << vt->vertex()->point() << std::endl;
	    HalfEdgeHandle su = vt->next()->opposite()->prev();
	    removeFacet( toErase, su );
	    removeFacet( toFlip, su );
	    HalfEdgeHandle nsu = P().join_facet( su );
	    vt = P().erase_center_vertex( vt );
	    P().erase_center_vertex( vt );
	  }
	else
	  {
	    std::cout << "Erasing vertex: " << st->vertex()->point() << std::endl;
	    P().erase_center_vertex( st );
	  }
	changes = true;
        toErase.clear();
        toFlip.clear();
      }
    DGtal::trace.info() << "Processing flips..." << std::endl;
    while ( ! toFlip.empty() )
      {
        typename MarkHE::iterator itbest = toFlip.begin();
        double a = angleFlippedEdge( *itbest );
        for ( typename MarkHE::iterator it = itbest, ite = toFlip.end();
              it != ite; ++it )
          {
            double b = angleFlippedEdge( *it );
            if ( b < a ) 
              { a = b; itbest = it; }
          }
	HalfEdgeHandle he = *itbest;
	toFlip.erase( itbest );
	// HalfEdgeHandle he = *( toFlip.begin() );
	// toFlip.erase( toFlip.begin() );
    	if ( checkFlipGeometry( he ) && checkFlipGeometry( he->opposite() ) )
    	  {
    	    removeFacet( toFlip, he );
    	    removeFacet( toFlip, he->opposite() );
    	    std::cout << "Flipping edge: " << he->vertex()->point() 
    		      << " -> " << he->opposite()->vertex()->point() 
                      << " a=" << a << std::endl;
    	    P().flip_edge( he );
    	    changes = true;
            toFlip.clear();
    	  }
    	else
    	  std::cout << "Edge has changed: " << he->vertex()->point() 
    		    << " -> " << he->opposite()->vertex()->point() << std::endl;
      }
    // Process concavities.
    return changes;
  }

  void markFacet( MarkHE & markedHE, HalfEdgeHandle he )
  {
    HalfEdgeHandle e = he;
    do {
      markedHE.insert( e );
      e = e->next();
    } while ( e != he );
  }
  void removeFacet( MarkHE & markedHE, HalfEdgeHandle he )
  {
    HalfEdgeHandle e = he;
    do {
      typename MarkHE::iterator it = markedHE.find( e );
      if ( it != markedHE.end() ) markedHE.erase( it );
      e = e->next();
    } while ( e != he );
  }

  template <typename Viewer>
  void view( Viewer & viewer, double retract ) const
  {
    DGtal::Color colBasicFacet1( 0, 255, 0, 255 );
    for ( FacetConstIterator it = P().facets_begin(), ite = P().facets_end();
	  it != ite; ++it )
      {
	Facet f = *it;
	// HalfedgeFacetConstCirculator ci = it->facet_begin();
	HalfEdgeHandle he = f.halfedge(); // *ci; //->halfedge();
	Point a( he->vertex()->point() );
	Point b( he->next()->vertex()->point() );
	Point c( he->next()->next()->vertex()->point() );
	displayTriangle( viewer, toDGtal, a, b, c,
			 colBasicFacet1, retract );
      }
  }
  
};


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

/**
   Main function.

   @param argc the number of parameters given on the line command.

   @param argv an array of C-string, such that argv[0] is the name of
   the program, argv[1] the first parameter, etc.

  "29*z+47*y+23*x-5"
  "10*z-x^2+y^2-100"
  "z*x*y+x^4-5*x^2+2*y^2*z-z^2-1000"
  "(15*z-x^2+y^2-100)*(x^2+y^2+z^2-1000)" nice
  "(x^2+y^2+(z+5)^2-100)^2+10*(x^2+y^2+(z-3)^2)-2000" bowl
  "x^2+y^2+2*z^2-x*y*z+z^3-100" dragonfly
  "(x^2+y^2+(z+5)^2)^2-x^2*(z^2-x^2-y^2)-100" joli coeur
  "0.5*(z^2-4*4)^2+(x^2-7*7)^2+(y^2-7*7)^2-7.2*7.2*7.2*7.2" convexites et concavites
  "(1-(x^2+y^2))^2*(1+(x^2+y^2))^2-z" "cup"
*/
int main( int argc, char ** argv ) {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_3                          Point3;
  // typedef CGAL::Simple_cartesian<double>     Kernel;
  typedef CGAL::Polyhedron_3<K>               Polyhedron;
  typedef CGAL::Triangulation_3<K>            Triangulation;
  typedef CGAL::Delaunay_triangulation_3<K>   Delaunay;
  typedef Polyhedron::HalfedgeDS              HalfedgeDS;
  typedef Delaunay::Finite_cells_iterator     FiniteCellsIterator;
  
  using namespace DGtal;

  typedef KhalimskySpaceND<3,DGtal::int64_t> K3;
  typedef Z3::Vector Vector;
  typedef Z3::RealPoint RealPoint;
  typedef K3::SCell SCell;
  typedef K3::SCellSet SCellSet;
  typedef SCellSet::const_iterator SCellSetConstIterator;
  // typedef ImplicitRoundedHyperCube<Z3> Shape;
  typedef RealPoint::Coordinate Ring;
  typedef toCGALFunctor<K>                    ToCGAL;
  typedef toDGtalFunctor<K>                   ToDGtal;
  typedef DigitalCore<K,HalfedgeDS>           DigCore;
  typedef DigCore::FacetSet                   FacetSet;
  typedef DigCore::Facet                      Facet;

  QApplication application(argc,argv); // remove Qt arguments.

  po::options_description general_opt("Specific allowed options (for Qt options, see Qt official site) are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("vol,v", po::value<std::string>(), "specifies the shape as some subset of a .vol file [arg]" )
    ("min,m", po::value<int>()->default_value( 1 ), "the minimum threshold in the .vol file: voxel x in shape iff m < I(x) <= M" )
    ("max,M", po::value<int>()->default_value( 255 ), "the maximum threshold in the .vol file: voxel x in shape iff m < I(x) <= M" )
    ("poly,p", po::value<std::string>(), "specifies the shape as the zero-level of the multivariate polynomial [arg]" )
    ("gridstep,g", po::value<double>()->default_value( 1.0 ), "specifies the digitization grid step ([arg] is a double) when the shape is given as a polynomial." )
    ("bounds,b", po::value<int>()->default_value( 20 ), "specifies the diagonal integral bounds [-arg,-arg,-arg]x[arg,arg,arg] for the digitization of the polynomial surface." )
    ("point-list,l", po::value<std::string>(), "specifies the shape as a list of digital points given in the filename [arg], x y z per line." )
    ("retract,r", po::value<double>()->default_value( 0.0 ), "retracts facets with factor [arg] at display time." )
    ("geometry,G", po::value<int>()->default_value( 0 ), "specifies how digital points are defined from the digital surface: arg=0: inner voxels, arg=1: inner and outer voxels, arg=2: pointels" )
    ("view,V", po::value<int>()->default_value( 1 ), "specifies what is displayed. The bits of [arg] indicates what is displayed:  0: digital core, 1: empty unreachable tetrahedra, 2: retracted tetrahedra, 3: digital points" )
    ; 
  
  // parse command line ----------------------------------------------
  bool parseOK=true;
  po::variables_map vm;
  try {
    po::command_line_parser clp( argc, argv );
    clp.options( general_opt );
    po::store( clp.run(), vm );
  } catch( const std::exception& ex ) {
    parseOK = false;
    trace.info() << "Error checking program options: "<< ex.what() << std::endl;
  }
  po::notify( vm );    
  if( !parseOK || vm.count("help")||argc<=1
      || !( vm.count("poly") || vm.count("vol") || vm.count("point-list") ) )
    {
      std::cout << "Usage: " << argv[0] << " [options] {--vol <vol-file> || --poly <polynomial-string>}\n"
		<< "Computes the linear reconstruction of the given digital surface, specified either as a thresholded .vol file, or a zero-level of a multivariate polynomial, or a list of digital points.\n"
		<< general_opt << "\n\n";
      std::cout << "Example:\n"
		<< argv[0] << " -p \"x^2+y^2+2*z^2-x*y*z+z^3-100\" -g " << (double) 0.5 << std::endl;
      return 0;
    }
  
  
  // Construction of the digital surface ----------------------------------------------
  trace.beginBlock("Construction of the digital surface");
  K3 ks;
  SCellSet boundary;
  //  SurfelAdjacency<3> sAdj( true );
  int ok = 4;
  if ( vm.count( "vol" ) )
    { // .vol file
      ok = makeSpaceAndBoundaryFromVolFile( ks, boundary,
					    vm["vol"].as<std::string>(),
					    vm["min"].as<int>(),
					    vm["max"].as<int>() );
     }
  else if ( vm.count( "poly" ) )
    {
      ok = makeSpaceAndBoundaryFromPolynomialString( ks, boundary,
						     vm[ "poly" ].as<std::string>(),
						     vm[ "gridstep" ].as<double>(),
						     vm[ "bounds" ].as<int>() );
    }
  else if ( vm.count( "point-list" ) )
    {
      ok = makeSpaceAndBoundaryFromPointList( ks, boundary,
                                              vm[ "point-list" ].as<std::string>() );
    }
  trace.endBlock();
  if ( ok != 0 ) return ok; // error.


  trace.beginBlock("Constructing the set of points");
  std::vector<PointZ3> digital_points;
  int geometry = vm["geometry"].as<int>();
  if ( geometry == 0 )
    getInnerVoxelCoordinates( digital_points, ks, boundary.begin(), boundary.end() );
  else if ( geometry == 1 )
    getInnerAndOuterVoxelCoordinates( digital_points, ks, boundary.begin(), boundary.end() );
  else if ( geometry == 2 )
    getPointelCoordinates( digital_points, ks, boundary.begin(), boundary.end() );

  trace.beginBlock("Shuffle points.");
  random_shuffle( digital_points.begin(), digital_points.end() );
  trace.endBlock();
  trace.endBlock();

  trace.beginBlock("Creating the Delaunay complex.");
  ToCGAL toCGAL;
  ToDGtal toDGtal;
  Delaunay T;
  double setsize = (double) digital_points.size()-1;
  trace.info() << "Vertices to process: " << setsize << std::endl;
  double step = 0.0;
  int i = 0;
  for ( std::vector<PointZ3>::const_iterator it = digital_points.begin(), 
	  itE = digital_points.end();
	it != itE; ++it, ++step, ++i )
    {
      if ( i % 1000 == 0 ) trace.progressBar( step, setsize );
      T.insert( toCGAL( *it ) );
    }
  trace.endBlock();

  DigCore core( T );
  // core.duplicateBoundary();
  core.extend();

  // start viewer
  int view = vm["view"].as<int>();
  double retract = vm["retract"].as<double>();

  Viewer3D<> viewerCore;
  viewerCore.show();
  Color colBasicFacet2( 0, 255, 255, 255 );
  Color colBasicFacet1( 0, 255, 0, 255 );
  Color colSpuriousTetrahedra( 255, 0, 0, 100 );
  if ( view & 0x1 ) { // View digital core.
    for ( FacetSet::const_iterator it = core.boundary().begin(), itend = core.boundary().end();
          it != itend; ++it )
      {
        // we display it.
        // Triangle triangle = t.triangle( *it );
	PointZ3 a( toDGtal( it->first->vertex( (it->second+1)%4 )->point() ) );
	PointZ3 b( toDGtal( it->first->vertex( (it->second+2)%4 )->point() ) );
	PointZ3 c( toDGtal( it->first->vertex( (it->second+3)%4 )->point() ) );
        Facet f2 = T.mirror_facet( *it );
        if ( core.isFacetBoundary( f2 ) )
          { // the mirror facet is also in the triangulation. 
            // We need to move vertices a little bit when two triangles are at the same position.
            PointZ3 n = (b-a)^(c-a);
            double norm = n.norm(PointZ3::L_2);
            double dx[ 3 ];
            for ( unsigned int j = 0; j < 3; ++j )
              dx[ j ] = 0.001*((double) n[j])/norm;
            Point3 pa( (double) a[ 0 ] + dx[ 0 ], (double) a[ 1 ] +  dx[ 1 ], (double) a[ 2 ] + dx[ 2 ] );
            Point3 pb( (double) b[ 0 ] + dx[ 0 ], (double) b[ 1 ] +  dx[ 1 ], (double) b[ 2 ] + dx[ 2 ] );
            Point3 pc( (double) c[ 0 ] + dx[ 0 ], (double) c[ 1 ] +  dx[ 1 ], (double) c[ 2 ] + dx[ 2 ] );
            displayTriangle( viewerCore, toDGtal, pa, pc, pb, colBasicFacet2, retract );
          }
        else
          {
            Point3 pa( (double) a[ 0 ], (double) a[ 1 ], (double) a[ 2 ] );
            Point3 pb( (double) b[ 0 ], (double) b[ 1 ], (double) b[ 2 ] );
            Point3 pc( (double) c[ 0 ], (double) c[ 1 ], (double) c[ 2 ] );
            displayTriangle( viewerCore, toDGtal, pa, pb, pc, colBasicFacet1, retract );
          }
      }
  } //  if ( view & 0x1 ) {

  // start viewer
  // int view = vm["view"].as<int>();
  // Viewer3D viewerRCH;
  // viewerRCH.show();
  // Color colBasicFacet1( 0, 255, 0, 255 );
  // if ( view & 0x1 ) { // View digital core.
  //   for ( FiniteCellsIterator it = T.finite_cells_begin(), 
  //           itend = T.finite_cells_end(); it != itend; ++it )
  //     {
  //       // we display it.
  // 	PointZ3 a( toDGtal( it->vertex( 0 )->point() ) );
  // 	PointZ3 b( toDGtal( it->vertex( 1 )->point() ) );
  // 	PointZ3 c( toDGtal( it->vertex( 2 )->point() ) );
  // 	PointZ3 d( toDGtal( it->vertex( 3 )->point() ) );
  // 	if ( countLatticePointsInTetrahedraIf4( a, b, c, d ) == 4 )
  // 	  displayCell( viewerRCH, toDGtal, it, colBasicFacet1 );
  //     }
  // }
  viewerCore << Viewer3D<>::updateDisplay;
  application.exec();

  std::cout << "number of vertices :  " ;
  std::cout << T.number_of_vertices() << std::endl;
  std::cout << "number of edges :  " ;
  std::cout << T.number_of_edges() << std::endl;
  std::cout << "number of facets :  " ;
  std::cout << T.number_of_facets() << std::endl;
  std::cout << "number of cells :  " ;
  std::cout << T.number_of_cells() << std::endl;

  while ( core.pruneBoundary() ) 
    ;

  PolyhedronFiller<K> PF;
  Polyhedron & P = PF.P();
  P.delegate( core );
  i = 0;
  while ( PF.expand() )
    {
      if ( i++ % 1000 == 0 )
        {
          Viewer3D<> viewer3d;
          viewer3d.show();
          PF.view( viewer3d, retract );
          viewer3d << Viewer3D<>::updateDisplay;
          application.exec();
          std::ostringstream sname;
          sname << "tmp/bdry-" << std::setfill( '0' ) << std::setw( 2 ) << (i/1000) << ".off";
          std::ofstream file_off( sname.str().c_str() );
          file_off << P;
          file_off.close();
        }
    }
  std::cout << "Polyhedron number of vertices :  " ;
  std::cout << P.size_of_vertices() << std::endl;

  Point3 low = toCGAL( ks.lowerBound() );
  Point3 high = toCGAL( ks.upperBound() );
  std::ofstream file_off( "bdry.off" );
  // CGAL::Geomview_stream gv( CGAL::Bbox_3( low.x(), low.y(), low.z(), high.x(), high.y(), high.z() ) );
  // gv << P;
  // std::cout << "Enter a key to finish" << std::endl;
  // char ch;
  // std::cin >> ch;
  file_off << P;
  file_off.close();

  // Build_triangle<HalfedgeDS> triangle;
  // P.delegate( triangle);
  // CGAL_assertion( P.is_triangle( P.halfedges_begin()));
  return 0;
}
