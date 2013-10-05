///////////////////////////////////////////////////////////////////////////////
// New version of 3d-triangulation which works as follows
// (1) T <- Delaunay triangulation( P ), Result is R
// (2) B <- Basic facets of T
// (3) expend from B by putting in R tetrahedra touching B and containing no integer points
// (4) during expansion, associates to tetrahedra the sets of integer points (increased at each expansion)
// (5) obtain boundary polyhedral surface (with sets of points for each face.
// (6) on the surface, merge concave regions if the union of points is still a plane.


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/pending/disjoint_sets.hpp>

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
// #include "DGtal/geometry/surfaces/COBAGenericNaivePlaneComputer.h"
// #include "DGtal/geometry/surfaces/ChordGenericNaivePlaneComputer.h"
#include "DGtal/geometry/surfaces/ChordGenericStandardPlaneComputer.h"

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

typedef DGtal::SpaceND<3, DGtal::int32_t> Z3;
typedef Z3::Point PointZ3;
typedef Z3::RealPoint RealPointZ3;
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
  viewer.addTriangle( a.x(), a.y(), a.z(),
		      c.x(), c.y(), c.z(),
		      b.x(), b.y(), b.z(),
		      col );
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
  RealPointZ3 A( a.x(), a.y(), a.z() );
  RealPointZ3 B( b.x(), b.y(), b.z() );
  RealPointZ3 C( c.x(), c.y(), c.z() );
  viewer.setFillColor( col );
  viewer.addTriangle( A, C, B );
}

template <typename Viewer, typename ToDGtal, typename Cell>
void displayCell( Viewer & viewer, const ToDGtal & toDGtal,
		  const Cell & c, const DGtal::Color & col,
                  double retract )
{
  for ( int i = 0; i < 4; ++i )
    displayFacet( viewer, toDGtal, std::make_pair( c, i ), col, retract );
}

template <typename Viewer, typename ToDGtal>
void displayPoint( Viewer & viewer, const ToDGtal & toDGtal,
		   typename ToDGtal::Point3 a,
		   const DGtal::Color & col, 
		   double retract )
{
  PointZ3 p = toDGtal( a );
  viewer << DGtal::CustomColors3D( col, col ) << p;
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
    typedef unsigned int                      Index;
    typedef typename std::vector<Index>       IndexList;
    typedef typename std::map<Facet,IndexList> MapFacet2Points;
    typedef typename HDS::Face_handle         PolyhedronFaceHandle;
    typedef typename std::map<Facet,PolyhedronFaceHandle> MapFacet2PolyhedronFaceHandle;
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

    // Stores the digital points.
    std::vector<Point> myPoints;
    // Stores the map that associates the points to each facet (more
    // precisely, the index of the points in myPoints).
    MapFacet2Points myFacet2Points;

    // Stores the correspondence between boundary facets of the core
    // and facets of the constructed polyhedron (after a call to operator()).
    MapFacet2PolyhedronFaceHandle myFacet2PolyhedronFaceHandle;

    //--------------------------------------------------------------------------
  public:

    inline 
    DigitalCore( Alias<Delaunay> t ) 
      : myDelaunay( t ), myTH( *myDelaunay )
    {
      countLatticePoints();
      computeBasicFacets();
      std::map<VertexHandle,Index> v2i;
      Index idx = 0;
      Index j = 0;
      for ( typename FacetSet::const_iterator it = myBoundary.begin(), 
	      itend = myBoundary.end();
            it != itend; ++it )
	{
	  Facet f = *it;
	  IndexList indices;
	  for ( int i = 0; i < 3; ++i )
	    {
	      VertexHandle v = f.first->vertex( (f.second+1+i)%4 );
	      typename std::map<VertexHandle,Index>::iterator itv = v2i.find( v );
	      if ( itv == v2i.end() )
		{
		  j = v2i[ v ] = idx++;
		  myPoints.push_back( v->point() );
		}
	      else
		j = itv->second;
	      indices.push_back( j );
	    }
	  std::sort( indices.begin(), indices.end() );
	  myFacet2Points[ f ] = indices;
	}
    }

    //--------------------------------------------------------------------------
  public:
    
    inline
    const Point & getPoint( Index i ) const
    { return myPoints[ i ]; }

    const IndexList & getPointIndices( const Facet & f ) const
    {
      typename MapFacet2Points::const_iterator it = myFacet2Points.find( f );
      ASSERT( it != myFacet2Points.end() && "[DigitalCore::getPoints] Invalid facet." );
      return it->second;
    }

    inline
    PolyhedronFaceHandle polyhedronFacet( const Facet & f ) const
    {
      typename MapFacet2PolyhedronFaceHandle::const_iterator it 
	= myFacet2PolyhedronFaceHandle.find( f );
      ASSERT( it != myFacet2PolyhedronFaceHandle.end() 
	      && "[DigitalCore::polyhedronFacet] Invalid facet." );
      return it->second;
    }

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
      PolyhedronFaceHandle fh = B.begin_facet();
      myFacet2PolyhedronFaceHandle[ f ] = fh;
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
      //--vvv Taking care of points within facet.
      IndexList points;
      for ( unsigned int i = 0; i < 4; ++i ) {
	if ( isBoundary[ i ] ) // adds the points of those facets to the following ones
	  {
	    IndexList result;
	    IndexList & points_facet_i = myFacet2Points[ Facet( cell, i ) ];
	    std::set_union( points.begin(), points.end(),
			    points_facet_i.begin(), points_facet_i.end(),
			    std::back_inserter( result ) );
	    std::swap( points, result );
	  }
      }
      //--^^^ Taking care of points within facet.
      for ( unsigned int i = 0; i < 4; ++i ) {
        if ( ! isBoundary[ i ] )
          {
            Facet nfacet = myDelaunay->mirror_facet( Facet( cell, i ) );
	    if ( myInterior.find( nfacet.first ) == myInterior.end() )
	      {
		priorityQ.insert( nfacet.first );
		myBoundary.insert( nfacet );
		//--- Taking care of points within facet.
		myFacet2Points[ nfacet ] = points;
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
class PolyhedronHull
{
  typedef PolyhedronHull<TKernel>              Self;
  typedef TKernel                              Kernel;
  typedef typename CGAL::Polyhedron_3<Kernel>  Polyhedron;
  typedef toCGALFunctor<Kernel>       ToCGAL;
  typedef toDGtalFunctor<Kernel>      ToDGtal;
  typedef typename Polyhedron::Vertex          Vertex;
  typedef typename Polyhedron::Halfedge        HalfEdge;
  typedef typename Polyhedron::Facet           Facet;
  typedef typename Polyhedron::Vertex_handle   VertexHandle;
  typedef typename Polyhedron::Halfedge_handle HalfEdgeHandle;
  typedef typename Polyhedron::Halfedge_const_handle HalfEdgeConstHandle;
  typedef typename Polyhedron::Facet_handle    FacetHandle;
  typedef typename Polyhedron::Facet_const_handle FacetConstHandle;
  typedef typename Polyhedron::Vertex_iterator VertexIterator;
  typedef typename Polyhedron::Edge_iterator   EdgeIterator;
  typedef typename Polyhedron::Facet_iterator  FacetIterator;
  typedef typename Polyhedron::Facet_const_iterator FacetConstIterator;
  typedef typename Polyhedron::Halfedge_around_facet_const_circulator 
  HalfedgeFacetConstCirculator;
  typedef typename Kernel::Point_3             Point;
  typedef typename Kernel::Vector_3            Vector;
  typedef typename Kernel::Plane_3             Plane;
  typedef typename Kernel::FT                  Component;

  typedef typename std::set<HalfEdgeHandle>    MarkHE;

  typedef typename DGtal::SpaceND<3,int>       DSpace;
  typedef typename DSpace::Point               DPoint;
  typedef typename DSpace::RealVector          DRealVector;
  typedef typename DGtal::ChordGenericStandardPlaneComputer<DSpace,DPoint,DGtal::int64_t> DigitalPlane;
  struct FacetInformation {
    DigitalPlane plane;
  };

  typedef std::map< HalfEdgeHandle, unsigned int> RankMap;
  typedef std::map< HalfEdgeHandle, HalfEdgeHandle> ParentMap;
  
  // Adapter to map: Vertex -> unsigned int, used by disjoint_sets.
  // Adapter to map: Vertex -> Vertex, used by disjoint_sets.
  typedef boost::associative_property_map<RankMap> RankAssociativePropertyMap;
  typedef boost::associative_property_map<ParentMap> ParentAssociativePropertyMap;
  typedef boost::disjoint_sets< RankAssociativePropertyMap, ParentAssociativePropertyMap > DisjointSets;


  /// The mapping from a facet (defined by its smallest half-edge) to plane information.
  typedef typename std::map<HalfEdgeHandle, FacetInformation> HE2Info;
  typedef typename HE2Info::iterator HE2InfoIterator;

private:
  Polyhedron myP;
  ToDGtal toDGtal;
  ToCGAL toCGAL;
  HE2Info myPlaneMap;
  RankMap myRankMap;
  RankAssociativePropertyMap myRankAssociativePropertyMap;
  ParentMap myParentMap;
  ParentAssociativePropertyMap myParentAssociativePropertyMap;
  DisjointSets myDisjointSets;
  std::map<HalfEdgeHandle,DGtal::Color> myColormap;
  
public:
  inline PolyhedronHull() 
    : myRankMap(), myRankAssociativePropertyMap( myRankMap ),
      myParentMap(), myParentAssociativePropertyMap( myParentMap ),
      myDisjointSets( myRankAssociativePropertyMap, myParentAssociativePropertyMap )
  {}

  inline Polyhedron & P() { return myP; }
  inline const Polyhedron & P() const { return myP; }

  inline HalfEdgeHandle nextCWAroundVertex( HalfEdgeHandle he ) const
  { return he->opposite()->next(); }
  inline HalfEdgeConstHandle nextCWAroundVertex( HalfEdgeConstHandle he ) const
  { return he->opposite()->next(); }
  inline unsigned int degreeBase( HalfEdgeConstHandle e ) const
  { 
    HalfEdgeConstHandle s = e;
    unsigned int j = 0;
    do {
      ++j;
      s = nextCWAroundVertex( s );
    } while ( s != e );
    return j;
  }

  inline HalfEdgeConstHandle smallestHalfEdge( HalfEdgeConstHandle e ) const
  {
    HalfEdgeConstHandle start = e;
    HalfEdgeConstHandle m = e;
    do {
      e = e->next();
      if ( e < m ) m = e;
    } while( e != start );
    return m;
  }

  inline HalfEdgeConstHandle smallestHalfEdge( FacetConstHandle f ) const
  { 
    HalfEdgeConstHandle h = f->halfedge();
    return smallestHalfEdge( h ); 
  }

  inline HalfEdgeHandle smallestHalfEdge( HalfEdgeHandle e ) 
  {
    HalfEdgeHandle start = e;
    HalfEdgeHandle m = e;
    do {
      e = e->next();
      if ( e < m ) m = e;
    } while( e != start );
    return m;
  }

  inline HalfEdgeHandle smallestHalfEdge( FacetHandle f ) 
  { 
    HalfEdgeHandle h = f->halfedge();
    return smallestHalfEdge( h ); 
  }

  template <typename Core>
  void initFromCore( Core & core )
  {
    // Get surface
    myP.delegate( core ); 
    // Initiate planes. Assume same order for facets in Core and facets in P.
    unsigned int j = 0;
    for ( typename Core::FacetSet::const_iterator it = core.boundary().begin(), 
	    ite = core.boundary().end();
	  it != ite; ++it, ++j )
      {
	typename Core::Facet f = *it;
	FacetHandle fh = core.polyhedronFacet( f );
	HalfEdgeHandle facet_he = smallestHalfEdge( fh );
	myDisjointSets.make_set( facet_he );
	typename Core::IndexList l = core.getPointIndices( f );
	FacetInformation & info = myPlaneMap[ facet_he ];
	// info.plane.init( 200, 3, 1 );
	info.plane.init( 2, 1 );
	std::vector<DPoint> dpoints( l.size() );
	for ( unsigned int i = 0; i < l.size(); ++i )
	  dpoints[ i ] = toDGtal( core.getPoint( l[ i ] ) );
	// std::cout << "Facet " << j << " (" << dpoints.size() << ") " << std::flush;
	info.plane.extend( dpoints[ 0 ] );
	bool ok = info.plane.extend( dpoints.begin()+1, dpoints.end() );
	if ( ! ok ) 
	  DGtal::trace.warning() << "[PolyhedronHull::initFromCore] Invalid plane." << std::endl;
	// std::cout << info.plane << std::endl;
      }

  }

  // Still may have errors.
  bool expand()
  {
    unsigned int nb_join = 0;
    unsigned int nb_edges = 0;
    MarkHE toJoin;
    std::cout << "Searching concavities: " << std::flush;
    for ( EdgeIterator it = myP.edges_begin(), ite = myP.edges_end();
	  it != ite; ++it, ++nb_edges )
      {
	HalfEdgeHandle he = it;
	// Get two neighboring facets
	// if ( signAngle( he ) > 0.0 ) continue; // do not touch convex pieces.
	toJoin.insert( he );
      }
    std::cout << toJoin.size() << " found." << std::endl;
    while ( ! toJoin.empty() )
      {
	HalfEdgeHandle he = *( toJoin.begin() );
	toJoin.erase( toJoin.begin() );
	// Check fusion
	HalfEdgeHandle he_face1 = smallestHalfEdge( he );
	HalfEdgeHandle he_face2 = smallestHalfEdge( he->opposite() );
	ASSERT( myPlaneMap.find( he_face1 ) != myPlaneMap.end() );
	ASSERT( myPlaneMap.find( he_face2 ) != myPlaneMap.end() );
	FacetInformation & info1 = myPlaneMap.find( he_face1 )->second;
	FacetInformation & info2 = myPlaneMap.find( he_face2 )->second;
	ASSERT( ! info1.plane.empty() );
	ASSERT( ! info2.plane.empty() );
	if ( info1.plane.size() == 1 ) continue; // invalid initial plane
	if ( info2.plane.size() == 1 ) continue; // invalid initial plane
	// std::cout << "---------------------------------------------------" << std::endl;
	// std::cout << "Plane 1: " << info1.plane << std::endl;
	// std::cout << "Plane 2: " << info2.plane << std::endl;
	DigitalPlane planeTester = info1.plane;
	bool ok = planeTester.extend( info2.plane.begin(), info2.plane.end() );
	if ( ! ok ) continue;
	if ( ( CGAL::circulator_size( he->vertex_begin() ) <= 3 )
	     || ( CGAL::circulator_size( he->opposite()->vertex_begin() ) <= 3 ) )
	  {
	    // std::cout << "Join_facet would create a non-manifold." << std::endl;
	    continue;
	  }
	removeFacet( toJoin, he_face1 );
	removeFacet( toJoin, he_face2 );
        HalfEdgeHandle nhe = myP.join_facet( he );
        HalfEdgeHandle nhe_face = smallestHalfEdge( nhe );
	myPlaneMap.erase( he_face1 );
	myPlaneMap.erase( he_face2 );
	typename HE2Info::iterator itnew = myPlaneMap.find( nhe_face );
	if ( itnew == myPlaneMap.end() )
	  {
	    FacetInformation info;
	    info.plane = planeTester;
	    myPlaneMap[ nhe_face ] = info;
	  }
	else
	  {
	    ASSERT( false );
	    itnew->second.plane = planeTester;
	  }
	++nb_join;
      }
    std::cout << "**********************************************************************"
	      << std::endl;
    std::cout << "Joining facets: " << nb_join << " / " << nb_edges << std::endl;
    return nb_join != 0;
  }  

  // Does not seg. fault anymore. However some natural merges are
  // blocked by the fact that all degrees stay >= 4.
  bool expandByVertices()
  {
    unsigned int nb_join = 0;
    unsigned int nb_vertices = 0;
    MarkHE toErase;
    std::cout << "Searching concavities: " << std::flush;
    for ( VertexIterator it = myP.vertices_begin(), ite = myP.vertices_end();
	  it != ite; ++it, ++nb_vertices )
      {
	HalfEdgeHandle he = it->halfedge();
	HalfEdgeHandle he_around_vtx = he;
	bool convex = true;
	do {
	  if ( signAngle( he_around_vtx ) < 0.0 ) convex = false;
	  he_around_vtx = nextCWAroundVertex( he_around_vtx );
	} while ( ( he_around_vtx != he ) );
	// Get two neighboring facets
	// if ( ! convex )
	  toErase.insert( he );
      }
    std::cout << toErase.size() << " found." << std::endl;
    while ( ! toErase.empty() )
      {
	HalfEdgeHandle he = *( toErase.begin() );
	toErase.erase( toErase.begin() );
	if ( degreeBase( he ) <= 2 ) continue;
	// Check fusion
	std::list<HalfEdgeHandle> faces;
	std::list<DigitalPlane*>  planes;
	std::vector<unsigned int> degrees;
	unsigned int degree = 0;
	HalfEdgeHandle he_around_vtx = he;
	do {
	  HalfEdgeHandle he_face = smallestHalfEdge( he_around_vtx );
	  ASSERT( myPlaneMap.find( he_face ) != myPlaneMap.end() );
	  faces.push_back( he_face );
	  planes.push_back( &( myPlaneMap.find( he_face )->second.plane ) );
	  degrees.push_back( degreeBase( he_around_vtx->opposite() ) );
	  he_around_vtx = nextCWAroundVertex( he_around_vtx );
	  ++degree;
	} while ( he_around_vtx != he );

	typename std::list<DigitalPlane*>::const_iterator itplanes = planes.begin();
	DigitalPlane planeTester = *( *itplanes++ );
	ASSERT( ! planeTester.empty() );
	if ( planeTester.size() == 1 ) continue; // invalid initial plane
	bool ok = true;
	for ( ; ok && ( itplanes != planes.end() ); ++itplanes )
	  {
	    const DigitalPlane & otherPlane = *( *itplanes );
	    ASSERT( ! otherPlane.empty() );
	    if ( otherPlane.size() == 1 ) // invalid initial plane
	      { ok = false; continue; } 
	    ok = planeTester.extend( otherPlane.begin(), otherPlane.end() );
	  }
	
	if ( ! ok ) continue;
	std::cout << "--- Erasing vertex " << he->opposite()->vertex()->point() 
		  << " of degree " << degree << ".";
	for ( unsigned int i = 0; i < degrees.size(); ++i )
	  std::cout << " " << degrees[ i ];
	bool bad_degree = false;
	for ( unsigned int i = 0; i < degrees.size(); ++i )
	  if ( degrees[ i ] <= 3 )
	    { bad_degree = true; break; }
	if ( bad_degree )
	  {
	    std::cout << " Cancelled." << std::endl;
	    continue;
	  }
	std::cout << std::endl;
	for ( typename std::list<HalfEdgeHandle>::const_iterator itf = faces.begin(),
		itfend = faces.end(); itf != itfend; ++itf )
	  {
	    removeFacet( toErase, *itf );
	    myPlaneMap.erase( *itf );
	  }
	// Look for polyhedron case.
	if ( P().is_tetrahedron( he ) )
	  {
	    DGtal::trace.warning() << "--- Erasing tetrahedron: " 
				   << he->vertex()->point() << std::endl;
	    P().erase_connected_component( he );
	    ++nb_join;
	    continue;
	  }
	// Look for almost-polyhedron case.
	HalfEdgeHandle nhe = P().erase_center_vertex( he->opposite() );
	HalfEdgeHandle nhe_face = smallestHalfEdge( nhe );
	typename HE2Info::iterator itnew = myPlaneMap.find( nhe_face );
	if ( itnew == myPlaneMap.end() )
	  {
	    FacetInformation info;
	    info.plane = planeTester;
	    myPlaneMap[ nhe_face ] = info;
	  }
	else
	  {
	    DGtal::trace.warning() << "--- Weird: faces has already a plane." 
				   << std::endl;
	    // ASSERT( false );
	    itnew->second.plane = planeTester;
	  }
	++nb_join;
      }
    std::cout << "**********************************************************************"
	      << std::endl;
    std::cout << "Erase vertices: " << nb_join << " / " << nb_vertices << std::endl;
    return nb_join != 0;
  }  

  struct HalfEdgeHandleComparator {

    Self & PH;
    
    HalfEdgeHandleComparator( Self & polyhedron_hull )
      : PH( polyhedron_hull )
    {}
    
    bool operator()( const HalfEdgeHandle he1, const HalfEdgeHandle he2 ) const
    {
	HalfEdgeHandle rep11 = PH.myDisjointSets.find_set( PH.smallestHalfEdge( he1 ) );
	HalfEdgeHandle rep12 = PH.myDisjointSets.find_set( PH.smallestHalfEdge( he1->opposite() ) );
	if ( rep11 == rep12 ) return true;
	HalfEdgeHandle rep21 = PH.myDisjointSets.find_set( PH.smallestHalfEdge( he2 ) );
	HalfEdgeHandle rep22 = PH.myDisjointSets.find_set( PH.smallestHalfEdge( he2->opposite() ) );
	if ( rep21 == rep22 ) return false;
	// Check size
	HalfEdgeHandle he_face11 = PH.smallestHalfEdge( rep11 );
	HalfEdgeHandle he_face12 = PH.smallestHalfEdge( rep12 );
	HalfEdgeHandle he_face21 = PH.smallestHalfEdge( rep21 );
	HalfEdgeHandle he_face22 = PH.smallestHalfEdge( rep22 );
	ASSERT( PH.myPlaneMap.find( he_face11 ) != PH.myPlaneMap.end() );
	ASSERT( PH.myPlaneMap.find( he_face12 ) != PH.myPlaneMap.end() );
	ASSERT( PH.myPlaneMap.find( he_face21 ) != PH.myPlaneMap.end() );
	ASSERT( PH.myPlaneMap.find( he_face22 ) != PH.myPlaneMap.end() );
	FacetInformation & info11 = PH.myPlaneMap.find( he_face11 )->second;
	FacetInformation & info12 = PH.myPlaneMap.find( he_face12 )->second;
	FacetInformation & info21 = PH.myPlaneMap.find( he_face21 )->second;
	FacetInformation & info22 = PH.myPlaneMap.find( he_face22 )->second;
	ASSERT( ! info11.plane.empty() );
	ASSERT( ! info12.plane.empty() );
	ASSERT( ! info21.plane.empty() );
	ASSERT( ! info22.plane.empty() );
        typename DigitalPlane::Size s1 = info11.plane.size() + info12.plane.size();
        typename DigitalPlane::Size s2 = info21.plane.size() + info22.plane.size();
        bool convex1 = PH.signAngle( he_face11, he_face12 ) > 0.0;
        bool convex2 = PH.signAngle( he_face21, he_face22 ) > 0.0;
        if ( convex1 ) 
          return convex2 ? s1 < s2 : false;
        else 
          return convex2 ? true : s2 > s1;
    }
  };

  typedef std::priority_queue< HalfEdgeHandle,
                               std::vector<HalfEdgeHandle>,
                               HalfEdgeHandleComparator > HEPriorityQueue;
  bool expandByUnionFind()
  {
   
    unsigned int nb_join = 0;
    unsigned int nb_edges = 0;
    MarkHE toJoin;
    HEPriorityQueue toJoinPQ( HalfEdgeHandleComparator( *this ) );
    std::cout << "Searching concavities: " << std::flush;
    for ( EdgeIterator it = myP.edges_begin(), ite = myP.edges_end();
	  it != ite; ++it, ++nb_edges )
      {
	HalfEdgeHandle he = it;
	HalfEdgeHandle rep1 = myDisjointSets.find_set( smallestHalfEdge( he ) );
	HalfEdgeHandle rep2 = myDisjointSets.find_set( smallestHalfEdge( he->opposite() ) );
	if ( rep1 != rep2 ) 
          {
            toJoin.insert( he );
            toJoinPQ.push( he ); 
          }
      }
    std::cout << toJoinPQ.size() << " found." << std::endl;
    while ( ! toJoinPQ.empty() )
      {
	HalfEdgeHandle he = toJoinPQ.top(); 
	toJoinPQ.pop(); 
        typename MarkHE::iterator it_he = toJoin.find( he );
        if ( it_he == toJoin.end() ) // if not found, invalid half edge.
          continue;
        toJoin.erase( it_he );
	HalfEdgeHandle rep1 = myDisjointSets.find_set( smallestHalfEdge( he ) );
	HalfEdgeHandle rep2 = myDisjointSets.find_set( smallestHalfEdge( he->opposite() ) );
	if ( rep1 == rep2 ) continue;
	// Check fusion
	HalfEdgeHandle he_face1 = smallestHalfEdge( rep1 );
	HalfEdgeHandle he_face2 = smallestHalfEdge( rep2 );
	ASSERT( myPlaneMap.find( he_face1 ) != myPlaneMap.end() );
	ASSERT( myPlaneMap.find( he_face2 ) != myPlaneMap.end() );
	FacetInformation & info1 = myPlaneMap.find( he_face1 )->second;
	FacetInformation & info2 = myPlaneMap.find( he_face2 )->second;
	ASSERT( ! info1.plane.empty() );
	ASSERT( ! info2.plane.empty() );
	if ( info1.plane.size() == 1 ) continue; // invalid initial plane
	if ( info2.plane.size() == 1 ) continue; // invalid initial plane
	// std::cout << "---------------------------------------------------" << std::endl;
	// std::cout << "Plane 1: " << info1.plane << std::endl;
	// std::cout << "Plane 2: " << info2.plane << std::endl;
	DigitalPlane planeTester = info1.plane;
	bool ok = planeTester.extend( info2.plane.begin(), info2.plane.end() );
	if ( ! ok ) continue;
	removeFacet( toJoin, he_face1 );
	removeFacet( toJoin, he_face2 );
	myPlaneMap.erase( he_face1 );
	myPlaneMap.erase( he_face2 );
	myDisjointSets.link( rep1, rep2 );
	HalfEdgeHandle nrep = myDisjointSets.find_set( rep1 );
	HalfEdgeHandle nhe_face = smallestHalfEdge( nrep );
	typename HE2Info::iterator itnew = myPlaneMap.find( nhe_face );
	if ( itnew == myPlaneMap.end() )
	  {
	    FacetInformation info;
	    info.plane = planeTester;
	    myPlaneMap[ nhe_face ] = info;
	  }
	else
	  {
	    DGtal::trace.warning() << "--- Weird: faces has already a plane." 
				   << std::endl;
	    itnew->second.plane = planeTester;
	  }
	++nb_join;
      }
    std::cout << "**********************************************************************"
	      << std::endl;
    std::cout << "Joining facets: " << nb_join << " / " << nb_edges << std::endl;
    return nb_join != 0;
  }  



  double signAngle( HalfEdgeConstHandle he ) const
  {
    DRealVector vn1, vn2;
    HalfEdgeConstHandle he_face1 = smallestHalfEdge( he );
    HalfEdgeConstHandle he_face2 = smallestHalfEdge( he->opposite() );
    const FacetInformation & info1 = myPlaneMap.find( he_face1 )->second;
    const FacetInformation & info2 = myPlaneMap.find( he_face2 )->second;
    vn1 = info1.plane.primitive().normal();
    vn2 = info2.plane.primitive().normal();
    Vector n1( vn1[ 0 ], vn1[ 1 ], vn1[ 2 ] );
    Vector n2( vn2[ 0 ], vn2[ 1 ], vn2[ 2 ] );
    Vector c = cross( n1, n2 );
    return c * ( he->vertex()->point() - he->opposite()->vertex()->point() ); 
    // positive is convex
  }

  double signAngle( HalfEdgeHandle he1, HalfEdgeHandle he2 )
  {
    DRealVector vn1, vn2;
    HalfEdgeHandle he_face1 = smallestHalfEdge( he1 );
    HalfEdgeHandle he_face2 = smallestHalfEdge( he2 );
    const FacetInformation & info1 = myPlaneMap.find( he_face1 )->second;
    const FacetInformation & info2 = myPlaneMap.find( he_face2 )->second;
    vn1 = info1.plane.primitive().normal();
    vn2 = info2.plane.primitive().normal();
    Vector n1( vn1[ 0 ], vn1[ 1 ], vn1[ 2 ] );
    Vector n2( vn2[ 0 ], vn2[ 1 ], vn2[ 2 ] );
    Vector c = cross( n1, n2 );
    double sign1 = c * ( he1->vertex()->point() - he1->opposite()->vertex()->point() ); 
    double sign2 = c * ( he2->vertex()->point() - he2->opposite()->vertex()->point() ); 
    if ( ( ( sign1 >= 0.0 ) && ( sign2 > 0.0 ) ) 
         || ( ( sign1 > 0.0 ) && ( sign2 >= 0.0 ) ) 
         || ( ( sign1 <= 0.0 ) && ( sign2 < 0.0 ) ) 
         || ( ( sign1 < 0.0 ) && ( sign2 <= 0.0 ) ) )
      return sign1+sign2;
    else
      return 0.0;
    // positive is convex
  }

  inline Vector cross( const Vector & a, const Vector & b ) const
  {
    return Vector( a.y() * b.z() - a.z() * b.y(),
		   a.z() * b.x() - a.x() * b.z(),
		   a.x() * b.y() - a.y() * b.x() );
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
	DGtal::Color col( random() % 256, random() % 256, random() % 256, 255 );

	// HalfedgeFacetConstCirculator ci = it->facet_begin();
	HalfEdgeHandle he = f.halfedge(); // *ci; //->halfedge();
	Point a( he->vertex()->point() );
	for ( HalfEdgeHandle he2 = he->next(); he2->next() != he; he2 = he2->next() )
	  {
	    Point b( he2->vertex()->point() );
	    Point c( he2->next()->vertex()->point() );
	    displayTriangle( viewer, toDGtal, a, b, c,
			     col, retract );
	  }
      }
  }

  template <typename Viewer>
  void viewUnionFind( Viewer & viewer, double retract )
  {
    DGtal::Color colBasicFacet1( 0, 255, 0, 255 );
    for ( FacetConstIterator it = P().facets_begin(), ite = P().facets_end();
	  it != ite; ++it )
      {
	Facet f = *it;
	HalfEdgeHandle he = f.halfedge();
	DGtal::Color col( random() % 256, random() % 256, random() % 256, 255 );
	myColormap[ myDisjointSets.find_set( smallestHalfEdge( he ) ) ] = col;
      }
    for ( FacetConstIterator it = P().facets_begin(), ite = P().facets_end();
	  it != ite; ++it )
      {
	Facet f = *it;
	HalfEdgeHandle he = f.halfedge();
	DGtal::Color col = myColormap[ myDisjointSets.find_set( smallestHalfEdge( he ) ) ];
	Point a( he->vertex()->point() );
	Point b( he->next()->vertex()->point() );
	Point c( he->next()->next()->vertex()->point() );
	displayTriangle( viewer, toDGtal, a, b, c,
			 col, retract );
      }
  }

  typedef typename std::set<HalfEdgeHandle> Representatives;
  typedef typename std::map<Point,Representatives> MapPoint2Planes;

  template <typename Viewer>
  void viewUnionFindPoints( Viewer & viewer, double retract )
  {
    MapPoint2Planes p2p;
    for ( FacetConstIterator it = P().facets_begin(), ite = P().facets_end();
	  it != ite; ++it )
      {
	Facet f = *it;
	HalfEdgeHandle he = f.halfedge();
	HalfEdgeHandle rep1 = myDisjointSets.find_set( smallestHalfEdge( he ) );
	HalfEdgeHandle he_face1 = smallestHalfEdge( rep1 );
	ASSERT( myPlaneMap.find( he_face1 ) != myPlaneMap.end() );
	FacetInformation & info1 = myPlaneMap.find( he_face1 )->second;
	for ( typename DigitalPlane::ConstIterator itp = info1.plane.begin(),
		itpend = info1.plane.end(); itp != itpend; ++itp )
	  {
	    Point p = toCGAL( *itp );
	    p2p[ p ].insert( rep1 );
	  }
      }
    for ( typename MapPoint2Planes::const_iterator it = p2p.begin(), ite = p2p.end();
	  it != ite; ++it )
      {
	Point p = it->first; 
	const Representatives & reps = it->second;
	DGtal::Color c = pointColor( reps.begin(), reps.end() );
	displayPoint( viewer, toDGtal, p, c, retract );
      }
  }

  template <typename HalfEdgeHandleIterator>
  DGtal::Color pointColor( HalfEdgeHandleIterator itb, HalfEdgeHandleIterator ite )
  {
    int r = 0, g = 0, b = 0;
    int n = 0;
    for ( ; itb != ite; ++itb, ++n )
      {
	DGtal::Color c = myColormap[ *itb ];
	r += c.red();
	g += c.green();
	b += c.blue();
      }
    return DGtal::Color( r / n, g / n, b / n );
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

  typedef KhalimskySpaceND<3,DGtal::int32_t> K3;
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
  typedef DigCore::IndexList                  IndexList;
  typedef DigCore::Index                      Index;

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

  Viewer3D<> viewerCore( ks );
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

  viewerCore << Viewer3D<>::updateDisplay;
  application.exec();

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

  // // Display facet points
  // i = 0;
  // for ( FacetSet::const_iterator it = core.boundary().begin(), 
  // 	  itend = core.boundary().end(); it != itend; ++it, ++i )
  //   {
  //     std::cout << "Facet " << i << ":";
  //     IndexList l = core.getPointIndices( *it );
  //     for ( unsigned int j = 0; j < l.size(); ++j )
  // 	std::cout << " " << l[ j ];
  //     std::cout << std::endl;
  //   }

  
  while ( core.pruneBoundary() ) 
    ;

  // Init polyhedron hull.
  PolyhedronHull<K> PH;
  PH.initFromCore( core );
  // PH.expand();
  while ( PH.expandByUnionFind() ) // ( PH.expandByVertices() )
    {
      Viewer3D<> viewer3d( ks );
      viewer3d.show();
      // PH.view( viewer3d, retract );
      PH.viewUnionFind( viewer3d, retract );
      viewer3d << Viewer3D<>::updateDisplay;
      application.exec();
      // PH.expand();
    }
  
  {
    Viewer3D<> viewer3d( ks );
    viewer3d.show();
    PH.viewUnionFind( viewer3d, retract );
    viewer3d << Viewer3D<>::updateDisplay;
    Viewer3D<> viewer3d_2( ks );
    viewer3d_2.show();
    PH.viewUnionFindPoints( viewer3d_2, retract );
    viewer3d_2 << Viewer3D<>::updateDisplay;
    application.exec();
  }

  std::cout << "number of vertices :  " ;
  std::cout << T.number_of_vertices() << std::endl;
  std::cout << "number of edges :  " ;
  std::cout << T.number_of_edges() << std::endl;
  std::cout << "number of facets :  " ;
  std::cout << T.number_of_facets() << std::endl;
  std::cout << "number of cells :  " ;
  std::cout << T.number_of_cells() << std::endl;


  Point3 low = toCGAL( ks.lowerBound() );
  Point3 high = toCGAL( ks.upperBound() );
  std::ofstream file_off( "bdry.off" );
  // CGAL::Geomview_stream gv( CGAL::Bbox_3( low.x(), low.y(), low.z(), high.x(), high.y(), high.z() ) );
  // gv << P;
  // std::cout << "Enter a key to finish" << std::endl;
  // char ch;
  // std::cin >> ch;
  file_off << PH.P();
  file_off.close();

  // Build_triangle<HalfedgeDS> triangle;
  // P.delegate( triangle);
  // CGAL_assertion( P.is_triangle( P.halfedges_begin()));
  return 0;
}

/*
  A few nice examples:

  ./3d-triangulation-3 -p "(15*x+21*y+37*z-6)*(31*x-18*y-28*z-5)" -b 10 -g 0.25 -G 2
  ./3d-triangulation-3 -p "x^2+y^2+2*z^2-x*y*z+z^3-100" -b 20 -g 0.5 -G 1
  ./3d-triangulation-3 -p "0.5*(z^2-4*4)^2+(x^2-7*7)^2+(y^2-7*7)^2-7.2*7.2*7.2*7.2" -b 15 -g 0.5 -G 1
  ./3d-triangulation-3 -p "(1-(x^2+y^2))^2*(1+(x^2+y^2))^2-z" -b 2 -g 0.05 -G 1

 */
