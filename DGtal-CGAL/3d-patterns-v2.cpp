#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <QtGui/qapplication.h>

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
#include <CGAL/Cartesian.h>
#include <CGAL/CORE/Expr.h>

#include "Auxiliary.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Cartesian<CORE::Expr> K;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
//typedef Delaunay::Vertex_circulator Vertex_circulator;
//typedef Delaunay::Edge_iterator  Edge_iterator;
typedef Delaunay::Cell_iterator     Cell_iterator;
typedef Delaunay::Facet             Facet;
typedef Delaunay::Facet_iterator    Facet_iterator;
typedef Delaunay::Edge              Edge;
typedef Delaunay::Edge_iterator     Edge_iterator;
typedef Delaunay::Vertex_iterator   Vertex_iterator;
typedef Delaunay::Vertex_handle     Vertex_handle;
typedef Delaunay::Point             CGALPoint;
typedef Delaunay::Cell_handle       Cell_handle;
typedef Delaunay::Triangle          Triangle;
typedef DGtal::SpaceND<3, DGtal::int32_t> Z3;
typedef Z3::Point Point;
typedef Z3::RealPoint RealPoint;
typedef DGtal::HyperRectDomain<Z3> Domain;
typedef Domain::ConstIterator DomainConstIterator;

Point toDGtal(const CGALPoint &p)
{
  return Point ( p.x(),
		 p.y(),
		 p.z() );
}
RealPoint toDGtalR(const CGALPoint &p)
{
  return RealPoint ( p.x(),
                     p.y(),
                     p.z() );
}
// Point toDGtal(const CGALPoint &p)
// {
//   return Point ( p.x().longValue(),
// 		 p.y().longValue(),
// 		 p.z().longValue() );
// }
CGALPoint toCGAL(const Point &p)
{
  return CGALPoint( p[0], p[1], p[2] );
}

inline
Point operator^( const Point & a, const Point & b )
{
  BOOST_STATIC_ASSERT( Point::dimension == 3 );
  return Point( a[ 1 ] * b[ 2 ] - a[ 2 ] * b[ 1 ],
		a[ 2 ] * b[ 0 ] - a[ 0 ] * b[ 2 ],
		a[ 0 ] * b[ 1 ] - a[ 1 ] * b[ 0 ] );
}


/**
   Handy structure to hold a CGAL edge of a Triangulaion_3. It is
   unambiguous (meaning an edge has only one representant) and permits
   comparison (for use as index in std::map).
*/
struct VEdge {
public:
  Vertex_handle first;
  Vertex_handle second;
  inline VEdge( const Edge & e )
  {
    first = (e.first)->vertex( e.second );
    second = (e.first)->vertex( e.third );
    if ( second < first ) std::swap( first, second );
  }
  inline VEdge( Vertex_handle v1, Vertex_handle v2 )
  {
    if ( v1 < v2 ) { first = v1; second = v2; }
    else           { first = v2; second = v1; }
  }
  bool operator<( const VEdge & other ) const
  {
    return ( first < other.first )
      || ( ( first == other.first )
           && ( second < other.second ) );
  }
};

/**
   Handy structure to hold a CGAL facet of a Triangulaion_3. It is
   unambiguous (meaning a facet has only one representant) and permits
   comparison (for use as index in std::map).
*/
struct VFacet {
public:
  Vertex_handle first;
  Vertex_handle second;
  Vertex_handle third;
  inline VFacet( const Facet & f )
  {
    int i = f.second;
    first = (f.first)->vertex( (i+1)%4 );
    second = (f.first)->vertex( (i+2)%4 );
    third = (f.first)->vertex( (i+3)%4 );
    sort();
  }
  inline VFacet( Vertex_handle v1, Vertex_handle v2, Vertex_handle v3 )
    : first( v1 ), second( v2 ), third( v3 )
  {
    sort();
  }
  
  inline void sort()
  {
    if ( second < first ) std::swap( first, second );
    if ( third < first )  std::swap( first, third );
    if ( third < second ) std::swap( second, third );
  }
  bool operator<( const VFacet & other ) const
  {
    return ( first < other.first )
      || ( ( first == other.first )
	   && ( ( second < other.second ) 
		|| ( ( second == other.second )  
		     && ( third < other.third ) ) ) );
  }
};

struct OFacet {
  Facet facet;
  double criterion;

  inline OFacet() : criterion( 0 ) {}

  inline OFacet( const Facet & f )
    : facet( f ) 
  {
    computeCriterion();
  }
  inline OFacet( const Facet & f, DGtal::int64_t nl1 )
    : facet( f ), criterion( nl1 )
  {}

  inline OFacet( const OFacet & other )
    : facet( other.facet ), criterion( other.criterion )
  {}

  inline
  OFacet & operator=( const OFacet & other )
  {
    if ( this != &other )
      {
	facet = other.facet;
	criterion = other.criterion;
      }
    return *this;
  }

  void computeCriterion()
  {
    criterion = 0.0;
    for ( int i = 0; i < 4; ++i )
      {
	Facet g = facet;
	g.second = i;
	VFacet vf( g );
	criterion = std::max( criterion, crossNormL1( vf ) );
      }
  }
  double crossNormL1( const VFacet & vf)
  {
    Point a( toDGtal( vf.first->point() ) );
    Point b( toDGtal( vf.second->point() ) );
    Point c( toDGtal( vf.third->point() ) );
    Point n( (b-a)^(c-a) );
    return n.norm( Point::L_1 );
  }
};

struct OFacetLessComparator{
  inline
  bool operator()( const OFacet & f1, const OFacet & f2 ) const 
  {
    return ( f1.criterion < f2.criterion )
      || ( ( f1.criterion == f2.criterion )
	   && ( f1.facet < f2.facet ) );
  }
};


///////////////////////////////////////////////////////////////////////////////
// OrderedCell
///////////////////////////////////////////////////////////////////////////////
struct OrderedCell {
  Cell_handle cell;
  double criterion;

  inline OrderedCell() : criterion( 0 ) {}

  inline OrderedCell( const Cell_handle & f )
    : cell( f ) 
  {
    computeCriterion();
  }
  inline OrderedCell( const Cell_handle & f, DGtal::int64_t nl1 )
    : cell( f ), criterion( nl1 )
  {}

  inline OrderedCell( const OrderedCell & other )
    : cell( other.cell ), criterion( other.criterion )
  {}

  inline
  OrderedCell & operator=( const OrderedCell & other )
  {
    if ( this != &other )
      {
	cell = other.cell;
	criterion = other.criterion;
      }
    return *this;
  }

  void computeCriterion()
  {
    criterion = 0.0;
    for ( int i = 0; i < 4; ++i )
      {
	criterion = std::max( criterion, crossNormL1( Facet( cell, i ) ) );
      }
  }

  double crossNormL1( const Facet & f )
  {
    int i = f.second;
    Point a( toDGtal( f.first->vertex( (i+1)%4 )->point() ) );
    Point b( toDGtal( f.first->vertex( (i+2)%4 )->point() ) );
    Point c( toDGtal( f.first->vertex( (i+3)%4 )->point() ) );
    Point n( (b-a)^(c-a) );
    return n.norm( Point::L_1 );
  }
};

struct OrderedCellLessComparator{
  inline
  bool operator()( const OrderedCell & f1, const OrderedCell & f2 ) const 
  {
    return ( f1.criterion < f2.criterion )
      || ( ( f1.criterion == f2.criterion )
	   && ( f1.cell < f2.cell ) );
  }
};


typedef std::map< Cell_iterator, DGtal::int64_t > MapCell2Int;
typedef std::map< VFacet, DGtal::int64_t > MapFacet2Int;
typedef std::map< VEdge, DGtal::int64_t > MapEdge2Int;
typedef std::map< Vertex_iterator, DGtal::int64_t > MapVertex2Int;
typedef std::set< Facet > FacetSet;
typedef std::set< OFacet, OFacetLessComparator > OFacetSet;

DGtal::uint64_t 
markBasicEdges( MapEdge2Int & basicMap, const Delaunay & t )
{
  DGtal::uint64_t nb = 0;
  for( Edge_iterator it = t.edges_begin(), itend = t.edges_end();
       it != itend; ++it)
    {
      VEdge edge( *it );
      if ( ! t.is_infinite( *it ) )
        {
          Point a( toDGtal( edge.first->point() ) );
          Point b( toDGtal( edge.second->point() ) );
          if ( (b-a).norm( Point::L_infty ) == 1 )
            {
              basicMap[ edge ] = 1;
              ++nb;
            }
          else
            {
              basicMap[ edge ] = 0;
            }
        }
      else
        {
          basicMap[ edge ] = -1;
        }
    }
  return nb;
}

/**
   @return true if the facet is planar and elementary.
*/
bool 
checkPlanarFacet( const Delaunay & t, const Facet & f )
{
  if ( ! t.is_infinite( f ) )
    { 
      VFacet vf( f );
      Point a( toDGtal( vf.first->point() ) );
      Point b( toDGtal( vf.second->point() ) );
      Point c( toDGtal( vf.third->point() ) );
      Point n( (b-a)^(c-a) );
      return n.norm( Point::L_infty ) == 1;
    }
  return false;
}
/**
   @return true if the facet is planar and elementary.
*/
bool 
checkPlanarVFacet( const VFacet & vf )
{
  Point a( toDGtal( vf.first->point() ) );
  Point b( toDGtal( vf.second->point() ) );
  Point c( toDGtal( vf.third->point() ) );
  Point n( (b-a)^(c-a) );
  return n.norm( Point::L_infty ) == 1;
}


DGtal::uint64_t 
markBasicFacets( FacetSet & bFacets, OFacetSet & qFacets,
		 const Delaunay & t, MapEdge2Int & basicEdgeMap )
{
  DGtal::uint64_t nb = 0;
  for( Facet_iterator it = t.facets_begin(), itend = t.facets_end();
       it != itend; ++it)
    {
      if ( ! t.is_infinite( *it ) )
        {
	  VFacet f( *it ); 
          Cell_iterator itCell = it->first; int i = it->second;
          VEdge e1( f.first, f.second );
          VEdge e2( f.second, f.third );
          VEdge e3( f.third, f.first );
          unsigned int n = 0;
          n += basicEdgeMap[ e1 ] == 1 ? 1 : 0;
          n += basicEdgeMap[ e2 ] == 1 ? 1 : 0;
          n += basicEdgeMap[ e3 ] == 1 ? 1 : 0;
          if ( ( n == 3 ) ) // || checkPlanarVFacet( f ) )
            {
	      OFacet f1( *it );
	      OFacet f2( t.mirror_facet( f1.facet ), f1.criterion );
	      bFacets.insert( f1.facet );
	      bFacets.insert( f2.facet );
	      qFacets.insert( f1 );
	      qFacets.insert( f2 );
              ++nb;
            }
        }
    }
  return nb;
}

DGtal::int64_t
countLatticePointsInTetrahedra( const Point & a, const Point & b, const Point & c, const Point & d )
{
  DGtal::int64_t nb = 0;
  Point ab = b - a;
  Point bc = c - b;
  Point cd = d - c;  
  Point da = a - d;
  Point abc = ab^bc;
  Point bcd = bc^cd;
  Point cda = cd^da;
  Point dab = da^ab;
  Point::Component abc_shift = abc.dot( a );
  Point::Component bcd_shift = bcd.dot( b );
  Point::Component cda_shift = cda.dot( c );
  Point::Component dab_shift = dab.dot( d );
  if ( abc.dot( d ) < abc_shift )  { abc = -abc; abc_shift = -abc_shift; }
  if ( bcd.dot( a ) < bcd_shift )  { bcd = -bcd; bcd_shift = -bcd_shift; }
  if ( cda.dot( b ) < cda_shift )  { cda = -cda; cda_shift = -cda_shift; }
  if ( dab.dot( c ) < dab_shift )  { dab = -dab; dab_shift = -dab_shift; }
  Point inf = a.inf( b ).inf( c ).inf( d );
  Point sup = a.sup( b ).sup( c ).sup( d );
  Domain domain( inf, sup );
  for ( DomainConstIterator it = domain.begin(), itE = domain.end();
	it != itE; ++it )
    {
      Point p = *it;
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
countLatticePointsInTetrahedraIf4( const Point & a, const Point & b, const Point & c, const Point & d )
{
  DGtal::int64_t nb = 0;
  Point ab = b - a;
  Point bc = c - b;
  Point cd = d - c;  
  Point da = a - d;
  Point abc = ab^bc;
  Point bcd = bc^cd;
  Point cda = cd^da;
  Point dab = da^ab;
  Point::Component abc_shift = abc.dot( a );
  Point::Component bcd_shift = bcd.dot( b );
  Point::Component cda_shift = cda.dot( c );
  Point::Component dab_shift = dab.dot( d );
  if ( abc.dot( d ) < abc_shift )  { abc = -abc; abc_shift = -abc_shift; }
  if ( bcd.dot( a ) < bcd_shift )  { bcd = -bcd; bcd_shift = -bcd_shift; }
  if ( cda.dot( b ) < cda_shift )  { cda = -cda; cda_shift = -cda_shift; }
  if ( dab.dot( c ) < dab_shift )  { dab = -dab; dab_shift = -dab_shift; }
  Point inf = a.inf( b ).inf( c ).inf( d );
  Point sup = a.sup( b ).sup( c ).sup( d );
  Domain domain( inf, sup );
  for ( DomainConstIterator it = domain.begin(), itE = domain.end();
	it != itE; ++it )
    {
      Point p = *it;
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

void 
getFacets( std::vector<Facet> & facets, const Cell_iterator & it )
{
  for ( int i = 0; i < 4; ++i )
    facets.push_back( Facet( it, i ) );
}

double
computeDihedralAngle( const Delaunay & t, const Facet & f1, const Facet & f2 )
{
  typedef CGAL::Vector_3<K> Vector;
  ASSERT( f1.first == f2.first );
  Vector n1 = t.triangle(f1).supporting_plane().orthogonal_vector();
  n1 = n1 / sqrt( n1.squared_length() );
  Vector n2 = t.triangle(f2).supporting_plane().orthogonal_vector();
  n2 = n2 / sqrt( n2.squared_length() );
  return acos( (double) (n1*n2) );
}

namespace DGtal {
  class DigitalCore {
  public:
    typedef Delaunay                      Triangulation;
    typedef Triangulation::Cell_iterator  CellIterator;
    typedef Triangulation::Cell_handle    CellHandle;
    typedef Triangulation::Facet          Facet;
    typedef Triangulation::Facet_iterator FacetIterator;
    typedef Triangulation::Vertex_handle  VertexHandle;
    typedef DGtal::uint64_t               Size;
    typedef std::map< CellHandle, Size >  MapCell2Size;
    typedef std::set< Facet >             FacetSet;
    typedef std::set< Cell_handle >       CellSet;
    typedef std::set< OrderedCell, 
                      OrderedCellLessComparator > OrderedCellSet;
    
  private:
    /// The Delaunay triangulation that carries the digital core.
    const Triangulation* myTriangulation;
    /// Area factor when extending the boundary.
    double myAreaFactorExtend;
    /// Area factor when retracting the boundary.
    double myAreaFactorRetract;
    /// Stores the number of lattice of points within each triangulation Cell.
    MapCell2Size myNbLatticePoints;
    /// Current boundary as a set of facets.
    FacetSet myBoundary;
    /// Current interior cells.
    CellSet myInterior;
    //--------------------------------------------------------------------------
  public:

    inline 
    DigitalCore() : myTriangulation( 0 ) 
    {}
    
    inline 
    DigitalCore( ConstAlias<Triangulation> t ) 
      : myTriangulation( &t )
    {
      countLatticePoints();
      computeBasicFacets();
    }

    //--------------------------------------------------------------------------
  public:

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
      if ( ! myTriangulation->is_infinite( f ) )
        {
          int i = f.second;
          Point a( toDGtal( f.first->vertex( (i+1)%4 )->point() ) );
          Point b( toDGtal( f.first->vertex( (i+2)%4 )->point() ) );
          Point c( toDGtal( f.first->vertex( (i+3)%4 )->point() ) );
          return ( ( (b-a).norm( Point::L_infty ) == 1 )
                   && ( (c-b).norm( Point::L_infty ) == 1 )
                   && ( (a-c).norm( Point::L_infty ) == 1 ) );
        }
      else
        return false;
    }

    inline 
    Size nbLatticePoints( const CellHandle & cell ) const
    {
      MapCell2Size::const_iterator it = myNbLatticePoints.find( cell );
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
      myTriangulation->incident_cells( vh, std::back_inserter( incident_cells ) );
      for ( std::vector<CellHandle>::const_iterator it = incident_cells.begin(),
              itend = incident_cells.end(); it != itend; ++it )
        {
          if ( ! isCellInterior( *it ) ) return false;
        }
      return true;
    }

    inline
    double area( const Facet & f ) const
    {
      return sqrt( myTriangulation->triangle( f ).squared_area() );
    }

    void extend( double area_factor_extend = 1.733 )
    {
      ASSERT( myTriangulation != 0 );
      trace.beginBlock("[DigitalCore] Extend boundary.");
      myAreaFactorExtend = area_factor_extend;

      // Queue for computing the digital core.
      OrderedCellSet priorityQ;
      
      // Prepare queue.
      for ( FacetSet::const_iterator it = myBoundary.begin(), itend = myBoundary.end();
            it != itend; ++it )
        priorityQ.insert( OrderedCell( it->first ) );

      // Start extension
      bool isMarked[ 4 ];
      unsigned int nb = 0;
      while ( ! priorityQ.empty() )
        {
          // trace.progressBar( 1.0, 1.0+log( (double) priorityQ.size() ) );
          OrderedCellSet::iterator it = priorityQ.begin();
          CellHandle cell_h = it->cell;
          priorityQ.erase( it );
          // infinite cells are not processed.
          if ( myTriangulation->is_infinite( cell_h ) ) continue;
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

    void retract( double area_factor_retract = 1.0 )
    {
      ASSERT( myTriangulation != 0 );
      trace.beginBlock("[DigitalCore] Retract boundary.");
      myAreaFactorRetract = area_factor_retract;

      // Queue for computing the digital core.
      OrderedCellSet priorityQ;
      CellSet processedCells;
      // Prepare queue.
      for ( FacetSet::const_iterator it = myBoundary.begin(), itend = myBoundary.end();
            it != itend; ++it )
        {
          Facet mf = myTriangulation->mirror_facet( *it );
          priorityQ.insert( OrderedCell( mf.first ) );
        }

      // Start retraction
      bool isMarked[ 4 ];
      unsigned int nb = 0;
      while ( ! priorityQ.empty() )
        {
          // trace.progressBar( 1.0, 1.0+log( (double) priorityQ.size() ) );
          OrderedCellSet::iterator it = priorityQ.begin();
          CellHandle cell_h = it->cell;
          priorityQ.erase( it );
          if ( processedCells.find( cell_h ) != processedCells.end() )
            continue;
          // // infinite cells are not processed.
          // if ( myTriangulation->is_infinite( cell_h ) ) continue;
          // // cells containing integer points in the interior are not processed.
          // if ( nbLatticePoints( cell_h ) != 4 )         continue;
          // checking cell faces for retraction.
          if ( checkCellRetraction( isMarked, cell_h ) )
            {
              retractCell( priorityQ, isMarked, cell_h );
              ++nb;
            }
          processedCells.insert( cell_h );
        }
      trace.info() << "- cells flipped = " << nb << std::endl;
      trace.info() << "- boundary has " << myBoundary.size() << " facets." << std::endl;
      trace.endBlock();
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
      bool propagate = n >= 3;
      if ( n == 2 )
        { // 2 are marked. We check their area.
          double area_boundary = 0.0;
          double area_exterior = 0.0;
          for ( int i = 0; i < 4; ++i ) {
            if ( isBoundary[ i ] ) 
              area_boundary += area( Facet( cell, i ) );
            else
              area_exterior += area( Facet( cell, i ) );
          }
          if ( area_exterior < myAreaFactorExtend * area_boundary )
            propagate = true;
        }
      return propagate;
    }
    
    void extendCell( OrderedCellSet & priorityQ, bool isBoundary[ 4 ], const CellHandle & cell )
    { // Switch cell to interior
      myInterior.insert( cell );
      for ( unsigned int i = 0; i < 4; ++i ) {
        if ( ! isBoundary[ i ] )
          {
            Facet nfacet = myTriangulation->mirror_facet( Facet( cell, i ) );
            priorityQ.insert( nfacet.first );
            myBoundary.insert( nfacet );
          }
        else
          {
            myBoundary.erase( Facet( cell, i ) );
          }
      }
    }

    bool checkCellRetraction( bool isBoundary[ 4 ], const CellHandle & cell ) const
    {
      unsigned int n = 0;
      std::vector<Facet> mirror_facets;
      for ( int i = 0; i < 4; ++i )
        {
          mirror_facets.push_back( myTriangulation->mirror_facet( Facet( cell, i ) ) ); 
          isBoundary[ i ] = isFacetBoundary( mirror_facets[ i ] );
          if ( isBoundary[ i ] ) ++n;
        }
      bool retract = false;
      if ( n == 2 )
        {
          // We check their area.
          double area_boundary = 0.0;
          double area_interior = 0.0;
          for ( int i = 0; i < 4; ++i ) {
            if ( isBoundary[ i ] ) 
              area_boundary += area( mirror_facets[ i ] );
            else
              area_interior += area( mirror_facets[ i ] );
          }
          if ( area_interior < myAreaFactorRetract * area_boundary )
            retract = true;
        }
      return retract;
    }
    
    void retractCell( OrderedCellSet & priorityQ, bool isBoundary[ 4 ], const CellHandle & cell )
    { // Switch cell to exterior
      myInterior.erase( cell );
      for ( unsigned int i = 0; i < 4; ++i ) {
        if ( ! isBoundary[ i ] )
          {
            Facet nfacet = Facet( cell, i );
            priorityQ.insert( myTriangulation->mirror_facet( nfacet ).first );
            myBoundary.insert( nfacet );
          }
        else
          {
            myBoundary.erase( myTriangulation->mirror_facet( Facet( cell, i ) ) );
          }
      }
    }



    void computeBasicFacets()
    {
      ASSERT( myTriangulation != 0 );
      trace.beginBlock("[DigitalCore] Compute basic facets.");
      DGtal::uint64_t nb = 0;
      for( FacetIterator it = myTriangulation->facets_begin(), itend = myTriangulation->facets_end();
           it != itend; ++it )
        {
          if ( isFacetBasic( *it ) ) // || checkPlanarVFacet( f ) )
            {
	      myBoundary.insert( *it );
	      myBoundary.insert( myTriangulation->mirror_facet( *it ) );
              ++nb;
            }
        }
      trace.info() << "- nb basic facets = " << nb << std::endl;
      trace.endBlock();
    }

    void countLatticePoints()
    {
      ASSERT( myTriangulation != 0 );
      trace.beginBlock("[DigitalCore] Counting lattice points.");
      double maxbar = myTriangulation->number_of_cells();
      double bar = 0.0;
      int i = 0;
      for( CellIterator it = myTriangulation->cells_begin(), 
             itend = myTriangulation->cells_end();
           it != itend; ++it, ++bar, ++i )
        {
          if ( i % 1000 == 0 ) trace.progressBar( bar, maxbar );
          if ( myTriangulation->is_infinite( it ) ) 
            myNbLatticePoints[ it ] = -1;
          else
            {
              Point a( toDGtal(it->vertex(0)->point())),
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
int main (int argc, char** argv )
{
  using namespace DGtal;

  typedef KhalimskySpaceND<3,DGtal::int32_t> K3;
  typedef Z3::Vector Vector;
  typedef Z3::RealPoint RealPoint;
  typedef K3::SCell SCell;
  typedef K3::SCellSet SCellSet;
  typedef SCellSet::const_iterator SCellSetConstIterator;
  // typedef ImplicitRoundedHyperCube<Z3> Shape;
  typedef RealPoint::Coordinate Ring;

  typedef Viewer3D<Z3,K3> MViewer3D;

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
    ("prune,P", po::value<double>(), "prunes the resulting the complex of the tetrahedra were the length of the dubious edge is [arg] times the length of the opposite edge." )
    ("geometry,G", po::value<int>()->default_value( 0 ), "specifies how digital points are defined from the digital surface: arg=0: inner voxels, arg=1: inner and outer voxels, arg=2: pointels" )
    ("area-factor-extend,A", po::value<double>()->default_value( sqrt(3.0) ), "specifies the convexity/concavity area ratio when extending the core." )
    ("area-factor-retract,a", po::value<double>()->default_value( 1.0/sqrt(3.0) ), "specifies the convexity/concavity area ratio when retracting the core." )
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

  SCellSet boundary2 (boundary);

  trace.beginBlock("Constructing the set of points");
  std::vector<Point> digital_points;
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
  Delaunay t;
  double setsize = (double) digital_points.size()-1;
  trace.info() << "Vertices to process: " << setsize << std::endl;
  double step = 0.0;
  int i = 0;
  for ( std::vector<Point>::const_iterator it = digital_points.begin(), itE = digital_points.end();
	it != itE; ++it, ++step, ++i )
    {
      if ( i % 1000 == 0 ) trace.progressBar( step, setsize );
      t.insert( toCGAL( *it ) );
    }
  trace.endBlock();

  // Testing Core
  double area_factor_extend = vm["area-factor-extend"].as<double>();
  double area_factor_retract = vm["area-factor-retract"].as<double>();
  DigitalCore core( t );
  core.extend( area_factor_extend );
  core.retract( area_factor_retract );

  // trace.beginBlock( "Couting interior and boundary vertices." );
  // Delaunay t2;
  // unsigned int nb_int = 0;
  // unsigned int nb_bd = 0;
  // for ( Vertex_iterator it = t.vertices_begin(), itend = t.vertices_end();
  //       it != itend; ++it )
  //   if ( core.isVertexInterior( it ) ) ++nb_int;
  //   else {
  //     t2.insert( it->point() );
  //     ++nb_bd;
  //   }
  // trace.info() << "- nb interior vertices = " << nb_int << std::endl;
  // trace.info() << "- nb boundary vertices = " << nb_bd << std::endl;
  // trace.endBlock();


  // start viewer
  int view = vm["view"].as<int>();
  MViewer3D viewerCore;
  viewerCore.show();
  Color colBasicFacet2( 0, 255, 255, 255 );
  Color colBasicFacet1( 0, 255, 0, 255 );
  Color colSpuriousTetrahedra( 255, 0, 0, 100 );
  if ( view & 0x1 ) { // View digital core.
    for ( FacetSet::const_iterator it = core.boundary().begin(), itend = core.boundary().end();
          it != itend; ++it )
      {
        // we display it.
        Triangle triangle = t.triangle( *it );
        Point a( toDGtal( triangle.vertex( 0 ) ) );
        Point b( toDGtal( triangle.vertex( 1 ) ) );
        Point c( toDGtal( triangle.vertex( 2 ) ) );
        Facet f2 = t.mirror_facet( *it );
        if ( core.isFacetBoundary( f2 ) )
          { // the mirror facet is also in the triangulation. 
            // We need to move vertices a little bit when two triangles are at the same position.
            Point n = (b-a)^(c-a);
            double norm = n.norm(Point::L_2);
            double dx[ 3 ];
            for ( unsigned int j = 0; j < 3; ++j )
              dx[ j ] = 0.001*((double) n[j])/norm;
            RealPoint A( (double) a[ 0 ] + dx[ 0 ], (double) a[ 1 ] +  dx[ 1 ], (double) a[ 2 ] + dx[ 2 ] );
            RealPoint B( (double) b[ 0 ] + dx[ 0 ], (double) b[ 1 ] +  dx[ 1 ], (double) b[ 2 ] + dx[ 2 ] );
            RealPoint C( (double) c[ 0 ] + dx[ 0 ], (double) c[ 1 ] +  dx[ 1 ], (double) c[ 2 ] + dx[ 2 ] );
            viewerCore.setFillColor( colBasicFacet2 );
            viewerCore.addTriangle( A, C, B );
          }
        else
          {
            RealPoint A( a[ 0 ], a[ 1 ], a[ 2 ] );
            RealPoint B( b[ 0 ], b[ 1 ], b[ 2 ] );
            RealPoint C( c[ 0 ], c[ 1 ], c[ 2 ] );
            viewerCore.setFillColor( colBasicFacet1 );
            viewerCore.addTriangle( A, B, C );
          }
      }
  } //  if ( view & 0x1 ) {

  if ( view & 0x2 ) { // View spurious tetrahedra
    for ( Cell_iterator it = t.cells_begin(), itend = t.cells_end();
          it != itend; ++it )
      {
        bool draw = false;
        Color col;
        if ( ( ! t.is_infinite( it ) )
             && ( core.nbLatticePoints(  it ) == 4 )
             && ( core.interior().find( it ) == core.interior().end() ) )
          {
            draw = true;
            col = colSpuriousTetrahedra;
          }
        if ( draw )
          {
            RealPoint a( toDGtalR(it->vertex(0)->point()));
            RealPoint b( toDGtalR(it->vertex(1)->point()));
            RealPoint c( toDGtalR(it->vertex(2)->point()));
            RealPoint d( toDGtalR(it->vertex(3)->point()));
            viewerCore.setFillColor( col );
            viewerCore.addTriangle( a, b, c );
            viewerCore.addTriangle( c, b, d );
            viewerCore.addTriangle( c, d, a );
            viewerCore.addTriangle( d, b, a );
          }
      }
  }
 
  viewerCore << MViewer3D::updateDisplay;
  application.exec();

  // start viewer
  MViewer3D viewer1, viewer2, viewer3;
  viewer1.show();
  viewer2.show();
  viewer3.show();
  // for ( SCellSetConstIterator it=boundary.begin(), itend=boundary.end(); it != itend; ++it )
  //    viewer1 << *it;
  
  trace.beginBlock("Counting lattice points.");
  MapCell2Int nbLatticePoints;
  for(Cell_iterator it = t.cells_begin(), itend=t.cells_end();
      it != itend; ++it)
    {
      if ( t.is_infinite( it ) ) 
        {
          nbLatticePoints[ it ] = -1;
        }
      else
        {
          Point a( toDGtal(it->vertex(0)->point())),
            b(toDGtal(it->vertex(1)->point())),
            c(toDGtal(it->vertex(2)->point())),
            d(toDGtal(it->vertex(3)->point()));
          nbLatticePoints[ it ] = countLatticePointsInTetrahedraIf4( a, b, c, d );
        }
    }
  trace.endBlock();

  trace.beginBlock("Computing basic edges");
  MapEdge2Int bEdge;
  uint64_t nbBE = markBasicEdges( bEdge, t );
  trace.info() << "Nb basic edges :" << nbBE << std::endl;

  Color colBasicEdge( 0, 0, 255, 255 );
  for( Edge_iterator it = t.edges_begin(), itend = t.edges_end();
       it != itend; ++it)
    {
      VEdge e( *it );
      if ( bEdge[ e ] == 1 )
        {
          Cell_handle itC = it->first; 
          int i = it->second;
          int j = it->third;
          RealPoint a( toDGtalR(itC->vertex( i )->point()));
          RealPoint b( toDGtalR(itC->vertex( j )->point()));
          viewer1.setLineColor( colBasicEdge );
          viewer1.addLine( a, b, 1.0 );
        }
    }
  trace.endBlock();

  trace.beginBlock("Extracting Min-Polyhedron");

  std::set<Cell_handle> markedCells;
  FacetSet basicFacet;
  FacetSet markedFacet;
  OFacetSet priorityQ;
  FacetSet elementQ;
  std::set<Cell_handle> weirdCells;
  uint64_t nbBF = markBasicFacets( basicFacet, priorityQ, t, bEdge );
  trace.info() << "Nb basic facets:" << nbBF << std::endl;

  // At the beginning, the queue is equal to marked facets.
  markedFacet = basicFacet;
  elementQ = markedFacet;
  bool isMarked[ 4 ];
  while ( ! priorityQ.empty() )
    {
      OFacetSet::iterator it = priorityQ.begin();
      OFacet ofacet = *it;
      priorityQ.erase( it );
      elementQ.erase( ofacet.facet );
      // infinite cells are not processed.
      if ( t.is_infinite( ofacet.facet.first ) ) continue;
      std::vector<Facet> facets;
      getFacets( facets, ofacet.facet.first );
      bool found_neighbors_in_queue = false;
      for ( unsigned int i = 0; i < facets.size(); ++i )
        {
          if ( elementQ.find( facets[ i ] ) != elementQ.end() )
            {
              found_neighbors_in_queue = true;
              break;
            }
        }
      // a later facet will take care of the possible fusion.
      if ( found_neighbors_in_queue ) continue;
      if ( nbLatticePoints[ ofacet.facet.first ] == 4 )
        { // potential fusion
          std::vector<unsigned int> f_marked;
          std::vector<unsigned int> f_unmarked;
          for ( unsigned int i = 0; i < facets.size(); ++i )
            {
              isMarked[ i ] = ( markedFacet.find( facets[ i ] ) != markedFacet.end() );
              if ( isMarked[ i ] )  f_marked.push_back( i );
              else                  f_unmarked.push_back( i );
            }
          // // Add basic planar facets when encountered.
          // for ( std::vector<unsigned int>::iterator it_fu = f_unmarked.begin();
          //       it_fu != f_unmarked.end(); )
          //   {
          //     if ( checkPlanarFacet( t, facets[ *it_fu ] ) )
          //       {
          //         Facet mf = t.mirror_facet( facets[ *it_fu ] );
          //         markedFacet.insert( facets[ *it_fu ] );
          //         markedFacet.insert( mf );
          //         isMarked[ *it_fu ] = true;
          //         f_marked.push_back( *it_fu );
          //         f_unmarked.erase( it_fu );
          //       }
          //     else 
          //       ++it_fu;
          //   }
          unsigned int n = f_marked.size();
          // At least 2 are marked, by convexity, we can close the gap
          // to the further faces of the cell.
          // std::cout << " " << n;
          bool propagate = n >= 3;
          if ( n == 2 )
            {
              // We must check that further tetrahedra do not contain integer points.
              // Facet h0 = t.mirror_facet( facets[ f_marked[ 0 ] ] );
              // Facet h1 = t.mirror_facet( facets[ f_marked[ 1 ] ] );
              // std::cout << "(" << nbLatticePoints[ h0.first ] 
              // 		<< "," << nbLatticePoints[ h1.first ] << ")";
              double area_interior = 
                sqrt( t.triangle( facets[ f_marked[ 0 ] ] ).squared_area() )
                + sqrt( t.triangle( facets[ f_marked[ 1 ] ] ).squared_area() );
              double area_exterior = 
                sqrt( t.triangle( facets[ f_unmarked[ 0 ] ] ).squared_area() )
                + sqrt( t.triangle( facets[ f_unmarked[ 1 ] ] ).squared_area() );
              if ( area_exterior > area_factor_extend * area_interior )
                weirdCells.insert( ofacet.facet.first );
              else
                propagate = true;
              // if ( ( nbLatticePoints[ h0.first ] != 4 )
              // 	   && ( nbLatticePoints[ h1.first ] != 4 ) )
              // 	propagate = false;
            }
          if ( propagate ) 
            {
              // marked cell as interior
              // std::cout << "+";
              markedCells.insert( ofacet.facet.first );
              for ( unsigned int i = 0; i < facets.size(); ++i )
                {
                  if ( ! isMarked[ i ] )
                    {
                      OFacet nfacet1( facets[ i ] );
                      OFacet nfacet2( t.mirror_facet( nfacet1.facet ), nfacet1.criterion );
                      // markedFacet.insert( nfacet1.facet );
                      // put everything into queue.
                      priorityQ.insert( nfacet2 );
                      markedFacet.insert( nfacet2.facet );
                      elementQ.insert( nfacet2.facet );
                    }
                }
            }
        }
    }
  std::cout << std::endl;
  trace.endBlock();

  trace.beginBlock("Prune weird cells");
  trace.info() << "weird cells: " << weirdCells.size() << std::endl;
  std::set<Cell_handle> removedWeirdCells;
  bool prune = vm.count( "prune" );
  if ( prune )
    {
      double factor = vm[ "prune" ].as<double>();
      double changed = true;
      while ( changed ) {
        changed = false;
        for ( std::set<Cell_handle>::const_iterator it = weirdCells.begin(), itend = weirdCells.end();
              it != itend; ++it )
          {
            if ( removedWeirdCells.find( *it ) == removedWeirdCells.end() ) 
              { 
                std::vector<Facet> facets;
                getFacets( facets, *it );
                std::vector<unsigned int> f_exterior;
                std::vector<unsigned int> f_interior;
                std::vector<unsigned int> f_basic;
                for ( unsigned int i = 0; i < facets.size(); ++i )
                  { // count exterior facets
                    Facet m = t.mirror_facet( facets[ i ] );
                    if ( ( markedFacet.find( m ) != markedFacet.end() )
                         && ( markedCells.find( m.first ) == markedCells.end() ) )
                      f_exterior.push_back( i );
                    else if ( markedFacet.find( facets[ i ] ) != markedFacet.end() )
                      f_interior.push_back( i );
                    if ( basicFacet.find( m ) != basicFacet.end() )
                      f_basic.push_back( i );
                  }
                if ( ( f_exterior.size() >= 3 ) ) //&& ( f_interior.size() == 2 ) )
                  {
                    //if ( f_basic.size() == 0 )
                    std::cerr << "Weird !"
                              << " ext=" << f_exterior.size()
                              << " int=" << f_interior.size()
                              << " bas=" << f_basic.size()
                              << std::endl;
                    // while ( f_basic.empty() )
                    //   {
                    //     f_exterior.erase
                    //   }
                  }
                if ( ( f_exterior.size() == 2 ) && ( f_interior.size() == 2 ) )
                  {
                    // (i,j) is the edge common to the two other faces.
                    int i = facets[ f_exterior[ 0 ] ].second; 
                    int j = facets[ f_exterior[ 1 ] ].second;
                    // (k,l) is the edge common to the two exterior faces.
                    int k = ( (i+1)%4 == j ) ? (i+2)%4 : (i+1)%4;
                    int l = (k+1)%4;
                    for ( ; (l == i ) || ( l == j ); l = (l+1)%4 ) ;
                    Point A( toDGtal( (*it )->vertex( i )->point() ) );
                    Point B( toDGtal( (*it )->vertex( j )->point() ) );
                    Point C( toDGtal( (*it )->vertex( k )->point() ) );
                    Point D( toDGtal( (*it )->vertex( l )->point() ) );

                    double angle = fabs( M_PI - computeDihedralAngle( t, facets[ f_exterior[ 0 ] ], 
                                                                      facets[ f_exterior[ 1 ] ] ) );
                    if ( angle < factor )
                      // if ( (C-D).norm( Point::L_1 ) > factor*(A-B).norm( Point::L_1 ) )
                      { // preference to shorter edge.
                        // std::cout << "(" << f_exterior.size() << "," << f_interior.size() << std::endl;
                        markedFacet.erase( t.mirror_facet( facets[ f_exterior[ 0 ] ] ) );
                        markedFacet.erase( t.mirror_facet( facets[ f_exterior[ 1 ] ] ) );
                        markedCells.erase( *it );
                        removedWeirdCells.insert( *it );
                        markedFacet.insert( facets[ f_interior[ 0 ] ] );
                        markedFacet.insert( facets[ f_interior[ 1 ] ] );
                        changed = true;
                      }
                  }
              }
          }
      } //  while ( changed ) {

    }
  trace.info() << "weird cells removed : " << removedWeirdCells.size() << std::endl;
  trace.endBlock();

  // Color colBasicFacet2( 0, 255, 255, 255 );
  // Color colBasicFacet1( 0, 255, 0, 255 );
  for ( FacetSet::const_iterator it = markedFacet.begin(), itend = markedFacet.end();
        it != itend; ++it )
    {
      Cell_handle cell = it->first;
      if ( markedCells.find( cell ) == markedCells.end() )
        { // we display it.
          Triangle triangle = t.triangle( *it );
          Point a( toDGtal( triangle.vertex( 0 ) ) );
          Point b( toDGtal( triangle.vertex( 1 ) ) );
          Point c( toDGtal( triangle.vertex( 2 ) ) );
          // int i = it->second;
          // Point a( toDGtal( cell->vertex( (i+1)%4 )->point() ) );
          // Point b( toDGtal( cell->vertex( (i+2)%4 )->point() ) );
          // Point c( toDGtal( cell->vertex( (i+3)%4 )->point() ) );
          Facet f2 = t.mirror_facet( *it );
          if ( ( markedFacet.find( f2 ) != markedFacet.end() )
               && ( markedCells.find( f2.first ) == markedCells.end() ) )
            { // the mirror facet is also in the triangulation. 
              // We need to move vertices a little bit when two triangles are at the same position.
              Point n = (b-a)^(c-a);
              double norm = n.norm(Point::L_2);
              double dx[ 3 ];
              for ( unsigned int j = 0; j < 3; ++j )
                dx[ j ] = 0.001*((double) n[j])/norm;
              RealPoint A( (double) a[ 0 ] + dx[ 0 ], (double) a[ 1 ] +  dx[ 1 ], (double) a[ 2 ] + dx[ 2 ] );
              RealPoint B( (double) b[ 0 ] + dx[ 0 ], (double) b[ 1 ] +  dx[ 1 ], (double) b[ 2 ] + dx[ 2 ] );
              RealPoint C( (double) c[ 0 ] + dx[ 0 ], (double) c[ 1 ] +  dx[ 1 ], (double) c[ 2 ] + dx[ 2 ] );
              viewer2.setFillColor( colBasicFacet2 );
              viewer2.addTriangle( A, C, B );
            }
          else
            {
              RealPoint A( a[ 0 ], a[ 1 ], a[ 2 ] );
              RealPoint B( b[ 0 ], b[ 1 ], b[ 2 ] );
              RealPoint C( c[ 0 ], c[ 1 ], c[ 2 ] );
              viewer2.setFillColor( colBasicFacet1 );
              viewer2.addTriangle( A, B, C );
            }
        }
    }

  //////////////////////////////////////////////////////////////////////
  for ( Cell_iterator it = t.cells_begin(), itend = t.cells_end();
        it != itend; ++it )
    {
      bool draw = false;
      Color col;
      if ( weirdCells.find( it ) != weirdCells.end() )
        {
          col = Color( 255, 100, 255, 100 );
          draw = true;
        }
      if ( prune && ( removedWeirdCells.find( it ) != removedWeirdCells.end() ) )
        {
          col = Color( 255, 255, 0, 100 );
          draw = true;
        }
      else
        if ( ( ! t.is_infinite( it ) )
             && ( nbLatticePoints[ it ] == 4 )
             && ( markedCells.find( it ) == markedCells.end() ) )
          {
            draw = true;
            col = Color( 255, 0, 0, 100 );
          }
      //if ( n > 4 ) trace.info() << " " << n;
      // if ( n <= 3 ) trace.error() << " Not enough lattice points" << std::endl;
      if ( draw )
        {
          RealPoint a( toDGtalR(it->vertex(0)->point()));
          RealPoint b( toDGtalR(it->vertex(1)->point()));
          RealPoint c( toDGtalR(it->vertex(2)->point()));
          RealPoint d( toDGtalR(it->vertex(3)->point()));
          viewer1.setFillColor( col );
          viewer1.addTriangle( a, b, c );
          viewer1.addTriangle( c, b, d );
          viewer1.addTriangle( c, d, a );
          viewer1.addTriangle( d, b, a );
        }
    }

  std::cout << "number of vertices :  " ;
  std::cout << t.number_of_vertices() << std::endl;
  std::cout << "number of edges :  " ;
  std::cout << t.number_of_edges() << std::endl;
  std::cout << "number of facets :  " ;
  std::cout << t.number_of_facets() << std::endl;
  std::cout << "number of cells :  " ;
  std::cout << t.number_of_cells() << std::endl;

  for ( SCellSetConstIterator it=boundary2.begin(), itend=boundary2.end(); it != itend; ++it )
    viewer3 << ks.sDirectIncident( *it, ks.sOrthDir( *it ) );

  //viewer1 << MViewer3D::updateDisplay;
  // viewer2 << MViewer3D::updateDisplay;
  viewer3 << MViewer3D::updateDisplay;
  application.exec();
  
  return 0;
}
