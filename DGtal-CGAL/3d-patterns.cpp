#include <iostream>
#include <vector>
#include <string>

#include <QtGui/qapplication.h>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/Shapes.h>
#include <DGtal/shapes/ShapeFactory.h>
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/MPolynomialReader.h"

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CORE/Expr.h>

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
typedef DGtal::SpaceND<3, DGtal::int64_t> Z3;
typedef Z3::Point Point;
typedef Z3::RealPoint RealPoint;
typedef DGtal::HyperRectDomain<Z3> Domain;
typedef Domain::ConstIterator DomainConstIterator;

struct EdgeIteratorComparator
{
  inline
  bool operator()( const Edge_iterator & e1, const Edge_iterator & e2 ) const
  {
    return ( e1->first < e2->first )
      || ( ( e1->first == e2->first ) 
           && ( ( e1->second < e2->second )
                || ( ( e1->second == e2->second )
                     && ( e1->third < e2->third ) ) ) );
  }
};

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
    if ( second < first  ) std::swap( first, second );
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

typedef std::map< Cell_iterator, DGtal::int64_t > MapCell2Int;
typedef std::map< VFacet, DGtal::int64_t > MapFacet2Int;
typedef std::map< VEdge, DGtal::int64_t > MapEdge2Int;
typedef std::map< Vertex_iterator, DGtal::int64_t > MapVertex2Int;


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

DGtal::uint64_t 
markBasicFacets( MapFacet2Int & basicFacetMap, const Delaunay & t, MapEdge2Int & basicEdgeMap )
{
  DGtal::uint64_t nb = 0;
  for( Facet_iterator it = t.facets_begin(), itend = t.facets_end();
       it != itend; ++it)
    {
      VFacet f( *it ); 
      if ( ! t.is_infinite( *it ) )
        {
          Cell_iterator itCell = it->first; int i = it->second;
          VEdge e1( f.first, f.second );
          VEdge e2( f.second, f.third );
          VEdge e3( f.third, f.first );
          unsigned int n = 0;
          n += basicEdgeMap[ e1 ] == 1 ? 1 : 0;
          n += basicEdgeMap[ e2 ] == 1 ? 1 : 0;
          n += basicEdgeMap[ e3 ] == 1 ? 1 : 0;
          //ASSERT( n < 3 );
          if ( n == 3 )
            {
              basicFacetMap[ f ] = 1;
              ++nb;
            }
          else
            basicFacetMap[ f ] = 0;
        }
      else
        basicFacetMap[ f ] = -1;
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

void 
getFacets( std::vector<Facet> & facets, const Cell_iterator & it )
{
  for ( int i = 0; i < 4; ++i )
    facets.push_back( Facet( it, i ) );
}

/*
  "29*z+47*y+23*x-5"
  "10*z-x^2+y^2-100"
  "z*x*y+x^4-5*x^2+2*y^2*z-z^2-1000"
  "(15*z-x^2+y^2-100)*(x^2+y^2+z^2-1000)" nice
 */
int main (int argc, char** argv )
{
  using namespace DGtal;

  
  typedef KhalimskySpaceND<3,DGtal::int64_t> K3;
  typedef Z3::Vector Vector;
  typedef Z3::RealPoint RealPoint;
  typedef K3::SCell SCell;
  typedef K3::SCellSet SCellSet;
  typedef SCellSet::const_iterator SCellSetConstIterator;
  // typedef ImplicitRoundedHyperCube<Z3> Shape;
  typedef RealPoint::Coordinate Ring;
  typedef MPolynomial<3, Ring> Polynomial3;
  typedef MPolynomialReader<3, Ring> Polynomial3Reader;
  
  typedef ImplicitPolynomial3Shape<Z3> Shape;
  typedef GaussDigitizer<Z3,Shape> Digitizer;


  QApplication application(argc,argv); // remove Qt arguments.

  Delaunay t;
  
  trace.beginBlock("Construction of the shape");
  //Shape shape( RealPoint( 0.0, 0.0, 0.0 ), 14, 3.0 );
  Polynomial3 P;
  Polynomial3Reader reader;
  std::string poly_str = argc > 1 ? argv[ 1 ] : "x^4+y^4+z^4-50000";
  std::string::const_iterator iter 
    = reader.read( P, poly_str.begin(), poly_str.end() );
  if ( iter != poly_str.end() )
    {
      std::cerr << "ERROR: I read only <" 
                << poly_str.substr( 0, iter - poly_str.begin() )
                << ">, and I built P=" << P << std::endl;
      return 1;
    }
  trace.info() << "P( X_0, X_1, X_2 ) = " << P << std::endl;

  Shape shape( P );
  double h = 1.0; 
  Digitizer dig;
  dig.attach( shape );
  // dig.init( shape.getLowerBound()+Vector::diagonal(-1),
  //           shape.getUpperBound()+Vector::diagonal(1), h ); 
  dig.init( Vector::diagonal(-50),
            Vector::diagonal(50), h ); 
  K3 ks;
  ks.init( dig.getLowerBound(), dig.getUpperBound(), true );
  SurfelAdjacency<3> sAdj( true );
  SCell bel = Surfaces<K3>::findABel( ks, dig, 10000 );
  SCellSet boundary;
  Surfaces<K3>::trackBoundary( boundary, ks, sAdj, dig, bel );
  trace.endBlock();

  trace.beginBlock("Delaunay");
  std::set<SCell> inner_points;
  for( SCellSetConstIterator it = boundary.begin(), itE = boundary.end(); it != itE; ++it )
    // Get inner point.
    inner_points.insert( ks.sDirectIncident( *it, ks.sOrthDir( *it ) ) );
  for ( std::set<SCell>::const_iterator it = inner_points.begin(), itE = inner_points.end();
	it != itE; ++it )
    t.insert( toCGAL( ks.sCoords( *it ) ) );
  trace.endBlock();

  // start viewer
  Viewer3D<> viewer;
  viewer.show();
  // for ( SCellSetConstIterator it=boundary.begin(), itend=boundary.end(); it != itend; ++it )
  //    viewer << *it;
  // for ( SCellSetConstIterator it=boundary.begin(), itend=boundary.end(); it != itend; ++it )
  //   viewer << ks.sDirectIncident( *it, ks.sOrthDir( *it ) );
  
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
          nbLatticePoints[ it ] = countLatticePointsInTetrahedra( a, b, c, d );
        }
    }

  MapEdge2Int bEdge;
  uint64_t nbBE = markBasicEdges( bEdge, t );
  trace.info() << "Nb basic edges :" << nbBE << std::endl;
  MapFacet2Int bFacet;
  uint64_t nbBF = markBasicFacets( bFacet, t, bEdge );
  trace.info() << "Nb basic facets:" << nbBF << std::endl;

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
          RealPoint a( toDGtalR( itC->vertex( i )->point() ) );
          RealPoint b( toDGtalR( itC->vertex( j )->point() ) );
          viewer.setLineColor( colBasicEdge );
	  viewer.addLine( a, b, 1.0 );
        }
    }

  // Propagate basic facets.
  std::queue<Cell_iterator> Q;
  for ( Cell_iterator it = t.cells_begin(), itend = t.cells_end();
        it != itend; ++it )
    if ( nbLatticePoints[ it ] == 4 ) Q.push( it );
  while ( ! Q.empty() )
    {
      Cell_iterator it = Q.front(); Q.pop();
      std::vector<Facet> facets;
      getFacets( facets, it );
      int j = -2;
      unsigned int n = 0;
      for ( int i = 0; i < facets.size(); ++i )
        {
          VFacet f( facets[ i ] );
          if ( bFacet[ f ] > 0 ) ++n;
        }
      std::cout << " " << n;
      if ( ( n >= 2 ) && ( n <= 3 ) ) // 2 or 3 facets were "basic"
        {
          for ( int i = 0; i < facets.size(); ++i )
            {
              std::cout << "+"; 
              VFacet f( facets[ i ] );
              if ( bFacet[ f ] <= 0 ) // contaminate exterior facets
                {
                  bFacet[ f ] = 1;
                  Facet f2 = t.mirror_facet( facets[ i ] );
                  if ( nbLatticePoints[ f2.first ] == 4 )
                    Q.push( f2.first );
                }
              else 
                {
                  bFacet[ f ] = 2; // hide inner facets
                }
            }
        }
    }

  Color colBasicFacet2( 0, 255, 255, 255 );
  Color colBasicFacet1( 0, 255, 0, 255 );
  for( Facet_iterator it = t.facets_begin(), itend = t.facets_end();
       it != itend; ++it)
    {
      VFacet f( *it );
      //std::cout << " " << bFacet[ f ];
      if ( bFacet[ f ] > 0 )
        {
          RealPoint a( toDGtalR( f.first->point() ) );
          RealPoint b( toDGtalR( f.second->point() ) );
          RealPoint c( toDGtalR( f.third->point() ) );
          viewer.setFillColor( (bFacet[ f ] == 1) ? colBasicFacet1 : colBasicFacet2 ); 
	  viewer.addTriangle( a, b, c );
        }
    }
  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////
  for ( Cell_iterator it = t.cells_begin(), itend = t.cells_end();
        it != itend; ++it )
    {
      if ( ( ! t.is_infinite( it ) )
           && ( nbLatticePoints[ it ] == 4 ) )
        //if ( n > 4 ) trace.info() << " " << n;
      // if ( n <= 3 ) trace.error() << " Not enough lattice points" << std::endl;
        {
          RealPoint a( toDGtalR( it->vertex(0)->point() ) );
          RealPoint b( toDGtalR( it->vertex(1)->point() ) );
          RealPoint c( toDGtalR( it->vertex(2)->point() ) );
          RealPoint d( toDGtalR( it->vertex(3)->point() ) );
          viewer.setFillColor( Color( 255, 0, 0, 100 ) );
          viewer.addTriangle( a, b, c );
          viewer.addTriangle( c, b, d );
          viewer.addTriangle( c, d, a );
          viewer.addTriangle( d, b, a );
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

  viewer << Viewer3D<>::updateDisplay;
  application.exec();
  
  return 0;
}
