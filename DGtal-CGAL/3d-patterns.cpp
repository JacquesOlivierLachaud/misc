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
typedef Delaunay::Cell_iterator  Cell_iterator;
typedef Delaunay::Point             CGALPoint;
typedef Delaunay::Cell_handle  Cell_handle;

typedef DGtal::SpaceND<3, DGtal::int64_t> Z3;
typedef Z3::Point Point;
typedef DGtal::HyperRectDomain<Z3> Domain;
typedef Domain::ConstIterator DomainConstIterator;

Point toDGtal(const CGALPoint &p)
{
  return Point ( p.x(),
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


using namespace DGtal;

int main (int argc, char** argv )
{
  
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
  Viewer3D viewer;
  viewer.show();
  // for ( SCellSetConstIterator it=boundary.begin(), itend=boundary.end(); it != itend; ++it )
  //   viewer << *it;

  
  for(Cell_iterator it = t.cells_begin(), itend=t.cells_end();
      it != itend; ++it)
    {
      if ( t.is_infinite( it ) ) continue;
      // bool inf_v = false;
      // for ( unsigned int j = 0; j < 4; ++j )
      // 	if ( it->vertex( j ) == t.infinite_vertex() )
      // 	  inf_v = true;
      // if ( inf_v ) continue;

      Point a( toDGtal(it->vertex(0)->point())),
	b(toDGtal(it->vertex(1)->point())),
	c(toDGtal(it->vertex(2)->point())),
	d(toDGtal(it->vertex(3)->point()));

      // Vector ab( b - a ), ac( c - a );
      // int d = ab[ 0 ] * ac[ 1 ] - ab[ 1 ] * ac[ 0 ];
      DGtal::int64_t n = countLatticePointsInTetrahedra( a, b, c, d );
      if ( n > 4 ) trace.info() << " " << n;
      // if ( n <= 3 ) trace.error() << " Not enough lattice points" << std::endl;
      if ( n == 4 )
	{
	  Color col( 255, 0, 0, 255 );
	  viewer.addTriangle( a[ 0 ], a[ 1 ], a[ 2 ],
			      b[ 0 ], b[ 1 ], b[ 2 ],
			      c[ 0 ], c[ 1 ], c[ 2 ],
			      col );
	  viewer.addTriangle( c[ 0 ], c[ 1 ], c[ 2 ],
			      b[ 0 ], b[ 1 ], b[ 2 ],
			      d[ 0 ], d[ 1 ], d[ 2 ], 
			      col );
	  viewer.addTriangle( c[ 0 ], c[ 1 ], c[ 2 ],
			      d[ 0 ], d[ 1 ], d[ 2 ], 
			      a[ 0 ], a[ 1 ], a[ 2 ],
			      col );
	  viewer.addTriangle( d[ 0 ], d[ 1 ], d[ 2 ], 
			      b[ 0 ], b[ 1 ], b[ 2 ],
			      a[ 0 ], a[ 1 ], a[ 2 ],
			      col );
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

  viewer << Viewer3D::updateDisplay;
  application.exec();
  
  return 0;
}
