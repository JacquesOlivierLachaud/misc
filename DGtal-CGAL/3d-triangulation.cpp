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
#include "Triangulation3DHelper.h"
#include "SimplicialStrip3D.h"
#include "RelativeConvexHull.h"

typedef DGtal::SpaceND<3, DGtal::int64_t> Z3;
typedef Z3::Point PointZ3;
typedef Z3::Point RealPoint;
typedef DGtal::HyperRectDomain<Z3> Domain;
typedef Domain::ConstIterator DomainConstIterator;

template <typename TKernel3>
struct toCGALFunctor {
  typedef typename TKernel3::Point_3 Point3;
  inline Point3 operator()( const PointZ3 & p ) const
  {
    return Point3( p[0], p[1], p[2] );
  }
};

template <typename TKernel3>
struct toDGtalFunctor {
  typedef typename TKernel3::Point_3 Point3;
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
		   const Facet & f, const DGtal::Color & col )
{
  typedef typename ToDGtal::Point3 Point;
  Point a( f.first->vertex( (f.second+1)%4 )->point() );
  Point b( f.first->vertex( (f.second+2)%4 )->point() );
  Point c( f.first->vertex( (f.second+3)%4 )->point() );
  viewer.setFillColor( col );
  viewer.addTriangle( RealPoint( a.x(), a.y(), a.z() ),
		      RealPoint( c.x(), c.y(), c.z() ),
		      RealPoint( b.x(), b.y(), b.z() ) );
}

template <typename Viewer, typename ToDGtal, typename Cell>
void displayCell( Viewer & viewer, const ToDGtal & toDGtal,
		  const Cell & c, const DGtal::Color & col )
{
  for ( int i = 0; i < 4; ++i )
    displayFacet( viewer, toDGtal, std::make_pair( c, i ), col );
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

  typedef KhalimskySpaceND<3,DGtal::int64_t> K3;
  typedef Z3::Vector Vector;
  typedef Z3::RealPoint RealPoint;
  typedef K3::SCell SCell;
  typedef K3::SCellSet SCellSet;
  typedef SCellSet::const_iterator SCellSetConstIterator;
  // typedef ImplicitRoundedHyperCube<Z3> Shape;
  typedef RealPoint::Coordinate Ring;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_3<K>            Triangulation;
  typedef CGAL::Delaunay_triangulation_3<K>   Delaunay;
  typedef RelativeConvexHull<Triangulation,K> RCH;
  typedef RCH::CellHandle                     CellHandle;
  typedef RCH::Facet                          Facet;
  typedef RCH::Triangle                       Triangle;
  typedef RCH::FiniteFacetsIterator           FiniteFacetsIterator;
  typedef toCGALFunctor<K>                    ToCGAL;
  typedef toDGtalFunctor<K>                   ToDGtal;

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

  trace.beginBlock("Constructing the set of points");
  std::vector<PointZ3> digital_inside_points;
  std::vector<PointZ3> digital_outside_points;

  getInnerAndOuterVoxelCoordinates( digital_inside_points, digital_outside_points,
                                    ks, boundary.begin(), boundary.end() );
  trace.beginBlock("Shuffle points.");
  random_shuffle( digital_inside_points.begin(), digital_inside_points.end() );
  random_shuffle( digital_outside_points.begin(), digital_outside_points.end() );
  trace.endBlock();
  trace.endBlock();

  trace.beginBlock("Creating the Delaunay complex.");
  ToCGAL toCGAL;
  ToDGtal toDGtal;
  RCH rch;
  double setsize = (double) digital_inside_points.size()-1;
  trace.info() << "Vertices to process: " << setsize << std::endl;
  rch.add( toCGAL, digital_inside_points.begin(), digital_inside_points.end(), 1, 0 );
  double setsize2 = (double) digital_outside_points.size()-1;
  trace.info() << "Vertices to process: " << setsize2 << std::endl;
  rch.add( toCGAL, digital_outside_points.begin(), digital_outside_points.end(), 0, 1 );
  rch.terminate();
  std::cout << "number of vertices :  " ;
  std::cout << rch.T().number_of_vertices() << std::endl;
  std::cout << "number of edges :  " ;
  std::cout << rch.T().number_of_edges() << std::endl;
  std::cout << "number of facets :  " ;
  std::cout << rch.T().number_of_facets() << std::endl;
  std::cout << "number of cells :  " ;
  std::cout << rch.T().number_of_cells() << std::endl;

  trace.endBlock();

  Color colBasicFacet3( 255, 255, 0, 255 );
  Color colBasicFacet2( 0, 255, 255, 255 );
  Color colBasicFacet1( 0, 255, 0, 255 );
  Color colSpuriousTetrahedra( 255, 0, 0, 100 );
  trace.beginBlock("Computing the relative convex hull.");
  unsigned int i = 0;
  do 
    {
      trace.info() << "--------- step " << i << " ----------" << std::endl;
      Viewer3D<> viewerRCH;
      viewerRCH.show();
      for ( FiniteFacetsIterator it = rch.T().finite_facets_begin(), 
      	      itend = rch.T().finite_facets_end(); it != itend; ++it )
      	{
      	  // we display it.
      	  Facet f = *it;
      	  if ( rch.isCellExtended( f.first ) )
      	    displayCell( viewerRCH, toDGtal, f.first, colSpuriousTetrahedra );
      	  if ( rch.isFacetOnBoundary( f, 0 ) ) 
      	    displayFacet( viewerRCH, toDGtal, f, colBasicFacet1 );
      	  else if ( rch.isFacetOnRelativeBoundary( f, 0 ) ) 
      	    displayFacet( viewerRCH, toDGtal, f, colBasicFacet2 );
      	  else if ( rch.isFacetSubdivided( f, 0 ) ) 
      	    displayFacet( viewerRCH, toDGtal, f, colBasicFacet3 );
      	  else {
      	    f = rch.T().mirror_facet( f );
      	    if ( rch.isFacetOnBoundary( f, 0 ) ) 
      	      displayFacet( viewerRCH, toDGtal, f, colBasicFacet1 );
      	    else if ( rch.isFacetOnRelativeBoundary( f, 0 ) ) 
      	      displayFacet( viewerRCH, toDGtal, f, colBasicFacet2 );
      	    else if ( rch.isFacetSubdivided( f, 0 ) ) 
      	      displayFacet( viewerRCH, toDGtal, f, colBasicFacet3 );
      	  }
      	}
      viewerRCH << Viewer3D<>::updateDisplay;
      application.exec();
      ++i;
    }
  while ( rch.relativeHull( 0 ) );
  trace.endBlock();

  // start viewer
  int view = vm["view"].as<int>();
  Viewer3D<> viewerRCH;
  viewerRCH.show();
  if ( view & 0x1 ) { // View digital core.
    for ( FiniteFacetsIterator it = rch.T().finite_facets_begin(), 
            itend = rch.T().finite_facets_end(); it != itend; ++it )
      {
        // we display it.
        Facet f = *it;
	  if ( rch.isCellExtended( f.first ) )
	    displayCell( viewerRCH, toDGtal, f.first, colSpuriousTetrahedra );
	  if ( rch.isFacetOnBoundary( f, 0 ) ) 
	    displayFacet( viewerRCH, toDGtal, f, colBasicFacet1 );
	  else if ( rch.isFacetOnRelativeBoundary( f, 0 ) ) 
	    displayFacet( viewerRCH, toDGtal, f, colBasicFacet2 );
	  else if ( rch.isFacetSubdivided( f, 0 ) ) 
	    displayFacet( viewerRCH, toDGtal, f, colBasicFacet3 );
	  else {
	    f = rch.T().mirror_facet( f );
	    if ( rch.isFacetOnBoundary( f, 0 ) ) 
	      displayFacet( viewerRCH, toDGtal, f, colBasicFacet1 );
	    else if ( rch.isFacetOnRelativeBoundary( f, 0 ) ) 
	      displayFacet( viewerRCH, toDGtal, f, colBasicFacet2 );
	    else if ( rch.isFacetSubdivided( f, 0 ) ) 
	      displayFacet( viewerRCH, toDGtal, f, colBasicFacet3 );
	  }
      }
  }
  viewerRCH << Viewer3D<>::updateDisplay;
  application.exec();

  std::cout << "number of vertices :  " ;
  std::cout << rch.T().number_of_vertices() << std::endl;
  std::cout << "number of edges :  " ;
  std::cout << rch.T().number_of_edges() << std::endl;
  std::cout << "number of facets :  " ;
  std::cout << rch.T().number_of_facets() << std::endl;
  std::cout << "number of cells :  " ;
  std::cout << rch.T().number_of_cells() << std::endl;

  // for ( SCellSetConstIterator it=boundary.begin(), itend=boundary.end(); it != itend; ++it )
  //   viewer3 << ks.sDirectIncident( *it, ks.sOrthDir( *it ) );

  // viewer1 << Viewer3D::updateDisplay;
  // viewer2 << Viewer3D::updateDisplay;
  // viewer3 << Viewer3D::updateDisplay;
  // application.exec();
  
  return 0;
}
