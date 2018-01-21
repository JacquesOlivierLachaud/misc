#include <cfloat>
#include <iostream>
#include <vector>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <DGtal/base/Common.h>
#include <DGtal/base/ConstAlias.h>
#include <DGtal/base/Alias.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/shapes/TriangulatedSurface.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/GradientColorMap.h>
#include <DGtal/io/colormaps/GrayscaleColorMap.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/PPMWriter.h"

// #include <CGAL/Delaunay_triangulation_2.h>
// #include <CGAL/Constrained_Delaunay_triangulation_2.h>
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Triangulation_2.h>

#include "cairo.h"

// #include "Triangulation2DHelper.h"
// #include "UmbrellaPart2D.h"
// #include "Auxiliary.h"

namespace DGtal {
  struct ColorToRedFunctor {
    int operator()( int color ) const
    { return (color >> 16) & 0xff; }
  };
  struct ColorToGreenFunctor {
    int operator()( int color ) const
  { return (color >> 8) & 0xff; }
  };
  struct ColorToBlueFunctor {
    int operator()( int color ) const
    { return color & 0xff; }
  };
  struct GrayToGrayFunctor {
    int operator()( int color ) const
    { return color & 0xff; }
  };
  struct GrayToRedGreen {
    DGtal::Color operator()( int value ) const
    { 
      if ( value >= 128 ) return DGtal::Color( 0, (value-128)*2+1, 0 );
      else                return DGtal::Color( 0, 0, 255-value*2 );
    }
  };
  
  struct TVTriangulation
  {
    typedef Z2i::Integer               Integer;
    typedef Z2i::Point                 Point;
    typedef Z2i::Domain                Domain;
    typedef TriangulatedSurface<Point> Triangulation;
    typedef Triangulation::VertexIndex VertexIndex;
    typedef double                     Scalar;
    typedef PointVector< 3, Scalar >   Value;
    
    /// The domain triangulation
    Triangulation        T;
    /// The values at each vertex
    std::vector< Value > values;
    
    // Constructor from color image.
    template <typename Image, typename Fred, typename Fgreen, typename Fblue>
    TVTriangulation( const Image&  I,
		     const Fred&   fred,
		     const Fgreen& fgreen,
		     const Fblue&  fblue )
    {
      const Point taille = I.extent();
      // Creates vertices
      for ( auto p : I.domain() ) T.addVertex( p );
      // Creates triangles
      for ( Integer y = 0; y < taille[ 1 ] - 1; ++y ) {
	for ( Integer x = 0; x < taille[ 0 ] - 1; ++x ) {
	  VertexIndex v = y * taille[ 0 ] + x;
	  T.addTriangle( v, v + taille[ 0 ], v + taille[ 0 ] + 1 );
	  T.addTriangle( v, v + taille[ 0 ] + 1, v + 1 );
	}
      }
      bool ok = T.build();
      trace.info() << "Build triangulation: "
		   << ( ok ? "OK" : "ERROR" ) << std::endl;
      // Creates values
      VertexIndex v = 0;
      values.resize( I.size() );
      for ( auto val : I ) 
	values[ v ] = Value( (Scalar) fred  ( val ),
			     (Scalar) fgreen( val ),
			     (Scalar) fblue ( val ) );
    }
  };
}

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  using namespace DGtal;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "Specifies the input shape as a 2D image PPM filename.")
    ("random,r", po::value<double>()->default_value(1.0), "Keep only a proportion [arg] (between 0 and 1) of the input data point.")  
    ("limit,L", po::value<int>()->default_value(1000), "Gives the maximum number of passes (default is 10000).")  
    ("bitmap,b", po::value<double>()->default_value( 2.0 ), "Rasterization magnification factor [arg] for PNG export." )
    ("gouraud,g", "Displays faces with Gouraud-like shading.")  
    ("2ndorder,2", "Displays faces with 2nd order Gouraud-like shading.")  
    ;

  bool parseOK = true;
  po::variables_map vm;
  try {
    po::store( po::parse_command_line(argc, argv, general_opt), vm );  
  } catch ( const std::exception& ex ) {
    parseOK = false;
    trace.info() << "Error checking program options: " << ex.what() << std::endl;
  }
    
  po::notify(vm);    
  if( ! parseOK || vm.count("help") || argc <= 1 || (! vm.count( "input") ) )
    {
      trace.info()<< "Generate the best piecewise linear TV from a color image. It works by simply flipping edges of some initial triangulation to optimize the given energy." <<std::endl << "Basic usage: " << std::endl
		  << "\ttv-triangulation-color [options] -i <image.ppm> -b 4"<<std::endl
		  << general_opt << "\n";
      return 0;
    }

  // Useful types
  typedef DGtal::Z2i::Domain Domain;

  using namespace DGtal;

  trace.beginBlock("Construction of the triangulation");

  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned int> Image;
  ColorToRedFunctor c2r;
  ColorToGreenFunctor c2g;
  ColorToBlueFunctor c2b;
  GrayToGrayFunctor g2g;

  std::string imageFileName = vm[ "input" ].as<std::string>();
  Image image = GenericReader<Image>::import( imageFileName ); 

  TVTriangulation TVT( image, c2r, c2g, c2b );
  trace.info() << TVT.T << std::endl;
  trace.endBlock();

  bool gouraud = vm.count( "gouraud" );
  double b = vm[ "bitmap" ].as<double>();

  return 0;
}
