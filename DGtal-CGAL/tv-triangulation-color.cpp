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
#include <DGtal/images/ImageSelector.h>
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
///////////////////////////////////////////////////////////////////////////////
double randomUniform()
{
  return (double) random() / (double) RAND_MAX;
}

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
    typedef Triangulation::Vertex      Vertex;
    typedef Triangulation::Arc         Arc;
    typedef Triangulation::Face        Face;
    typedef Triangulation::VertexRange VertexRange;
    typedef double                     Scalar;
    typedef PointVector< 3, Scalar >   Value;
    
    /// The domain triangulation
    Triangulation        T;
    /// The values at each vertex
    std::vector< Value > values;
    /// The probability to flip an edge if energies are equal.
    Scalar              _equal_flip_prob;
    
    static Scalar square( Scalar x ) { return x*x; }
    
    /// The function that evaluates the energy at each triangle.
    Scalar energyTV( VertexIndex v1, VertexIndex v2, VertexIndex v3 ) const
    {
      const Value One = Value::diagonal( 1 );
      const Point& p1 = T.position( v1 );
      const Point& p2 = T.position( v2 );
      const Point& p3 = T.position( v3 );
      const Value& f1 = values[ v1 ];
      const Value& f2 = values[ v2 ];
      const Value& f3 = values[ v3 ];
      //trace.info() << p1 << p2 << p2 << f1 << f2 << f3 << std::endl;
      const Value   X = Value( p1[ 0 ], p2[ 0 ], p3[ 0 ] );
      const Value   Y = Value( p1[ 1 ], p2[ 1 ], p3[ 1 ] );
      const Value  Vr = Value( f1[ 0 ], f2[ 0 ], f3[ 0 ] );
      const Value  Vg = Value( f1[ 1 ], f2[ 1 ], f3[ 1 ] );
      const Value  Vb = Value( f1[ 2 ], f2[ 2 ], f3[ 2 ] );
      // TV on red
      const Scalar er = square( Vr.crossProduct( Y ).dot( One ) )
	+ square( X.crossProduct( Vr ).dot( One ) );
      const Scalar eg = square( Vg.crossProduct( Y ).dot( One ) )
	+ square( X.crossProduct( Vg ).dot( One ) );
      const Scalar eb = square( Vb.crossProduct( Y ).dot( One ) )
	+ square( X.crossProduct( Vb ).dot( One ) );
      return sqrt( er ) + sqrt( eg ) + sqrt( eb );
      // const Scalar er = fabs( Vr.crossProduct( Y ).dot( One ) )
      // 	+ fabs( X.crossProduct( Vr ).dot( One ) );
      // const Scalar eg = fabs( Vg.crossProduct( Y ).dot( One ) )
      // 	+ fabs( X.crossProduct( Vg ).dot( One ) );
      // const Scalar eb = fabs( Vb.crossProduct( Y ).dot( One ) )
      // 	+ fabs( X.crossProduct( Vb ).dot( One ) );
      // return sqrt( er*er+eg*eg+eb*eb );
    }
    
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
	  // T.addTriangle( v, v + taille[ 0 ], v + 1 );
	  // T.addTriangle( v + taille[ 0 ], v + taille[ 0 ] + 1, v + 1 );
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
      for ( unsigned int val : I ) {
	values[ v++ ] = Value( (Scalar) fred  ( val ),
			       (Scalar) fgreen( val ),
			       (Scalar) fblue ( val ) );
      }
      _equal_flip_prob = 0.5;
    }

    static Integer doesTurnLeft( const Point& p, const Point& q, const Point& r )
    {
      const Point pq = q - p;
      const Point qr = r - q;
      return pq[ 0 ] * qr[ 1 ] - pq[ 1 ] * qr[ 0 ];
    }
    static Integer doesTurnLeft( const Point& pq, const Point& qr )
    {
      return pq[ 0 ] * qr[ 1 ] - pq[ 1 ] * qr[ 0 ];
    }
    bool isConvex( const VertexRange& V ) const
    {
      Point P[] = { T.position( V[ 1 ] ) - T.position( V[ 0 ] ),
		    T.position( V[ 2 ] ) - T.position( V[ 1 ] ),
		    T.position( V[ 3 ] ) - T.position( V[ 2 ] ),
		    T.position( V[ 0 ] ) - T.position( V[ 3 ] ) };
      // trace.info() << P[ 0 ] << P[ 1 ] << P[ 2 ] << P[ 3 ] << std::endl;
      // trace.info() << " " << doesTurnLeft( P[ 0 ], P[ 1 ] )
      // 		   << " " << doesTurnLeft( P[ 1 ], P[ 2 ] )
      // 		   << " " << doesTurnLeft( P[ 2 ], P[ 3 ] )
      // 		   << " " << doesTurnLeft( P[ 3 ], P[ 0 ] ) << std::endl;
      bool cvx = ( doesTurnLeft( P[ 0 ], P[ 1 ] ) < 0 )
	&&       ( doesTurnLeft( P[ 1 ], P[ 2 ] ) < 0 )
	&&       ( doesTurnLeft( P[ 2 ], P[ 3 ] ) < 0 )
	&&       ( doesTurnLeft( P[ 3 ], P[ 0 ] ) < 0 );
      // if ( ! cvx )
      // 	trace.info() << "(" << V[ 0 ] << "," << V[ 1 ]
      // 		     << "," << V[ 2 ] << "," << V[ 3 ] << ")"
      // 		     << ( cvx ? "CVX" : "CCV" ) << std::endl;
      return cvx;
    }
    // NB: process only arcs (s,t) where ( t > s ).
    bool updateArc( const Arc a, Scalar& energy ) {
      energy        = 0;
      VertexRange P = T.verticesAroundArc( a );
      if ( P.size() != 4 )   return false;
      if ( P[ 0 ] < P[ 2 ] ) return false;
      if ( ! isConvex( P ) ) return false;
      // Computes energies
      const Scalar  E012 = energyTV( P[ 0 ], P[ 1 ], P[ 2 ] );
      const Scalar  E023 = energyTV( P[ 0 ], P[ 2 ], P[ 3 ] );
      const Scalar  E013 = energyTV( P[ 0 ], P[ 1 ], P[ 3 ] );
      const Scalar  E123 = energyTV( P[ 1 ], P[ 2 ], P[ 3 ] );
      const Scalar Ecurr = E012 + E023;
      const Scalar Eflip = E013 + E123;
      // trace.info() << "(" << P[ 0 ] << "," << P[ 1 ] << "," << P[ 2 ]
      // 		   << "," << P[ 3 ] << ") ";
      // trace.info() << "Ecurr=" << Ecurr << " Eflip=" << Eflip << std::endl;
      // @todo Does not take into account equality for now.a
      if ( ( Eflip < Ecurr )
	   || ( ( Eflip == Ecurr )
		&& ( randomUniform() < _equal_flip_prob ) ) )
	{
	  T.flip( a );
	  energy = Eflip;
	  return true;
	}
      else 
	{
	  energy = Ecurr;
	  return false;
	}
    }
    Integer onePass( Scalar& total_energy, Scalar equal_flip_prob )
    {
      _equal_flip_prob  = equal_flip_prob;
      Integer nbflipped = 0;
      total_energy      = 0;
      Scalar energy     = 0;
      for ( Arc a = 0; a < T.nbArcs(); ++a )
	{
	  if ( updateArc( a, energy ) )
	    nbflipped++;
	  total_energy += energy;
	}
      trace.info() << "Energy = " << total_energy
		   << " nbflipped=" << nbflipped
		   << "/" << T.nbEdges() << std::endl;
      return nbflipped;
    }
    
  };

  // Useful function for viewing triangulations.

  /**
     This class is intended for visualizing Affine Valued
     triangulation with CAIRO.
  */
  class CairoViewerTV
  {
  public:
    typedef TVTriangulation          TVT;
    typedef TVT::Value               Value;
    typedef TVT::Triangulation       Triangulation;
    typedef TVT::VertexIndex         VertexIndex;
    typedef TVT::Vertex              Vertex;
    typedef TVT::Arc                 Arc;
    typedef TVT::Face                Face;
    typedef TVT::VertexRange         VertexRange;
    typedef TVT::Scalar              Scalar;
    typedef TVT::Point               Point;
    typedef Z2i::RealPoint           RealPoint;

  private:
    const double _redf, _greenf, _bluef;
    int _x0, _y0;
    int _width, _height;
    double _xf, _yf;
    bool _gouraud;
    cairo_surface_t* _surface;
    cairo_t* _cr;

  public:

    // enum Mode { Gray, Red, Green, Blue };

    /**
       Constructor. 
    */
    CairoViewerTV( int x0, int y0, int width, int height, 
		   double xfactor = 1.0, double yfactor = 1.0,
		   bool gouraud = false )
      : _redf( 1.0/255.0f ), _greenf( 1.0/255.0f ), _bluef( 1.0/255.0f ),
	_x0( x0 ), _y0( y0 ), _width( width ), _height( height ),
	_xf( xfactor ), _yf( yfactor ), _gouraud( gouraud )
    {
      _surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32,
					     width, height );
      _cr = cairo_create ( _surface );
      // Fill the background with black
      cairo_set_source_rgba( _cr, 0.0, 0.0, 0.0, 1.0 );
      cairo_rectangle ( _cr, 0, 0, _width, _height );
      cairo_fill( _cr );
    }
    /// Destructor.
    ~CairoViewerTV()
    {
      cairo_destroy( _cr );
      cairo_surface_destroy( _surface );
    }
  
    void save( const char* file_name ) const
    {
      cairo_surface_write_to_png( _surface, file_name );
    }

    inline int i( double x ) const
    {
      return (int)round( (x+0.5) * _xf ) - _x0;
    }
    
    inline int j( double y ) const
    {
      return _height - (int)(round( (y+0.5) * _yf ) - _y0) - 1;
    }
    
    void viewGouraudTriangle( RealPoint a, RealPoint b, RealPoint c, 
			      Value val_a, Value val_b, Value val_c ) 
    {
      cairo_pattern_t * pattern = cairo_pattern_create_mesh();
      /* Add a Gouraud-shaded triangle */
      cairo_mesh_pattern_begin_patch (pattern);
      cairo_mesh_pattern_move_to (pattern, i( a[ 0 ] ), j( a[ 1 ] ) );
      cairo_mesh_pattern_line_to (pattern, i( b[ 0 ] ), j( b[ 1 ] ) );
      cairo_mesh_pattern_line_to (pattern, i( c[ 0 ] ), j( c[ 1 ] ) );
      cairo_mesh_pattern_set_corner_color_rgb (pattern, 0,
					       val_a[ 0 ] * _redf,
					       val_a[ 1 ] * _greenf,
					       val_a[ 2 ] * _bluef );
      cairo_mesh_pattern_set_corner_color_rgb (pattern, 1,
					       val_b[ 0 ] * _redf,
					       val_b[ 1 ] * _greenf,
					       val_b[ 2 ] * _bluef );
      cairo_mesh_pattern_set_corner_color_rgb (pattern, 2,
					       val_c[ 0 ] * _redf,
					       val_c[ 1 ] * _greenf,
					       val_c[ 2 ] * _bluef );
      cairo_mesh_pattern_end_patch (pattern);
      cairo_set_source( _cr, pattern );
      cairo_move_to( _cr, i( a[ 0 ] ), j( a[ 1 ] ) );
      cairo_line_to( _cr, i( b[ 0 ] ), j( b[ 1 ] ) );
      cairo_line_to( _cr, i( c[ 0 ] ), j( c[ 1 ] ) );
      cairo_close_path( _cr );
      cairo_fill( _cr );
      cairo_pattern_destroy( pattern );
    }

    void viewFlatTriangle( RealPoint a, RealPoint b, RealPoint c, 
			   Value val )
    {
      cairo_set_source_rgb( _cr,
			    val[ 0 ] * _redf,
			    val[ 1 ] * _greenf,
			    val[ 2 ] * _bluef );
      cairo_move_to( _cr, i( a[ 0 ] ), j( a[ 1 ] ) );
      cairo_line_to( _cr, i( b[ 0 ] ), j( b[ 1 ] ) );
      cairo_line_to( _cr, i( c[ 0 ] ), j( c[ 1 ] ) );
      cairo_close_path( _cr );
      cairo_fill( _cr );
    }

    void viewTVTGouraudTriangle( TVT & tvT, Face f )
    {
      VertexRange V = tvT.T.verticesAroundFace( f );
      Point a = tvT.T.position( V[ 0 ] );
      Point b = tvT.T.position( V[ 1 ] );
      Point c = tvT.T.position( V[ 2 ] );
      viewGouraudTriangle( RealPoint( a[ 0 ], a[ 1 ] ),
			   RealPoint( b[ 0 ], b[ 1 ] ),
			   RealPoint( c[ 0 ], c[ 1 ] ),
			   tvT.values[ V[ 0 ] ],
			   tvT.values[ V[ 1 ] ],
			   tvT.values[ V[ 2 ] ] );
    }
    void viewTVTFlatTriangle( TVT & tvT, Face f )
    {
      VertexRange V = tvT.T.verticesAroundFace( f );
      Point       a = tvT.T.position( V[ 0 ] );
      Point       b = tvT.T.position( V[ 1 ] );
      Point       c = tvT.T.position( V[ 2 ] );
      Value     val = tvT.values[ V[ 0 ] ]
	+ tvT.values[ V[ 1 ] ] + tvT.values[ V[ 2 ] ];
      val          /= 3.0;
      viewFlatTriangle( RealPoint( a[ 0 ], a[ 1 ] ),
			RealPoint( b[ 0 ], b[ 1 ] ),
			RealPoint( c[ 0 ], c[ 1 ] ), val );
    }

    /**
       Displays the AVT with flat or Gouraud shading.
    */
    void view( TVT & tvT )
    {
      cairo_set_operator( _cr,  CAIRO_OPERATOR_ADD );
      // _redf   = ( mode == Red )   || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
      // _greenf = ( mode == Green ) || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
      // _bluef  = ( mode == Blue )  || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
      cairo_set_line_width( _cr, 0.0 ); 
      cairo_set_line_cap( _cr, CAIRO_LINE_CAP_BUTT );
      cairo_set_line_join( _cr, CAIRO_LINE_JOIN_BEVEL );
      for ( Face f = 0; f < tvT.T.nbFaces(); ++f )
	{
	  if ( f % 1000 == 0 )
	    trace.info() << f << "/" << tvT.T.nbFaces() << std::endl;
	  if ( _gouraud ) viewTVTGouraudTriangle( tvT, f );
	  else            viewTVTFlatTriangle   ( tvT, f );
	}
      // cairo_operator_t op;
      // switch ( mode ) {
      // case Red: 
      // case Green:
      // case Blue: 
      //   op = CAIRO_OPERATOR_ADD; break;
      // case Gray:
      // default: op = CAIRO_OPERATOR_SOURCE; break;
      // };
      // op = CAIRO_OPERATOR_ADD;
      // cairo_set_operator( _cr, op );
    }

  };
  
  void viewTVTriangulationColor
  ( TVTriangulation& tvT, double b, double x0, double y0, double x1, double y1,
    bool gouraud, std::string fname )
  {
    CairoViewerTV cviewer
      ( (int) round( x0 ), (int) round( y0 ), 
	(int) round( (x1+1 - x0) * b ), (int) round( (y1+1 - y0) * b ), 
	b, b, gouraud );
    cviewer.view( tvT );
    cviewer.save( fname.c_str() );
  }
  
  void viewTVTriangulationColorAll
  ( TVTriangulation& tvT, double b, double x0, double y0, double x1, double y1,
    std::string fname )
  {
    viewTVTriangulationColor( tvT, b, x0, y0, x1, y1, false, fname + ".png" );
    viewTVTriangulationColor( tvT, b, x0, y0, x1, y1, true, fname + "-g.png" );
  }
  

} // namespace DGtal



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
  typedef ImageSelector < Z2i::Domain, unsigned int>::Type Image;
  ColorToRedFunctor c2r;
  ColorToGreenFunctor c2g;
  ColorToBlueFunctor c2b;
  GrayToGrayFunctor g2g;

  std::string imageFileName = vm[ "input" ].as<std::string>();
  Image image = GenericReader<Image>::import( imageFileName ); 

  TVTriangulation TVT( image, c2r, c2g, c2b );
  trace.info() << TVT.T << std::endl;
  trace.endBlock();

  trace.beginBlock("Flipping triangulation");
  int iter = 0;
  int nb   = 0;
  double equal_flip_prob = 0.0;
  do {
    if ( iter++ > 100 ) break;
    double energy = 0.0;
    nb = TVT.onePass( energy, equal_flip_prob );
    equal_flip_prob *= 0.9;
    TVT.T.isValid();
  } while ( nb > 0 );
  trace.endBlock();

  trace.beginBlock("Displaying triangulation");
  bool gouraud = vm.count( "gouraud" );
  double     b = vm[ "bitmap" ].as<double>();
  double    x0 = 0.0;
  double    y0 = 0.0;
  double    x1 = (double) image.domain().upperBound()[ 0 ];
  double    y1 = (double) image.domain().upperBound()[ 1 ];
  viewTVTriangulationColorAll( TVT, b, x0, y0, x1, y1, "avt-after" );
  trace.endBlock();

  return 0;
}
