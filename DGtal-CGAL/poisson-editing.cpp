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
#include <DGtal/kernel/BasicPointPredicates.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/PPMWriter.h>
#include <DGtal/io/writers/PGMWriter.h>
#include <DGtal/math/RealFFT.h>
#include "FourierPoisson.h"

// struct UnsignedInt2Color {
//   Color operator()( unsigned int val ) const { return Color( val ); }
// };
struct UnsignedInt2GrayLevel {
  unsigned char operator()( unsigned int val ) const {
    return static_cast<unsigned char>( val );
  }
};


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;
///////////////////////////////////////////////////////////////////////////////

double sqr( double x ) { return x*x; }

int main( int argc, char** argv )
{
  using namespace DGtal;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "Specifies the input shape as a 2D image PPM filename.")
    ("output,o", po::value<std::string>()->default_value("after-poisson.pgm"), "Specifies the basename of output displayed images.");
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
      trace.info()<< "Poisson editing." <<std::endl << "Basic usage: " << std::endl
		  << "\tpoisson-editing [options] -i <image.ppm> "<<std::endl
		  << general_opt << "\n";
      return 0;
    }

  // Useful types
  using namespace std;
  using namespace DGtal;

  typedef Z2i::Space      Space;
  typedef Z2i::Domain     Domain;
  typedef Z2i::Point      Point;
  typedef Z2i::RealVector RealVector;
  typedef ImageSelector < Z2i::Domain, unsigned int>::Type Image;
  typedef ImageSelector < Z2i::Domain, Color>::Type ColorImage;
  typedef ImageSelector < Z2i::Domain, RealVector>::Type RealVectorImage;
  
  trace.beginBlock("Loading image");
  std::string img_fname = vm[ "input" ].as<std::string>();
  std::string out_fname = vm[ "output" ].as<std::string>();
  Image image           = GenericReader<Image>::import( img_fname ); 
  std::string extension = img_fname.substr(img_fname.find_last_of(".") + 1);
  bool            color = false;
  if ( extension == "ppm" ) color = true;
  trace.info() << "Image <" << img_fname
	       << "> size=" << image.extent()[ 0 ]
	       << "x" << image.extent()[ 1 ]
	       << " color=" << ( color ? "True" : "False" ) << std::endl;
  //trace.info() << std::fixed;
  trace.endBlock();

  //DGtal::image::functions::bandFilter( image, 0.1, 1.0, 128.0 );
  PGMWriter<Image,UnsignedInt2GrayLevel>::exportPGM( "filtered.pgm", image );

  const Dimension d = color ? 3 : 1;
  Domain domain( image.domain() );
  Domain sym_domain( Domain( -domain.upperBound() - Point::diagonal( 1 ),
			     domain.upperBound() ) );
  RealFFT< Domain, double > U ( sym_domain );
  // Calcule l'image gradient.
  RealFFT< Domain, double > Vx( sym_domain );
  RealFFT< Domain, double > Vy( sym_domain );
  auto IU  = U.getSpatialImage();
  auto IVx = Vx.getSpatialImage();
  auto IVy = Vy.getSpatialImage();
  for ( unsigned int y = 0; y < image.extent()[ 1 ]; ++y )
     for ( unsigned int x = 0; x < image.extent()[ 0 ]; ++x )
       {
	 // Extend by symetry
	 Point p( x, y );
	 Point sym_px ( -1 - p[ 0 ], p[ 1 ] );
	 Point sym_py ( p[ 0 ], -1 - p[ 1 ] );
	 Point sym_pxy( -1 - p[ 0 ], -1 - p[ 1 ] );
	 // Point px( p[0] > 0 ? p[0]-1 : p[ 0 ], p[1] ); 
	 // Point py( p[0], p[1] > 0 ? p[1]-1 : p[ 1 ] );
	 // double gx = ( (double)image( p ) - (double)image( px ) ) / 1.0;
	 // double gy = ( (double)image( p ) - (double)image( py ) ) / 1.0;
	 std::array<Point,4> all_p = { p, sym_px, sym_py, sym_pxy }; 
	 for ( auto&&q : all_p )
	   IU.setValue(  q, (double) image( p ) );
       }
  for ( auto&& p : sym_domain )
    {
      Point px( p[0]+1, p[1] ); px = px.inf( sym_domain.upperBound() );
      Point py( p[0], p[1]+1 ); py = py.inf( sym_domain.upperBound() );
      Point bpx( p[0]-1, p[1] ); bpx = bpx.sup( sym_domain.lowerBound() );
      Point bpy( p[0], p[1]-1 ); bpy = bpy.sup( sym_domain.lowerBound() );
      //double gx = ( IU( px ) - IU( p ) ) / 1.0;
      //double gy = ( IU( py ) - IU( p ) ) / 1.0;
      double gx = ( IU( px ) - IU( bpx ) ) / 2.0;
      double gy = ( IU( py ) - IU( bpy ) ) / 2.0;
      if ( ( gx*gx + gy*gy ) >= 500.0 ) {
  	IVx.setValue( p, gx );
  	IVy.setValue( p, gy );
      } else {
  	IVx.setValue( p, 0.0 );
  	IVy.setValue( p, 0.0 );
      }
    }
  trace.info() << "Compute FFT[V]" << std::endl;
  trace.info() << "spatial   domain = " << Vx.getSpatialDomain() << std::endl;
  trace.info() << "frequency domain = " << Vx.getFreqDomain() << std::endl;
  U.forwardFFT();
  Vx.forwardFFT();
  Vy.forwardFFT();
  Domain freq_domain = Vx.getFreqDomain();
  auto FU  = U.getFreqImage();
  auto FVx = Vx.getFreqImage();
  auto FVy = Vy.getFreqImage();
  typedef RealFFT< Domain, double >::Complex Complex;
  const double two_pi = 2.0 * M_PI;
  const Complex     i = Complex( 0.0, 1.0 );
  const double      J = ( sym_domain.upperBound()[0] - sym_domain.lowerBound()[0] ) + 1;
  const double      L = ( sym_domain.upperBound()[1] - sym_domain.lowerBound()[1] ) + 1;
  for ( auto&& p : freq_domain )
    {
      auto  freq = U.calcScaledFreqCoords( p );
      const double  mj = freq[ 0 ];
      // const double smj = freq[ 0 ] >= 0.0 ? freq[ 0 ] : freq[ 0 ] + 1.0 ;
      // double  mj = ( p[ 0 ] >= (J/2) ) ? ( ( p[ 0 ] - J ) / J ) : ( p[ 0 ] / J );
      const double  nl = freq[ 1 ];
      // const double snl = freq[ 1 ] >= 0.0 ? freq[ 1 ] : freq[ 1 ] + 1.0 ;
      // double  nl = ( p[ 1 ] >= (L/2) ) ? ( ( p[ 1 ] - L ) / L ) : ( p[ 1 ] / L );
      Complex val = two_pi * i * ( mj * FVx( p ) + nl * FVy( p ) )
	/ ( sqr( two_pi * mj ) + sqr( two_pi * nl ) );
      // std::cout << "(" << mj << "," << nl << ")" << std::endl;
      //if ( mj == 0.0 && nl == 0.0 )
      if ( p == Point::zero ) {
	std::cout << "[0,0] = " << FU( p )/(256.*256.) << std::endl;
	// FU.setValue( p, 0.0 );
      }
      else if ( mj >= 0.5 )
	FU.setValue( p, -val );
      else
	FU.setValue( p, -val );
      //if ( p[ 0 ] == J/2 ) std::cout << freq << " " << mj << " " << nl << std::endl;
    }
  U.backwardFFT();
  for ( unsigned int y = 0; y < image.extent()[ 1 ]; ++y )
    for ( unsigned int x = 0; x < image.extent()[ 0 ]; ++x )
      {
	Point p( x, y );
	//std:: cout << IU( p ) << std::endl;
	image.setValue( p,
			(unsigned char)
			std::max( 0.0, std::min( 255.0, IU( p ) ) ) );
      }
  PGMWriter<Image,UnsignedInt2GrayLevel>::exportPGM( out_fname.c_str(), image );
  
  return 0;
}
  
  
