#include <cfloat>
#include <iostream>
#include <vector>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageSelector.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/PPMWriter.h"
#include "DGtal/io/writers/PGMWriter.h"
#include "ImageTVRegularization.h"

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
    ("input,i", po::value<std::string>(), "Specifies the input 2D image filename.")
    ("output,o", po::value<std::string>()->default_value( "tv-output.ppm" ), "Specifies the name of the output 2D image.")
    ("lambda,l", po::value<double>(), "The data fidelity term in TV denoising (10: very high, 0.01: very low" ) 
    ("dt", po::value<double>()->default_value( 0.248 ), "The time step in TV denoising (should be lower than 0.25)" ) 
    ("tolerance,t", po::value<double>()->default_value( 0.01 ), "The tolerance to stop the TV denoising." ) 
    ("tv-max-iter,N", po::value<int>()->default_value( 10 ), "The maximum number of iteration in TV's algorithm." )
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
      trace.info()<< "Computes the Total Variation regularization of a color ppm or gray-level pgm image. It is an implementation of Chambolle, Pock primal-dual algorithm." <<std::endl << "Basic usage: " << std::endl
		  << "\ttv-image [options] -i <image.ppm> -l 0.1 "<<std::endl
		  << general_opt << "\n";
      return 0;
    }

  // Useful types

  using namespace DGtal;

  trace.beginBlock("Reading input image");
  typedef Z2i::Space         Space;
  typedef Z2i::Domain        Domain;
  typedef ImageSelector <Domain, unsigned int>::Type Image;
  std::string img_fname = vm[ "input" ].as<std::string>();
  Image           image = GenericReader<Image>::import( img_fname ); 
  std::string input_ext = img_fname.substr(img_fname.find_last_of(".") + 1);
  bool            color = false;
  if ( input_ext == "ppm" ) color = true;
  trace.info() << "Image <" << img_fname
	       << "> size=" << image.extent()[ 0 ]
	       << "x" << image.extent()[ 1 ]
	       << " color=" << ( color ? "True" : "False" ) << std::endl;
  trace.endBlock();

  trace.info() << std::fixed;
  Image output_u( image.domain() );
  std::string out_fname = vm[ "output" ].as<std::string>();
  std::string   out_ext = out_fname.substr(out_fname.find_last_of(".") + 1);
  bool        out_color = false;
  if ( out_ext == "ppm" ) out_color = true;

  trace.beginBlock("TV regularization");
  double lambda = vm[ "lambda" ].as<double>();
  double     dt = vm[ "dt" ].as<double>();
  double    tol = vm[ "tolerance" ].as<double>();
  int  max_iter = vm[ "tv-max-iter" ].as<int>();
  if ( color ) {
    typedef ImageTVRegularization<Space, 3> ColorTV;
    ColorTV tv;
    tv.init( image, ColorTV::Color2ValueFunctor() );
    tv.optimize( lambda, dt, tol, max_iter );
    if ( out_color ) tv.outputU( output_u, ColorTV::Value2ColorFunctor() );
    else             tv.outputU( output_u, ColorTV::Value2GrayLevelFunctor() );
  } else {
    typedef ImageTVRegularization<Space, 1> GrayLevelTV;
    GrayLevelTV tv;
    tv.init( image, GrayLevelTV::GrayLevel2ValueFunctor() );
    tv.optimize( lambda, dt, tol, max_iter );
    if ( out_color ) tv.outputU( output_u, GrayLevelTV::Value2ColorFunctor() );
    else             tv.outputU( output_u, GrayLevelTV::Value2GrayLevelFunctor() );
  }
  trace.endBlock();
  
  trace.beginBlock("Output TV image (possibly quantified)");
  struct UnsignedInt2Color {
    Color operator()( unsigned int val ) const { return Color( val ); }
  };
  struct UnsignedInt2UnsignedChar {
    unsigned char operator()( unsigned int val ) const { return val & 0xff; }
  };
  if ( out_color )
    PPMWriter<Image, UnsignedInt2Color>::exportPPM( out_fname.c_str(), output_u );
  else
    PGMWriter<Image,UnsignedInt2UnsignedChar>::exportPGM( out_fname.c_str(), output_u );
  trace.endBlock();


  return 0;
}
