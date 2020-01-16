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
#include <DGtal/images/SimpleThresholdForegroundPredicate.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/PPMWriter.h>
#include <DGtal/io/writers/PGMWriter.h>
#include <DGtal/io/writers/GenericWriter.h>

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace DGtal;
typedef Z2i::Space  Space;
typedef Z2i::Domain Domain;
typedef Z2i::Point  Point;
typedef ImageSelector < Z2i::Domain, unsigned int>::Type GenericImage;
typedef ImageSelector < Z2i::Domain, Color>::Type ColorImage;
typedef ImageSelector < Z2i::Domain, double>::Type DoubleImage;

void
makeRGBChannels( DoubleImage& R, DoubleImage& G, DoubleImage& B, const GenericImage& I )
{
  R = G = B = DoubleImage( I.domain() );
  auto itR = R.begin();
  auto itG = G.begin();
  auto itB = B.begin();
  for ( auto it  = I.begin(), itE = I.end(); it != itE; ++it )
    {
      const unsigned int val = *it;
      *itR++ = (val >> 16) & 0xff;
      *itG++ = (val >>  8) & 0xff;
      *itB++ = val & 0xff;
    }
}

struct RColor {
  double r, g, b;
  RColor( double _r = 0.0f, double _g = 0.0f, double _b = 0.0f )
    : r( _r ), g( _g ), b( _b ) {}

  double dot( const RColor& other ) const
  {
    return r * other.r + g * other.g + b * other.b;
  }

  double squaredNorm() const
  {
    return this->dot( *this );
  }

  double norm() const
  {
    return sqrt( squaredNorm() );
  }

  RColor& operator+=( const RColor& other )
  {
    r += other.r;
    g += other.g;
    b += other.b;
    return *this;
  }

  RColor operator-( const RColor& other ) const
  {
    return RColor( r - other.r, g - other.g, b - other.b );
  }

  RColor operator/( double v ) const
  {
    return RColor( r / v, g / v, b / v );
  }
};

std::ostream& operator<<( std::ostream& out, RColor color )
{
  out << "(" << color.r << ", " << color.g << ", " << color.b << ")";
  return out;
}

double rand01()
{
  return (double) rand() / (double) RAND_MAX;
}
int irandk( int k )
{
  return (int) floor( rand01() * (double) k );
}

// Sliced transport along direction \a dir of input1 toward input2
// @return the total cost
double slicedTransport( std::vector<int>& indices1,
			std::vector<int>& indices2,
			const RColor dir,
			const std::vector<RColor>& input1,
			const std::vector<RColor>& input2 )
{
  const int n = input1.size();
  indices1.resize( n );
  indices2.resize( n );
  for ( int i = 0; i < n; ++i )
    indices1[ i ] = indices2[ i ] = i;
  std::vector<double> p1( n );
  std::vector<double> p2( n );
  for ( int i = 0; i < n; ++i )
    {
      p1[ i ] = input1[ i ].dot( dir );
      p2[ i ] = input2[ i ].dot( dir );
    }
  std::sort( indices1.begin(), indices1.end(),
	     [&p1](int i, int j) { return p1[ i ] < p1[ j ]; } );
  std::sort( indices2.begin(), indices2.end(),
	     [&p2](int i, int j) { return p2[ i ] < p2[ j ]; } );
  std::vector<RColor> sorted1( n );
  std::vector<RColor> sorted2( n );
  double cost = 0.0;
  for ( int i = 0; i < n; ++i )
    {
      sorted1[ indices1[ i ] ] = input1[ i ];
      sorted2[ indices2[ i ] ] = input2[ i ];
    }
  for ( int i = 0; i < n; ++i )
    cost += ( sorted1[ i ] - sorted2[ i ] ).squaredNorm();
  return cost;
}

double bestSlicedTransport( std::vector<int>& indices1,
			    std::vector<int>& indices2,
			    const std::vector<RColor>& input1,
			    const std::vector<RColor>& input2,
			    int nb = 100 )
{
  double best_cost = std::numeric_limits<double>::infinity();
  RColor best_dir;
  for ( int i = 0; i < nb; ++i )
    {
      RColor dir( rand01()*2.0 - 1.0, rand01()*2.0 - 1.0, rand01()*2.0 - 1.0 );
      dir = dir / dir.norm();
      double cost = slicedTransport( indices1, indices2, dir, input1, input2 );
      if ( cost < best_cost )
	{
	  best_dir  = dir;
	  best_cost = cost;
          trace.info() << "[" << i << "]"
                       << " dir=" << dir
                       << " cost=" << cost
                       << " best_cost=" << best_cost << std::endl;
	}
    }
  double cost = slicedTransport( indices1, indices2, best_dir, input1, input2 );
  return best_cost;
}


struct ImagePartition
{
  DoubleImage R;
  DoubleImage G;
  DoubleImage B;
  GenericImage L;
  std::vector<Point>  Rep;
  std::vector<RColor> RepColorSum;
  std::vector<int>    RepColorNb;
  
  typedef std::tuple<Point,double,int> WeightedPoint;
  struct WeightedPointComparator {
    bool operator()( WeightedPoint p, WeightedPoint q ) const
    {
      return std::get<1>( p ) < std::get<1>( q )
				|| ( std::get<1>( p ) == std::get<1>( q ) &&
				     std::get<0>( p ) < std::get<0>( q ) );
    }
  };
  
  ImagePartition( const GenericImage& I )
    : R( I.domain() ), G( I.domain() ), B( I.domain() ), L( I.domain() )
  {
    makeRGBChannels( R, G, B, I );
  }

  RColor getRepColor( int i ) const
  {
    return RepColorSum[ i ] / (double) RepColorNb[ i ];
  }
  RColor getColor( Point p ) const
  {
    return RColor( R( p ), G( p ), B( p ) );
  }

  double getWeight( Point p, int i ) const
  {
    return -( getColor( p ) - getRepColor( i ) ).squaredNorm() * RepColorNb[ i ];
  }
  void partition( int k )
  {
    choosePoints( k );
    std::priority_queue< WeightedPoint,
			 std::vector<WeightedPoint>,
			 WeightedPointComparator > Q;
    for ( int i = 1; i < Rep.size(); ++i )
      Q.emplace( make_tuple( Rep[ i ], 1.0, i ) );
    while ( ! Q.empty() )
      {
	auto wp = Q.top(); Q.pop();
	Point p = std::get<0>( wp );
	int   i = std::get<2>( wp );
	if ( L( p ) != 0 ) continue;
	L.setValue( p, i );
	RepColorSum[ i ] += getColor( p );
	RepColorNb[ i ]  += 1;
	std::array<Point,4> neighbors
	{ Point( p[ 0 ] - 1, p[ 1 ] ), Point( p[ 0 ] + 1, p[ 1 ] ),
	    Point( p[ 0 ], p[ 1 ] - 1 ), Point( p[ 0 ], p[ 1 ] + 1 ) };
	for ( auto q : neighbors )
	  if ( L.domain().isInside( q ) && ( L( q ) == 0 ) )
	    Q.emplace( std::make_tuple( q, getWeight( q, i ), i ) );
      }
  }

  std::vector<RColor> makeRepColorVector() const
  {
    std::vector<RColor> RepColor( Rep.size() - 1 );
    for ( int i = 0; i < RepColor.size(); ++i )
      {
	RepColor[ i ] = getRepColor( i+1 );
      }
    return RepColor;
  }
  
  void choosePoints( int k )
  {
    Rep.resize( k*k+1 );
    RepColorSum = std::vector<RColor>( Rep.size(), RColor() );
    RepColorNb  = std::vector<int>   ( Rep.size(), 0 );
    int w = L.extent()[ 0 ];
    int h = L.extent()[ 1 ];
    int xstep = w / k;
    int ystep = h / k;
    for ( int i = 0; i < k; ++i ) 
      for ( int j = 0; j < k; ++j )
	{
	  int x = xstep * j + irandk( xstep );
	  int y = ystep * i + irandk( ystep );
	  Rep[ i*k + j + 1 ] = Point( std::min( x, w-1 ), std::min( y, h-1 ) );
	}
    for ( auto& v : L ) v = 0;
  }

  void debug_info()
  {
    for ( int i = 1; i < Rep.size(); i++ )
      trace.info() << "R[" << i << "]=" << Rep[ i ] << " #=" << RepColorNb[ i ]
		   << " col=" << getRepColor( i ) << std::endl;
  }

};

void transport( const GenericImage& image1,
		const GenericImage& image2,
		int samples )
{
  ImagePartition IP1( image1 );
  ImagePartition IP2( image2 );
  IP1.partition( samples );
  IP2.partition( samples );
  IP1.debug_info();
  IP2.debug_info();
  auto input1 = IP1.makeRepColorVector();
  auto input2 = IP2.makeRepColorVector();
  std::vector<int> indices1;
  std::vector<int> indices2;
  bestSlicedTransport( indices1, indices2, input1, input2 );
}

int main( int argc, char** argv )
{

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input1,1", po::value<std::string>(), "Specifies the first input image PPM filename.")
    ("input2,2", po::value<std::string>(), "Specifies the second input image PPM filename.")
    ("output,o", po::value<std::string>()->default_value("out"), "Specifies the basename of output color transfered images." )
    ("samples,n", po::value<int>()->default_value( 20 ), "Defines the number of superpixels per row and per column" );

  bool parseOK = true;
  po::variables_map vm;
  try {
    po::store( po::parse_command_line(argc, argv, general_opt), vm );  
  } catch ( const std::exception& ex ) {
    parseOK = false;
    trace.info() << "Error checking program options: " << ex.what() << std::endl;
  }
    
  po::notify(vm);    
  if( ! parseOK || vm.count("help") || argc <= 1
      || (! vm.count( "input1") )
      || (! vm.count( "input2") ) )
    {
      trace.info()<< "Performs a color transfer between the two images."
		  <<std::endl << "Basic usage: " << std::endl
		  << "\tcolor-transfer [options] -1 <image.ppm> -2 <image.ppm> "
		  <<std::endl << general_opt << "\n";
      return 0;
    }

  trace.beginBlock("Loading images");
  std::string input1 = vm[ "input1" ].as<std::string>();
  std::string input2 = vm[ "input2" ].as<std::string>();
  const GenericImage image1 = GenericReader<GenericImage>::import( input1 ); 
  const GenericImage image2 = GenericReader<GenericImage>::import( input2 ); 
  std::string   ext1 = input1.substr(input1.find_last_of(".") + 1);
  std::string   ext2 = input2.substr(input2.find_last_of(".") + 1);
  bool        color1 = ( ext1 == "ppm" );
  bool        color2 = ( ext2 == "ppm" );
  trace.info() << std::fixed
	       << "Image1 <" << input1
	       << "> size=" << image1.extent()[ 0 ]
	       << "x" << image1.extent()[ 1 ]
	       << " color=" << ( color1 ? "True" : "False" ) << std::endl;
  trace.info() << "Image2 <" << input2
	       << "> size=" << image2.extent()[ 0 ]
	       << "x" << image2.extent()[ 1 ]
	       << " color=" << ( color2 ? "True" : "False" ) << std::endl;
  trace.endBlock();
  if ( ! ( color1 && color2 ) )
    {
      trace.info() << "Requires color images." << std::endl;
      return 0;
    }
  trace.beginBlock("Partition images");
  int samples = vm[ "samples" ].as<int>();
  transport( image1, image2, samples );
  trace.endBlock();
  return 0;
}
  
