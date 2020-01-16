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
makeRGBChannels( DoubleImage& R, DoubleImage& G, DoubleImage& B,
                 const GenericImage& I )
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
void
makeImageFromRGBChannels( GenericImage& I,
                          const DoubleImage& R, const DoubleImage& G, const DoubleImage& B )
{
  I = GenericImage( R.domain() );
  auto itR = R.begin();
  auto itG = G.begin();
  auto itB = B.begin();
  for ( auto it  = I.begin(), itE = I.end(); it != itE; ++it )
    {
      unsigned int r = (unsigned int) round( std::max( 0.0, std::min( 255.0, *itR++ ) ) );
      unsigned int g = (unsigned int) round( std::max( 0.0, std::min( 255.0, *itG++ ) ) );
      unsigned int b = (unsigned int) round( std::max( 0.0, std::min( 255.0, *itB++ ) ) );
      *it = ( r << 16 ) + ( g << 8 ) + b;
    }
}

struct RColor {
  double r, g, b;
  RColor( double _r = 0.0, double _g = 0.0, double _b = 0.0 )
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

  RColor operator+( const RColor& other ) const
  {
    return RColor( r + other.r, g + other.g, b + other.b );
  }
  RColor operator-( const RColor& other ) const
  {
    return RColor( r - other.r, g - other.g, b - other.b );
  }

  RColor operator/( double v ) const
  {
    return RColor( r / v, g / v, b / v );
  }
  RColor operator*( double v ) const
  {
    return RColor( r * v, g * v, b * v );
  }
};

RColor operator*( double v, RColor col )
{
  return RColor( col.r * v, col.g * v, col.b * v );
}

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

void interpolateColors( std::vector<RColor>&       output,
                        std::vector<int>&          in2out,
                        const std::vector<RColor>& input )
{
  const auto  n_in  = input.size();
  const auto  n_out = output.size();
  in2out.resize( n_in );
  // trace.info() << "IC #in=" << n_in << " #out=" << n_out << std::endl;
  const double coef = (double) n_in / (double) n_out;
  for ( int i = 0; i < n_out; ++i )
    {
      double x = (double) i * coef;
      int    k = (int) floor( x );
      double t = x - k;
      output[ i ] = ( 1.0 - t ) * input[ k ] + t * ( (k+1) < n_in ? input[ k+1 ] : input[ k ] );
      in2out[ k ] = i;
      // trace.info() << " " << output[ i ];
    }
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
  // debug
  // std::vector<double> sorted_p1( n );
  // std::vector<double> sorted_p2( n );
  // for ( int i = 0; i < n; ++i ) sorted_p1[ i ] = p1[ indices1[ i ] ];
  // for ( int i = 0; i < n; ++i ) sorted_p2[ i ] = p2[ indices2[ i ] ];
  // for ( int i = 0; i < n - 1; ++i ) {
  //   trace.info() << sorted_p1[ i ] << std::endl;
  //   if ( sorted_p1[ i ] > sorted_p1[ i+1 ] ) trace.error() << "bad sort p1" << std::endl;
  // }
  // for ( int i = 0; i < n - 1; ++i ) {
  //   trace.info() << sorted_p1[ i ] << std::endl;
  //   if ( sorted_p2[ i ] > sorted_p2[ i+1 ] ) trace.error() << "bad sort p2" << std::endl;
  // }
  std::vector<RColor> sorted1( n );
  std::vector<RColor> sorted2( n );
  double cost = 0.0;
  for ( int i = 0; i < n; ++i )
    {
      //sorted1[ indices1[ i ] ] = input1[ i ];
      //sorted2[ indices2[ i ] ] = input2[ i ];
      sorted1[ i ] = input1[ indices1[ i ] ];
      sorted2[ i ] = input2[ indices2[ i ] ];
    }
  for ( int i = 0; i < n; ++i )
    cost += ( sorted1[ i ] - sorted2[ i ] ).squaredNorm();
  //cost += pow( ( sorted1[ i ] - sorted2[ i ] ).dot( dir ), 2.0 );
  return cost;
}

// Unbalanced sliced transport along direction \a dir of input1 toward input2
// @return the total cost
double unbalancedSlicedTransport( std::vector<int>& indices1,
                                  std::vector<int>& indices2,
                                  std::vector<RColor>& output1,
                                  std::vector<RColor>& output2,
                                  const RColor dir,
                                  const std::vector<RColor>& input1,
                                  const std::vector<RColor>& input2 )
{
  const int n1 = input1.size();
  const int n2 = input2.size();
  indices1.resize( n1 );
  indices2.resize( n2 );
  for ( int i = 0; i < n1; ++i ) indices1[ i ] = i;
  for ( int i = 0; i < n2; ++i ) indices2[ i ] = i;
  std::vector<double> p1( n1 );
  std::vector<double> p2( n2 );
  for ( int i = 0; i < n1; ++i ) p1[ i ] = input1[ i ].dot( dir );
  for ( int i = 0; i < n2; ++i ) p2[ i ] = input2[ i ].dot( dir );
  std::sort( indices1.begin(), indices1.end(),
	     [&p1](int i, int j) { return p1[ i ] < p1[ j ]; } );
  std::sort( indices2.begin(), indices2.end(),
	     [&p2](int i, int j) { return p2[ i ] < p2[ j ]; } );
  const int FACTOR = 3;
  const int n      = FACTOR * std::max( n1, n2 );
  std::vector<RColor> sorted1( n1 );
  std::vector<RColor> sorted2( n2 );
  for ( int i = 0; i < n1; ++i ) sorted1[ indices1[ i ] ] = input1[ i ];
  for ( int i = 0; i < n2; ++i ) sorted2[ indices2[ i ] ] = input2[ i ];
  output1.resize( n );
  output2.resize( n );
  std::vector<int> in2out1;
  std::vector<int> in2out2;
  interpolateColors( output1, in2out1, sorted1 );
  interpolateColors( output2, in2out2, sorted2 );
  double cost = 0.0;
  for ( int i = 0; i < n; ++i )
    cost += ( output1[ i ] - output2[ i ] ).squaredNorm();
  for ( int i = 0; i < n1; ++i ) indices1[ i ] = in2out1[ indices1[ i ] ];
  for ( int i = 0; i < n2; ++i ) indices2[ i ] = in2out2[ indices2[ i ] ];
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
  return cost;
}

double bestUnbalancedSlicedTransport( std::vector<int>& indices1,
                                      std::vector<int>& indices2,
                                      std::vector<RColor>& output1,
                                      std::vector<RColor>& output2,
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
      double cost = unbalancedSlicedTransport( indices1, indices2,
                                               output1, output2,
                                               dir, input1, input2 );
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
  double cost = unbalancedSlicedTransport( indices1, indices2,
                                           output1, output2,
                                           best_dir, input1, input2 );
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
  std::vector< std::vector<Point> >  all_points;
  std::vector< std::vector<RColor> > all_colors;
  
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

  void getImage( GenericImage& I )
  {
    makeImageFromRGBChannels( I, R, G, B );
  }

  
  RColor getRepColor( int i ) const
  {
    return RepColorSum[ i ] / (double) RepColorNb[ i ];
  }

  RColor getColor( Point p ) const
  {
    return RColor( R( p ), G( p ), B( p ) );
  }

  void setColor( Point p, RColor col )
  {
    R.setValue( p, col.r );
    G.setValue( p, col.g );
    B.setValue( p, col.b );
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
    computeAllVectors();
  }

  void computeAllVectors()
  {
    all_points.resize( Rep.size() );
    all_colors.resize( Rep.size() );
    for ( int i = 0; i < Rep.size(); ++i )
      {
        all_points[ i ].clear();
        all_colors[ i ].clear();
      }
    for ( auto p : L.domain() )
      {
        int    i = L( p );
        RColor c = getColor( p );
        all_points[ i ].push_back( p );
        all_colors[ i ].push_back( c );
      }
  }


  std::vector<RColor> getColorVector() const
  {
    std::vector<RColor> colors;
    for ( auto p : L.domain() )
      colors.push_back( getColor( p ) );
    return colors;
  }

  void setColorVector( const std::vector<int>& indices,
                       const std::vector<RColor>& colors )
  {
    std::vector<RColor> sorted_colors( colors.size() );
    for ( int i = 0; i < indices.size(); ++i )
      sorted_colors[ i ] = colors[ indices[ i ] ];
    //sorted_colors[ indices[ i ] ] = colors[ i ];
    int i = 0;
    for ( auto p : L.domain() )
      setColor( p, sorted_colors[ i++ ] );
  }
  
  void setColorVector( const std::vector<int>& indices1,
                       const std::vector<int>& indices2,
                       const std::vector<RColor>& colors2 )
  {
    std::vector<RColor> sorted_colors2( colors2.size() );
    for ( int i = 0; i < indices2.size(); ++i )
      sorted_colors2[ i ] = colors2[ indices2[ i ] ];
    //sorted_colors2[ indices2[ i ] ] = colors2[ i ];
    std::vector<int> rev1( indices1.size() );
    for ( int i = 0; i < indices1.size(); ++i )
      rev1[ indices1[ i ] ] = i;
    int i = 0;
    for ( auto p : L.domain() )
      setColor( p, sorted_colors2[ rev1[ i++ ] ] );
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

  void applyColors( int idx, const std::vector<int>& local_indices,
                    const std::vector<RColor>&       new_colors )
  {
    int j = 0;
    for ( auto p : all_points[ idx ] )
        setColor( p, new_colors[ local_indices[ j++ ] ] );
  }
  
};

std::pair< GenericImage, GenericImage >
transport( const GenericImage& image1,
           const GenericImage& image2,
           int samples )
{
  ImagePartition IP1( image1 );
  ImagePartition IP2( image2 );
  // Computes super-pixels
  IP1.partition( samples );
  IP2.partition( samples );
  // IP1.debug_info();
  // IP2.debug_info();
  // Computes best sliced transport between super-pixels
  auto input1 = IP1.makeRepColorVector();
  auto input2 = IP2.makeRepColorVector();
  std::vector<int> indices1;
  std::vector<int> indices2;
  bestSlicedTransport( indices1, indices2, input1, input2 );
  // Computes best unbalanced sliced transport within super-pixels
  for ( int k = 0; k < indices1.size(); ++k )
    {
      const int i1 = indices1[ k ] + 1;
      const int i2 = indices2[ k ] + 1;
      std::vector<int> local_indices1;
      std::vector<int> local_indices2;
      std::vector<RColor> output1;
      std::vector<RColor> output2;
      double cost = bestUnbalancedSlicedTransport( local_indices1, local_indices2, 
                                                   output1, output2,
                                                   IP1.all_colors[ i1 ], IP2.all_colors[ i2 ] );
      IP1.applyColors( i1, local_indices1, output2 );
      IP2.applyColors( i2, local_indices2, output1 );
      trace.info() << "[" << k << "] " << i1 << " <-> " << i2 << " cost=" << cost << std::endl;
    }
  GenericImage Out1( image1.domain() );
  GenericImage Out2( image2.domain() );
  IP1.getImage( Out1 );
  IP2.getImage( Out2 );
  return std::make_pair( Out1, Out2 ); 
}

std::pair< GenericImage, GenericImage >
simpleTransport( const GenericImage& image1,
                 const GenericImage& image2 )
{
  auto D1 = image1.domain();
  auto D2 = image2.domain();
  auto D  = Domain( D1.lowerBound().sup( D2.lowerBound() ),
                    D1.upperBound().inf( D2.upperBound() ) );
  GenericImage image1bis( D );
  GenericImage image2bis( D );
  for ( auto p : D ) image1bis.setValue( p, image1( p ) );
  for ( auto p : D ) image2bis.setValue( p, image2( p ) );
  ImagePartition IP1( image1bis );
  ImagePartition IP2( image2bis );
  auto input1 = IP1.getColorVector();
  auto input2 = IP2.getColorVector();
  std::vector<int> indices1;
  std::vector<int> indices2;
  bestSlicedTransport( indices1, indices2, input1, input2 );
  IP1.setColorVector( indices1, indices2, input2 );
  IP2.setColorVector( indices2, indices1, input1 );

  GenericImage Out1( image1bis.domain() );
  GenericImage Out2( image2bis.domain() );
  IP1.getImage( Out1 );
  IP2.getImage( Out2 );
  return std::make_pair( Out1, Out2 ); 
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
  std::string output = vm[ "output" ].as<std::string>();
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
  auto images = ( samples == 0 )
    ? simpleTransport( image1, image2 )
    : transport( image1, image2, samples );
  struct UnsignedInt2Color {
    Color operator()( unsigned int val ) const { return Color( val ); }
  };
  PPMWriter<GenericImage, UnsignedInt2Color>::exportPPM( ( output + "-1.ppm" ).c_str(),
                                                         images.first );
  PPMWriter<GenericImage, UnsignedInt2Color>::exportPPM( ( output + "-2.ppm" ).c_str(),
                                                         images.second );
  trace.endBlock();
  return 0;
}
  
