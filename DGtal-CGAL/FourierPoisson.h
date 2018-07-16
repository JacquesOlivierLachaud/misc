#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/images/CImage.h>
#include <DGtal/io/Color.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/PPMWriter.h>
#include <DGtal/io/writers/PGMWriter.h>
#include <DGtal/math/RealFFT.h>

namespace DGtal {

  namespace image {
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
    struct ColorToGrayFunctor {
      int operator()( int color ) const
      {
	double val = 0.2126 * ( (double) ( ( color >> 16 ) & 0xff ) )
	  +          0.7152 * ( (double) ( ( color >> 8  ) & 0xff ) )
	  +          0.0722 * ( (double) ( ( color >> 0  ) & 0xff ) );
	return ( (int) round( val ) ) & 0xff;
      }
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
  
    struct Utils {
      typedef double                     Scalar;
      typedef PointVector< 3, Scalar >   Value;
      // typedef std::array< Value, 2 >     VectorValue;
      struct VectorValue {
	Value x;
	Value y;
      };
      typedef std::vector<Scalar>        ScalarForm;
      typedef std::vector<Value>         ValueForm;
      typedef std::vector<VectorValue>   VectorValueForm;
      
      template <typename Image>
      ScalarForm getScalarForm( const Image& image, bool color )
      {
	unsigned int v = 0;
	ScalarForm I( image.size() );
	if ( color ) {
	  auto F = ColorToGrayFunctor();
	  for ( unsigned int val : image )
	    I[ v++ ] = (Scalar) F( val );
	} else {
	  auto F = GrayToGrayFunctor();
	  for ( unsigned int val : image )
	    I[ v++ ] = (Scalar) F( val );
	}
	return I;
      }
      template <typename Image>
      ValueForm getValueForm( const Image& image, bool color )
      {
	unsigned int v = 0;
	ValueForm I( image.size() );
	if ( color ) {
	  auto R = ColorToRedFunctor();
	  auto G = ColorToGreenFunctor();
	  auto B = ColorToBlueFunctor();
	  for ( unsigned int val : image )
	    I[ v++ ] = Value( (Scalar) R( val ), (Scalar) G( val ), (Scalar) B( val ) );
	} else {
	  auto F = GrayToGrayFunctor();
	  for ( unsigned int val : image )
	    I[ v++ ] = Value( (Scalar) F( val ), (Scalar) F( val ), (Scalar) F( val ) );
	}
	return I;
      }
    };

    // ----------------------------------------------------------------------
    namespace functions {
      template <typename TImageDst, typename TImageSrc,
		typename TFunctor>
      void transform( TImageDst& dst, const TImageSrc& src,
		      TFunctor F )
      {
	std::transform( src.begin(), src.end(), dst.begin(), F );
      }

      template <typename TImage>
      void bandFilter( TImage& image,
		       double low_freq = 0.0, double high_freq = 0.5,
		       double shift = 0.0,
		       double output_min = 0.0,
		       double output_max = 255.0 )
      {
	typedef TImage Image;
	BOOST_CONCEPT_ASSERT(( DGtal::concepts::CImage< Image > ));

	typedef typename Image::Domain         Domain;
	typedef typename Domain::Space         Space;
	typedef typename Space::Point          Point;
	typedef typename Space::RealVector     RealVector;
	typedef typename RealVector::Component Scalar;
	typedef typename Image::Value          Value;

	Domain domain( image.domain() );
	RealFFT< Domain, Scalar > F ( domain );
	auto IF = F.getSpatialImage();
	auto FF = F.getFreqImage();
	// trace.info() << "spatial   domain = " << F.getSpatialDomain() << std::endl;
	// trace.info() << "frequency domain = " << F.getFreqDomain() << std::endl;
	auto IF_it = IF.begin();
	for ( auto I_it  = image.begin(), I_itE = image.end(); I_it != I_itE; I_it++ )
	  *IF_it++ = (Scalar) *I_it;
	F.forwardFFT();
	Scalar l2 = (Scalar) ( low_freq  * low_freq  );
	Scalar h2 = (Scalar) ( high_freq * high_freq );
	for ( auto&& fp : F.getFreqDomain() )
	  {
	    auto  freq = F.calcScaledFreqCoords( fp );
	    auto nfreq = freq.dot( freq );
	    if ( nfreq < l2 || nfreq > h2 )
	      FF.setValue( fp, (Scalar) 0 );
	  }
	F.backwardFFT();
	auto I_it  = image.begin();
	for ( auto IF_it = IF.begin(), IF_itE = IF.end(); IF_it != IF_itE; IF_it++ ) {
	  Value val = (Value) std::max( output_min,
					std::min( output_max, shift + (*IF_it) ) );
	  *I_it++ = val;
	}
      }      
    } // namespace functions
  } // namespace image
} // namespace DGtal
