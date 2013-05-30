#pragma once
#ifndef _DGTAL_AUXILIARY_H_
#define _DGTAL_AUXILIARY_H_

#include <string>
#include <sstream>
#include <fstream>

#include "DGtal/base/Common.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/io/readers/VolReader.h"


namespace DGtal 
{

  /**
     Creates a digital Khalimsky space and a set of surfels that is the boundary of the thresholded .vol image.
  */
  template <typename KSpace, typename SCellSet>
  int
  makeSpaceAndBoundaryFromVolFile
  ( KSpace & ks,
    SCellSet & boundary,
    std::string input, int minThreshold, int maxThreshold )
  {
    typedef typename KSpace::Space Z3;
    typedef HyperRectDomain<Z3> Domain;
    typedef typename ImageSelector < Domain, int>::Type MyImage;
    typedef typename DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet;
    
    trace.beginBlock( "Reading vol file into an image." );
    MyImage image = VolReader<MyImage>::importVol( input );
    DigitalSet set3d ( image.domain() );
    DGtal::SetFromImage<DigitalSet>
      ::append( set3d, image, minThreshold, maxThreshold );
    trace.endBlock();
    trace.beginBlock( "Creating space." );
    bool space_ok = ks.init( image.domain().lowerBound(),
			     image.domain().upperBound(), true );
    if (!space_ok)
      {
	trace.error() << "Error in the Khamisky space construction." << std::endl;
	return 2;
      }
    
    trace.endBlock();
    trace.beginBlock( "Extracting boundary by scanning the space. " );
    Surfaces<KSpace>::sMakeBoundary( boundary,
				 ks, set3d,
				 image.domain().lowerBound(),
				 image.domain().upperBound() );
    trace.info() << "Digital surface has " << boundary.size() << " surfels."
		 << std::endl;
    trace.endBlock();
    return 0;
  }

  /**
     Creates a digital Khalimsky space and a set of surfels that
     approximate the zero-level set of the given polynomial surface.
  */
  template <typename KSpace, typename SCellSet>
  int
  makeSpaceAndBoundaryFromPolynomialString
  ( KSpace & ks,
    SCellSet & boundary,
    std::string poly_str, double h, int b )
  {
    typedef typename KSpace::Space Z3;
    typedef typename KSpace::SCell SCell;
    typedef typename KSpace::Vector Vector;
    typedef typename Z3::RealPoint RealPoint;
    typedef typename RealPoint::Coordinate Ring;

    typedef MPolynomial<3, Ring> Polynomial3;
    typedef MPolynomialReader<3, Ring> Polynomial3Reader;
    typedef ImplicitPolynomial3Shape<Z3> Shape;
    typedef GaussDigitizer<Z3,Shape> Digitizer;

    trace.beginBlock( "Reading polynomial." );
    Polynomial3 P;
    Polynomial3Reader reader;
    std::string::const_iterator iter 
      = reader.read( P, poly_str.begin(), poly_str.end() );
    if ( iter != poly_str.end() )
      {
	trace.error() << "I read only <" 
		      << poly_str.substr( 0, iter - poly_str.begin() )
		      << ">, and I built P=" << P << std::endl;
	return 3;
      }
    trace.info() << "P( X_0, X_1, X_2 ) = " << P << std::endl;
    trace.endBlock();
    
    trace.beginBlock( "Extract polynomial surface by tracking." );
    SurfelAdjacency<3> sAdj( true );
    Shape shape( P );
    Digitizer dig;
    dig.attach( shape );
    dig.init( Vector::diagonal( -b ),
	      Vector::diagonal(  b ), h ); 
    bool space_ok = ks.init( dig.getLowerBound(), dig.getUpperBound(), true );
    if (!space_ok)
      {
	trace.error() << "Error in the Khalimsky space construction." << std::endl;
	return 2;
      }
    SCell bel = Surfaces<KSpace>::findABel( ks, dig, 100000 );
    trace.info() << "initial bel (Khalimsky coordinates): " << ks.sKCoords( bel ) << std::endl;
    Surfaces<KSpace>::trackBoundary( boundary, ks, sAdj, dig, bel );
    trace.info() << "Digital surface has " << boundary.size() << " surfels."
		 << std::endl;
    trace.endBlock();
    return 0;
  }


  /**
     Creates a digital Khalimsky space and a set of surfels that
     approximate the zero-level set of the given polynomial surface.
  */
  template <typename KSpace, typename SCellSet>
  int
  makeSpaceAndBoundaryFromPointList
  ( KSpace & ks,
    SCellSet & boundary,
    std::string point_list_filename )
  {
    typedef typename KSpace::Point Point;
    typedef typename KSpace::SCell SCell;
    typedef typename Point::Coordinate Coordinate;
    trace.beginBlock( "Reading point list from file." );
    std::ifstream inputfile( point_list_filename.c_str() );
    std::string str;
    std::vector<Point> points;
    Point lower, upper;
    while ( ( ! inputfile.eof() ) && inputfile.good() )
      {
        std::getline( inputfile, str );
        if ( ( ! str.empty() ) && ( str[ 0 ] != '#' ) )
          {
            Point p;
            std::istringstream istr( str );
            Coordinate x, y, z;
            //istr >> x >> y >> z;
            istr >> p[ 0 ] >> p[ 1 ] >> p[ 2 ];
            lower = points.empty() ? p : lower.inf( p );
            upper = points.empty() ? p : upper.sup( p );
            points.push_back( p );
          }
      }
    bool ok = inputfile.eof();
    trace.endBlock();
    lower -= Point::diagonal( 1 );
    upper += Point::diagonal( 1 );
    if ( ! ks.init( lower, upper, true ) )
      return 2;
    for ( typename std::vector<Point>::const_iterator it = points.begin(), itend = points.end();
          it != itend; ++it )
      {
        SCell v = ks.sSpel( *it );
        for ( Dimension i = 0; i < KSpace::dimension; ++i )
          {
            boundary.insert( ks.sIncident( v, i, true ) );
            boundary.insert( ks.sIncident( v, i, false ) );
          }
      }
    return 0;
  }

  template <typename KSpace, typename SCellSetConstIterator>
  int
  getInnerVoxelCoordinates( std::vector< typename KSpace::Point > & points,
			    const KSpace & ks,
			    SCellSetConstIterator it, SCellSetConstIterator ite )
  {
    typedef typename KSpace::Cell Cell;
    trace.beginBlock("Getting inner voxel coordinates.");
    std::set<Cell> inner_points;
    for( ; it != ite; ++it )
      // Get inner point.
      inner_points.insert( ks.unsigns( ks.sDirectIncident( *it, ks.sOrthDir( *it ) ) ) );
    for ( typename std::set<Cell>::const_iterator itS = inner_points.begin(), itSEnd = inner_points.end();
	  itS != itSEnd; ++itS )
      points.push_back( ks.uCoords( *itS ) );
    trace.endBlock();
    return 0;
  }

  template <typename KSpace, typename SCellSetConstIterator>
  int
  getInnerAndOuterVoxelCoordinates( std::vector< typename KSpace::Point > & points,
				    const KSpace & ks,
				    SCellSetConstIterator it, SCellSetConstIterator ite )
  {
    typedef typename KSpace::Cell Cell;
    trace.beginBlock("Getting inner and outer voxel coordinates.");
    std::set<Cell> inner_points;
    for( ; it != ite; ++it )
      { 
	// Get inner point.
	inner_points.insert( ks.unsigns( ks.sDirectIncident( *it, ks.sOrthDir( *it ) ) ) );
	// Get Outer point.
	inner_points.insert( ks.unsigns( ks.sIndirectIncident( *it, ks.sOrthDir( *it ) ) ) );
      }
    for ( typename std::set<Cell>::const_iterator itS = inner_points.begin(), itSEnd = inner_points.end();
	  itS != itSEnd; ++itS )
      points.push_back( ks.uCoords( *itS ) );
    trace.endBlock();
    return 0;
  }

  template <typename KSpace, typename SCellSetConstIterator>
  int
  getPointelCoordinates( std::vector< typename KSpace::Point > & points,
			 const KSpace & ks,
			 SCellSetConstIterator it, SCellSetConstIterator ite )
  {
    typedef typename KSpace::SCell SCell;
    typedef typename KSpace::Cell Cell;
    trace.beginBlock("Getting pointels coordinates.");
    std::set<Cell> pointels;
    for( ; it != ite; ++it )
      // Get 4 pointels
      {
	DGtal::Dimension k = ks.sOrthDir( *it );
	DGtal::Dimension i = (k+1)%3;
	DGtal::Dimension j = (k+2)%3;
	SCell linel_i0 = ks.sIncident( *it, i, false );
	SCell linel_i1 = ks.sIncident( *it, i, true );
	pointels.insert( ks.unsigns( ks.sIncident( linel_i0, j, false ) ) );
	pointels.insert( ks.unsigns( ks.sIncident( linel_i0, j, true ) ) );
	pointels.insert( ks.unsigns( ks.sIncident( linel_i1, j, false ) ) );
	pointels.insert( ks.unsigns( ks.sIncident( linel_i1, j, true ) ) );
      }
    for ( typename std::set<Cell>::const_iterator itS = pointels.begin(), itSEnd = pointels.end();
	  itS != itSEnd; ++itS )
      points.push_back( ks.uCoords( *itS ) );
    trace.endBlock();
    return 0;
  }

  /**
     Given a set of points \a P, computes the set of \f$ Z^d \setminus P \f$ that touches P.
  */
  template <typename Space>
  void computeBorder( std::set< typename Space::Point > & Q, 
                      const std::set< typename Space::Point > & P )
  {
    typedef typename Space::Point Point;
    typedef typename std::set< Point > PointSet;
    typedef typename std::set< Point >::const_iterator PointSetConstIterator;
    typedef typename Point::UnsignedComponent Size;
    typedef HyperRectDomain< Space > Domain;
    typedef typename Domain::ConstIterator DomainConstIterator;
    Point d1 = Point::diagonal( 1 );
    Q.clear();
    // Scan surroundings of all points of P.
    for ( PointSetConstIterator it = P.begin(), itend = P.end(); 
          it != itend; ++it )
      {
        Domain localDomain( *it - d1, *it + d1 );
        for ( DomainConstIterator itd = localDomain.begin(), itdend = localDomain.end();
              itd != itdend; ++itd )
          {
            Point p = *itd;
            if ( P.find( p ) == P.end() ) Q.insert( p );
          }
      }
  }

}


#endif //_DGTAL_AUXILIARY_H_
