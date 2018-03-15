/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file ImageTriangulation.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2018/01/26
 *
 * Header file for module ImageTriangulation.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ImageTriangulation_RECURSES)
#error Recursive header files inclusion detected in ImageTriangulation.h
#else // defined(ImageTriangulation_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ImageTriangulation_RECURSES

#if !defined ImageTriangulation_h
/** Prevents repeated inclusion of headers. */
#define ImageTriangulation_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <set>
#include <DGtal/base/Common.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/TriangulatedSurface.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class ImageTriangulation
  /**
     Description of class 'ImageTriangulation' <p> \brief 
  */
  template <typename TSpace, int M>
  class ImageTriangulation
  {
  public:
    typedef TSpace                            Space;
    typedef ImageTriangulation< TSpace, M> Self;
    BOOST_CONCEPT_ASSERT(( concepts::CSpace< Space > ));
    BOOST_STATIC_ASSERT (( Space::dimension <= 3 ));
    BOOST_STATIC_ASSERT (( M >= 1 ));

    typedef typename Space::Point                 Point;
    typedef typename Space::Vector                Vector;
    typedef typename Space::RealPoint             RealPoint;
    typedef typename Space::RealVector            RealVector;
    typedef typename Space::Size                  Size;
    typedef HyperRectDomain<Space>                Domain;
    typedef typename RealVector::Component        Scalar;
    typedef PointVector< M, Scalar >              Value;
    typedef std::array< Value, Space::dimension > VectorValue;
    ///static constants to store the dimension.
    static const Dimension N = Space::dimension;
    typedef std::vector<Scalar>                   ScalarForm;
    typedef std::vector<Value>                    ValueForm;
    typedef std::vector<VectorValue>              VectorValueForm;
    typedef TriangulatedSurface<RealPoint>        Triangulation;
    typedef typename Triangulation::VertexIndex   VertexIndex;

    struct Pixel {
      Point  xy;
      Value  v;
      Scalar prob;
      bool operator<( const Pixel& other ) const {
	return ( prob > other.prob )
		    || ( ( prob == other.prob ) && ( xy < other.xy ) );
      }
    };
    
    Domain          _domain;
    ValueForm       _I;
    std::set<Pixel> _subset;
    Triangulation   _T;
    
    // ----------------------- Standard services ------------------------------
  public:
  
    /**
     * Destructor.
     */
    ~ImageTriangulation() {}

    /// Default constructor. The object is invalid.
    ImageTriangulation() : _domain( Point(), Point() ) {}

    void init( Domain aDomain, const ValueForm& I, const ScalarForm& L,
	       Scalar compression = 0.05,
	       std::function< Scalar(Scalar) > P = [] ( Scalar x ) { return x; } )
    {
      trace.info() << "[ImageTriangulation::init] extract significant points."
		   << std::endl;
      _domain = aDomain;
      _subset.clear();
      std::vector< Pixel > pixels;
      VertexIndex vtx = 0;
      for ( auto p : _domain ) {
	pixels.push_back( Pixel { p, I[ vtx ], P( L[ vtx ] ) } );
	vtx++;
      }
      // Add extremal points
      Size w = ( _domain.upperBound() - _domain.lowerBound() )[ 0 ] + 1;
      Size h = ( _domain.upperBound() - _domain.lowerBound() )[ 1 ] + 1;
      _subset.insert( pixels[ 0 ] );
      _subset.insert( pixels[ w-1 ] );
      _subset.insert( pixels[ w*(h-1) ] );
      _subset.insert( pixels[ w*h-1 ] );
      std::sort( pixels.begin(), pixels.end() );
      // Compute cumulative density function
      ScalarForm cdf( pixels.size() );
      Scalar sum = 0.0;
      for ( Size i = 0; i < pixels.size(); ++i ) {
	sum += pow( pixels[ i ].prob, 1 ) + 1;
	cdf[ i ] = sum;
      }
      for ( Size i = 0; i < pixels.size(); ++i ) cdf[ i ] /= sum; 
      // Extracts a subset of points.
      Size nb_points = (Size) round( I.size() * compression );
      while ( _subset.size() < nb_points ) {
	Scalar x = rand() / (Scalar) RAND_MAX;
	auto  it = lower_bound( cdf.begin(), cdf.end(), x );
	Size   i = it - cdf.begin();
	_subset.insert( pixels[ i ] );
	// std::cout << "i=" << i << " nb=" << _subset.size()
	// 	  << "/" << nb_points << "/" << I.size() << std::endl;
      }
      // Compute Delaunay Triangulation
      trace.info() << "[ImageTriangulation::init] Compute Delaunay Triangulation."
		   << std::endl;
      typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
      typedef CGAL::Delaunay_triangulation_2<Kernel>           DelTriangulation;
      typedef DelTriangulation::Finite_faces_iterator          FiniteFacesIterator;
      typedef DelTriangulation::Face_handle                    FaceHandle;
      typedef DelTriangulation::Point                          CGALPoint2;
      typedef DelTriangulation::Vertex_handle                  VertexHandle;
      /// Build Delaunay triangulation.
      DelTriangulation DT;
      for ( auto pixel : _subset ) {
	CGALPoint2 pt( pixel.xy[ 0 ], pixel.xy[ 1 ] );
	DT.insert( pt );
	// VertexHandle vh = DT.insert( pt );
	// pt2vh[ pixel.xy ] = vh;
      }
      // Build triangulated surface
      _I.resize( _subset.size() );
      VertexIndex i = 0;
      std::map< Point, VertexIndex > pt2vtx;
      for ( auto pixel : _subset ) {
	_T.addVertex( pixel.xy );
	_I[ i ] = pixel.v;
	pt2vtx[ pixel.xy ] = i;
	i++;
      }
      for ( FiniteFacesIterator it = DT.finite_faces_begin(),
	      itend = DT.finite_faces_end(); it != itend; ++it ) {
	FaceHandle fh = it;
	CGALPoint2 a = fh->vertex( 0 )->point();
	CGALPoint2 b = fh->vertex( 1 )->point();
	CGALPoint2 c = fh->vertex( 2 )->point();
	Point pa( (int) a[ 0 ], (int) a[ 1 ] );
	Point pb( (int) b[ 0 ], (int) b[ 1 ] );
	Point pc( (int) c[ 0 ], (int) c[ 1 ] );
	_T.addTriangle( pt2vtx[ pa ], pt2vtx[ pc ], pt2vtx[ pb ] );
      }
      bool ok = _T.build();
      trace.info() << "Build triangulation: "
		   << ( ok ? "OK" : "ERROR" ) << std::endl;
      
    }
    
    // ------------------------- Private Datas --------------------------------
  private:

    
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ImageTriangulation

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ImageTriangulation_h

#undef ImageTriangulation_RECURSES
#endif // else defined(ImageTriangulation_RECURSES)
