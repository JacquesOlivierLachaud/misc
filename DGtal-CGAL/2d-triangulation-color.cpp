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
#include <DGtal/shapes/Shapes.h>
#include <DGtal/shapes/ShapeFactory.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/geometry/curves/GridCurve.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/GradientColorMap.h>
#include <DGtal/io/colormaps/GrayscaleColorMap.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/PPMWriter.h"
#include <DGtal/geometry/curves/ArithmeticalDSS.h>
#include <DGtal/geometry/curves/SaturatedSegmentation.h>
#include "DGtal/io/readers/PointListReader.h"

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

#include "cairo.h"

#include "Triangulation2DHelper.h"
#include "UmbrellaPart2D.h"
#include "Auxiliary.h"

static const double EPSILON = 0.0000001;
template <typename CGALPoint>
DGtal::Z2i::Point toDGtal( const CGALPoint & p )
{
  return DGtal::Z2i::Point ( p.x(),
                             p.y() );
}
template <typename CGALPoint>
CGALPoint toCGAL( const DGtal::Z2i::Point & p )
{
  return CGALPoint( p[ 0 ], p[ 1 ] );
}

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

template <typename Image, typename ToScalarFunctor>
class Derivatives
{
public:
  typedef typename Image::Domain Domain;
  typedef typename Domain::Point Point;
  const Image & myImage;
  const ToScalarFunctor & myFunctor;

  struct FFunctor {
    const Derivatives & _der;
    FFunctor( const Derivatives & der ) : _der( der ) {}
    inline double operator()( double x, double y ) const
    { return _der.f( x, y ); }
  };
  struct FxFunctor {
    const Derivatives & _der;
    FxFunctor( const Derivatives & der ) : _der( der ) {}
    inline double operator()( double x, double y ) const
    { return _der.dx( x, y ); }
  };
  struct FyFunctor {
    const Derivatives & _der;
    FyFunctor( const Derivatives & der ) : _der( der ) {}
    inline double operator()( double x, double y ) const
    { return _der.dy( x, y ); }
  };

  Derivatives( const Image & image, const ToScalarFunctor & f )
    : myImage( image ), myFunctor( f )
  {}
  Derivatives( const Derivatives & other )
    : myImage( other.myImage ), myFunctor( other.myFunctor )
  {}

  typedef double (Derivatives::*PointFunctorMethod)( Point ) const;

  double interpolate( double x, double y, 
		      PointFunctorMethod fct ) const
  {
    Point p( (int) floor( x ), (int) floor( y ) );
    double lx = x - p[ 0 ];
    double ly = y - p[ 1 ];
    return (1.0-lx) * (1.0-ly) * (this->*fct)( p ) 
      +     lx      * (1.0-ly) * (this->*fct)( p + Point(1,0) ) 
      +    (1.0-lx) * ly       * (this->*fct)( p + Point(0,1) ) 
      +     lx      * ly       * (this->*fct)( p + Point(1,1) );
  }

  double f( Point p ) const
  {
    if ( myImage.domain().isInside( p ) ) 
      return myFunctor( myImage( p ) );
    Point q = p;
    Point lo = myImage.domain().lowerBound();
    Point up = myImage.domain().upperBound();
    for ( unsigned int i = 0; i < 2; ++i ) {
      if ( p[ i ] < lo[ i ] ) {
	q[ i ] = lo[ i ];
	p[ i ] = 2*q[ i ] - p[ i ];
      }
      else if ( p[ i ] > up[ i ] ) {
	q[ i ] = up[ i ];
	p[ i ] = 2*q[ i ] - p[ i ];
      }
    }
    return 2 * myFunctor( myImage( q ) ) - myFunctor( myImage( p ) );
  }

  double f( double x, double y ) const
  {
    return interpolate( x, y, &Derivatives::f );
  }
  
  double dx( Point p ) const
  {
    double v1  = f( p + Point( 1, 0 ) );
    double vm1 = f( p - Point( 1, 0 ) );
    return (v1-vm1)/2.0;
  }

  double dx( double x, double y ) const
  {
    return interpolate( x, y, &Derivatives::dx );
  }

  double dy( Point p ) const
  {
    double v1  = f( p + Point( 0, 1 ) );
    double vm1 = f( p - Point( 0, 1 ) );
    return (v1-vm1)/2.0;
  }

  double dy( double x, double y ) const
  {
    return interpolate( x, y, &Derivatives::dy );
  }

  double dx2( Point p ) const
  {
    double v1 = f( p + Point( 1, 0 ) );
    double v0 = f( p );
    double vm1 = f( p - Point( 1, 0 ) );
    return (v1+vm1-2.0*v0)/2.0;
  }
  double dy2( Point p ) const
  {
    double v1 = f( p + Point( 1, 0 ) );
    double v0 = f( p );
    double vm1 = f( p - Point( 1, 0 ) );
    return (v1+vm1-2.0*v0)/2.0;
  }

  
};


/**
   An AVT is an abbreviation of affine valued triangulation. It is a
   triangulation T whose vertices are valued (such that the values are
   comparable. The remarkable property is that, for any value v, the
   subtriangulation (T >= v) is a triangulation of the convex hull of
   (T >= v) relative to R^2 \setminus (T < v).

   @tparam TTriangulation2 any kind of CGAL 2D triangulation.

   @tparam TKernel2 any kind of CGAL 2D Kernel, for instance
   CGAL::Exact_predicates_inexact_constructions_kernel
*/
template <typename TTriangulation2, typename TKernel2, typename TValue>
class AVT
{
public:
  typedef TTriangulation2                                   Triangulation2;
  typedef TKernel2                                          Kernel2;
  typedef TValue                                            Value;
  typedef AVT<Triangulation2, Kernel2, Value>               Self;
  typedef typename Triangulation2::Finite_vertices_iterator FiniteVerticesIterator;
  typedef typename Triangulation2::Vertex_circulator        VertexCirculator;
  typedef typename Triangulation2::Vertex_handle            VertexHandle;
  typedef typename Triangulation2::Edge_iterator            EdgeIterator;
  typedef typename Triangulation2::Finite_edges_iterator    FiniteEdgesIterator;
  typedef typename Triangulation2::All_edges_iterator       AllEdgesIterator;
  typedef typename Triangulation2::Edge                     Edge;
  typedef typename Triangulation2::Finite_faces_iterator    FiniteFacesIterator;
  typedef typename Triangulation2::Face_handle              FaceHandle;
  typedef typename Triangulation2::Point                    Point;
  typedef typename Triangulation2::Point                    Point2;
  typedef typename Kernel2::Vector_2                        Vector;
  typedef typename Kernel2::Vector_2                        Vector2;
  typedef typename Kernel2::FT                              Coordinate;
  typedef typename Kernel2::FT                              Component;
  typedef typename std::map<VertexHandle,Value>             VertexHandle2ValueMap;
  typedef typename std::map< std::pair<VertexHandle,VertexHandle>,
                             Value>                         Edge2ValueMap;
  typedef typename std::map<FaceHandle,Value>               FaceHandle2ValueMap;
  typedef DGtal::UmbrellaPart2D<Triangulation2,Kernel2>     Strip;
  typedef typename std::set<Edge>                           EdgeSet;
  typedef typename std::map<Point,VertexHandle>             Point2VertexHandleMap;

  /**
     This structure is used as a predicate for vertex/edge/face, such
     that this predicate returns 'true' whenever the value stored at
     the vertex/edge/face is no greater than \a myGreatestValue, given
     at instanciation.
  */
  struct NoGreaterThanValuePredicate {
    Self* myAVT;
    Value myGreatestValue;

    NoGreaterThanValuePredicate() : myAVT( 0 ) {}

    NoGreaterThanValuePredicate( Self & avt, const Value & v )
      : myAVT( &avt ), myGreatestValue( v ) {}

    inline
    bool operator()( VertexHandle v ) const
    {
      return myAVT->value( v ) <= myGreatestValue;
    }

    inline
    bool operator()( const Edge & e ) const
    {
      return myAVT->value( e ) <= myGreatestValue;
    }

    inline
    bool operator()( const FaceHandle & f ) const
    {
      return myAVT->value( f ) <= myGreatestValue;
    }
  };

  
protected:
  /// The current triangulation.
  Triangulation2 _T;
  const Value _invalid;
  /// A mapping VertexHandle -> Value that stores for each vertex its value.
  VertexHandle2ValueMap _vFct;
  /// A mapping Edge -> Value that stores for each edge its value (nb:
  /// any edge has the same value as its mirror edge).
  Edge2ValueMap _eFct;
  /// A mapping FaceHandle -> Value that stores for each face its value.
  FaceHandle2ValueMap _fFct;
  /// A mapping Point -> VertexHandle that stores for each point its corresponding vertex.
  Point2VertexHandleMap _p2vhMap
;
  DGtal::Board2D _debugBoard;

public:
  /// Provides useful methods on triangulation
  DGtal::Triangulation2DHelper<Triangulation2,Kernel2> TH;

  //------------------------------ basic services ---------------------------------------
public:

  /// Constructor.
  AVT( Value invalid = Value() )
    : TH( _T ), _invalid( invalid ) 
  {
    setInfiniteValue( invalid );
  }

  /// The object is reseted.
  void clear()
  {
    _T.clear();
    _vFct.clear();
    _eFct.clear();
    _fFct.clear();
    _p2vhMap.clear();
    setInfiniteValue( invalid );
  }

  /**
     @return a const reference to the current triangulation.
  */
  inline const Triangulation2 & T() const
  { return _T; }

  /**
     Adds the point \a pt with value \a val to the triangulation.
  */
  inline void add( const Point & pt, Value val )
  {
    VertexHandle vh = _T.insert( pt );
    _p2vhMap[ pt ] = vh;
    setValue( vh, val );
  }

  /** 
      Adds a sequence of (point,value) to the triangulation. Here, a
      range [itb,ite) is given for the points and a start iterator
      itvalb for the values.
      
      @param itb an iterator on the first point to add.
      @param ite an iterator after the last point to add.
      @param itvalv an iterator on the first value to add.

      The value type of PointIterator should be Point.
      The value type of ValueIterator should be Value.
  */
  template <typename PointIterator, typename ValueIterator>
  void add( PointIterator itb, PointIterator ite, ValueIterator itvalb )
  {
    for ( ; itb != ite; ++itb, ++itvalb )
      add( *itb, *itvalb );
  }

  template <typename PointValueIterator>
  void add( PointValueIterator itb, PointValueIterator ite )
  {
    for ( ; itb != ite; ++itb )
      add( itb->first, itb->second );
  }

  bool isQuadSaddle( const Edge & e ) const
  {
    VertexHandle v[ 4 ];
    Value val[ 4 ];
    v[ 0 ] = TH.target( e );
    v[ 1 ] = e.first->vertex( e.second ); 
    v[ 2 ] = TH.source( e );
    Edge me = T().mirror_edge( e );
    v[ 3 ] = me.first->vertex( me.second ); 
    for ( int i = 0; i < 4; ++i ) 
      val[ i ] = value( v[ i ] );
    return ( std::min( val[ 0 ], val[ 2 ] ) > std::max( val[ 1 ], val[ 3 ] ) )
      || ( std::max( val[ 0 ], val[ 2 ] ) < std::min( val[ 1 ], val[ 3 ] ) );
  }

  std::pair<Point,Value> middleQuad( const Edge & e ) const
  {
    VertexHandle v[ 4 ];
    v[ 0 ] = TH.target( e );
    v[ 1 ] = e.first->vertex( e.second ); 
    v[ 2 ] = TH.source( e );
    Edge me = T().mirror_edge( e );
    v[ 3 ] = me.first->vertex( me.second ); 
    Value mval = 0;
    double x = 0.0;
    double y = 0.0;
    for ( int i = 0; i < 4; ++i ) 
      {
	x += v[ i ]->point().x();
	y += v[ i ]->point().y();
	mval = mval + value( v[ i ] );
      }
    return std::make_pair( Point( x / 4.0, y / 4.0 ), mval / 4 );
  }

  void processSaddlePoints()
  {
    Edge e;
    std::vector< std::pair<VertexHandle, VertexHandle> > saddles;
    for ( FiniteEdgesIterator it = T().finite_edges_begin(), itend = T().finite_edges_end();
          it != itend; ++it )
      {
        e = *it; 
	if ( isQuadSaddle( e ) )
	  saddles.push_back( std::make_pair( TH.source( e ), TH.target( e ) ) );
      }
    for ( typename std::vector< std::pair<VertexHandle, VertexHandle> >::const_iterator 
	    it = saddles.begin(), itend = saddles.end();
	  it != itend; ++it )
      {
	VertexHandle v0 = it->first;
	VertexHandle v1 = it->second;
	if ( ! T().is_edge( v0, v1, e.first, e.second ) ) continue;
	if ( ! isQuadSaddle( e ) )                        continue;
	if ( TH.source( e ) != v0 )                       e = T().mirror_edge( e );
	std::pair<Point,Value> res = middleQuad( e );
	add( res.first, res.second );
      }
  }

  /**
     @param pt any point of the triangulation.
     @return the corresponding vertex.
  */
  VertexHandle vertex( const Point & pt ) const
  {
    typename Point2VertexHandleMap::const_iterator it = _p2vhMap.find( pt );
    ASSERT( it != _p2vhMap.end() );
    return it->second;
  }

  /**
     Used by fullRelativeHull to insert candidate concavities into the Queue.
     
     @param e any finite edge of the triangulation.
  */
  void insertQueue( std::set< std::pair< VertexHandle, VertexHandle > > & inQueue,
		    const Edge & e ) 
  {
    ASSERT( ! T().is_infinite( e ) );
    VertexHandle v0 = TH.source( e );
    VertexHandle v1 = TH.target( e );
    const Value val_v0 = value( v0 );
    const Value val_v1 = value( v1 );
    const Value val_e =  value( e );
    if ( std::min( val_v0, val_v1 ) > val_e )
      {
        DGtal::trace.warning() << "[AVT::insertQueue] Weird low value for edge:"
                               << " e=" << val_e << " v0=" << val_v0 << " v1=" << val_v1
                               << std::endl;
        return;
      }
    if ( std::max( val_v0, val_v1 ) < val_e )
      {
        DGtal::trace.warning() << "[AVT::insertQueue] Weird big value for edge:"
                               << " e=" << val_e << " v0=" << val_v0 << " v1=" << val_v1
                               << std::endl;
        return;
      }
    if ( ( val_v1 < val_v0 ) && ( val_e < val_v0 ) )
      {
	Strip strip( T() );
        NoGreaterThanValuePredicate predNoGreaterThanV1( *this, val_e );
        strip.init( predNoGreaterThanV1, predNoGreaterThanV1, e );
        if ( strip.isConcave() ) 
          inQueue.insert( std::make_pair( v0, v1 ) );
      }
    else if ( ( val_v0 < val_v1 ) && ( val_e < val_v1 ) )
      {
	Strip strip( T() );
        NoGreaterThanValuePredicate predNoGreaterThanV1( *this, val_e );
        strip.init( predNoGreaterThanV1, predNoGreaterThanV1, T().mirror_edge( e ) );
        if ( strip.isConcave() ) 
          inQueue.insert( std::make_pair( v1, v0 ) );
      }
  }

  /**
     This procedure computes the relative hull by locating and
     flipping concavities in the triangulation. Concavities are
     defined as lower part of an umbrella. Complexity is O(e), where e
     is the number of edges.

     This is a kind of simplest algorithm, where we do not try to sort
     concavities, or push again edges around
     concavities. Surprisingly, it is hard to do better without
     complexifying a lot the way concavities (v0,v1) are pushed into
     the queue.

     @return 'true' if some concavity was flipped, 'false' when no concavity was found.
  */
  bool fullRelativeHull()
  {
    bool changes = false;
    // static const Label rmark = 1;
    std::set< std::pair< VertexHandle, VertexHandle > > inQueue;
    Strip strip( T() );
    // DGtal::trace.beginBlock( "Searching concavities" );
    for ( FiniteEdgesIterator it = T().finite_edges_begin(), itend = T().finite_edges_end();
          it != itend; ++it )
      {
        Edge e = *it;
	insertQueue( inQueue, e );
      }
    DGtal::trace.info() << "- Found " << inQueue.size() << " potential concavities." 
			<< std::endl;
    // DGtal::trace.endBlock();
    // DGtal::trace.beginBlock( "Flipping concavities" );
    unsigned int nb_checked = 0;
    unsigned int nb_flipped = 0;
    unsigned int nb_concave = 0;
    Edge e;
    while ( ! inQueue.empty() )
      {
	VertexHandle v0 = inQueue.begin()->first;
	VertexHandle v1 = inQueue.begin()->second;
        ++nb_checked;
        if ( nb_checked % 1000 == 0 ) 
          DGtal::trace.info() << "- Queue=" << inQueue.size()
                              << ", flipped " << nb_flipped << "/" << nb_checked 
                              << " edges in a concavity." << std::endl;
	if ( ! T().is_edge( v0, v1, e.first, e.second ) )
	  {
	    inQueue.erase( inQueue.begin() );
	    continue; // (v0,v1) is not an edge anymore.
	  }
	if ( TH.source( e ) != v0 ) 
	  e = T().mirror_edge( e );
        Value val_V1 = value( e ); // value( TH.target( e ) );
        if ( val_V1 >= value( v0 ) ) // edge has already the value of the source.
	  {
	    inQueue.erase( inQueue.begin() );
	    continue; 
	  }
        NoGreaterThanValuePredicate predNoGreaterThanV1( *this, val_V1 );
	strip.init( predNoGreaterThanV1, predNoGreaterThanV1, e );
	ASSERT( strip.isValid() );
        if ( strip.isConcave() ) 
          {
            // Extract correct level/value for the concavity.
            Value val0 = std::max( value( strip.e0() ), // first edge
                                   value( T().mirror_edge( strip.e0() ).first ) ); // face just before.
            Value val1 = std::max( value( strip.en() ), // last edge
                                   value( strip.en().first ) ); // face just after.
            Value val = std::min( std::min( val0, val1 ),
                                  value( TH.source( e ) ) );
            ASSERT( val >= val_V1 );
            // We may already update the value of border edges.
            setValue( strip.e0(), std::max( value( strip.e0() ), val ) );
            setValue( strip.en(), std::max( value( strip.en() ), val ) );
	    unsigned int i = strip.getFlippableEdgeIndex();
	    if ( i != strip.size() ) // at least one edge is flippable.
	      {
                Edge fedge = strip.e( i );
                Edge fmedge = T().mirror_edge( fedge );
		Edge edge_quad = T().mirror_edge( strip.e( i-1 ) );
                FaceHandle f1 = fedge.first;
                FaceHandle f2 = fmedge.first;
                // remove value for f1, f2, fedge, fmedge
                eraseValue( f1 ); 
                eraseValue( f2 );
                eraseValue( fedge );
		Edge mirror_first = T().mirror_edge( strip.e( 0 ) );
                // DGtal::trace.info() << "  - Flipping " << TH.source( fedge )->point() 
                //                     << " -> " <<  TH.target( fedge )->point() << std::endl;
		_T.flip( fedge.first, fedge.second );
                fedge = T().mirror_edge( edge_quad );
                Edge new_edge = TH.nextCCWAroundFace( fedge );
                // Make sure all values are reseted.
                eraseValue( new_edge ); 
                eraseValue( new_edge.first );
                eraseValue( T().mirror_edge( new_edge ).first );
                Edge fedge2 = TH.nextCCWAroundFace( new_edge );
                // setValue( fedge, std::max( value( fedge ), val_V1+1 ) );
                // setValue( new_edge, std::max( value( new_edge ), val_V1+1 ) );
                // setValue( fedge2, std::max( value( fedge2 ), val_V1+1 ) );
		if ( strip.size() == 3 ) 
		  { // last flip made a triangle
                    // setValue( new_edge, std::max( value( new_edge ), val ) );
                    setValue( new_edge.first, std::max( value( new_edge.first ), val ) );
		    inQueue.erase( inQueue.begin() );
		  }
		changes = true; ++nb_flipped;
	      }
	    else // none are flippable. This is a concave/flat piece.
	      { // All faces and edges are set to value V1+1
                setValue( strip.f( 0 ), std::max( value( strip.f( 0 ) ), val ) );
                for ( unsigned int i = 1; i < strip.size() - 1; ++i )
                  {
                    setValue( strip.e( i ), std::max( value( strip.e( i ) ), val ) );
                    setValue( strip.f( i ), std::max( value( strip.f( i ) ), val ) );
                  }
                changes = true; ++nb_concave;
		inQueue.erase( inQueue.begin() );
	      }
          }
	else
	  inQueue.erase( inQueue.begin() );
      }
    ASSERT( inQueue.empty() );
    DGtal::trace.info() << "- Flipped " << nb_flipped 
                        << ", concave " << nb_concave
                        << " / " << nb_checked << " edges in a concavity." << std::endl;
    //if ( nb_flipped == 1 && nb_checked == 2 ) changes = false;
    // DGtal::trace.endBlock();
    return changes;
  }
  
  // ---------------------- services --------------------------------
public:
  
  /**
     Given a vertex specified by \a v and a value \val, computes the
     CCW oriented local contour for this value and returns 'true' if
     there is an isocontour of such value at this place. The pairs of
     consecutive edges are returned in \a in_edges and \a
     out_edges. The edges of the form (v1,v) are in \a in_edges, the
     edges of the form (v,v2) are in \a out_edges. Since the
     isocontour is oriented CCW, each cell of any edge has value
     greater or equal to \a val. There should be as many \a in_edges
     as \a out_edges.
     
     @param[out] in_edges the edges forming the isocontour around vertex \a v.
     @param[out] out_edges the edges forming the isocontour around vertex \a v.
     @param[in] v the vertex handle specifying where we look for the isocontour.
     @param[in] val the value of the isocontour (there may be several isocontours at this vertex).
  */

  bool localIsocontour( std::vector<Edge> & in_edges, 
			std::vector<Edge> & out_edges, 
			VertexHandle v, Value val )
  {
    bool first_edge_is_out = false;
    in_edges.clear();
    out_edges.clear();
    Edge first = *( T().incident_edges( v ) );
    if ( TH.source( first ) != v ) first = T().mirror_edge( first );
    Edge e = first;
    do {
      ASSERT( TH.source( e ) == v );
      Value val1 = value( e.first );
      Value val2 = value( T().mirror_edge( e ).first );
      if ( ( val1 >= val ) && ( val > val2 ) )
	{
	  if ( ! first_edge_is_out )
	    {
	      first_edge_is_out = true;
	      first = e; // starts at this position.
	    }
	  out_edges.push_back( e );
	}
      else if ( ( val1 < val ) && ( val <= val2 ) )
	{
	  if ( first_edge_is_out )
	    in_edges.push_back( T().mirror_edge( e ) );
	}
      e = TH.nextCCWAroundSourceVertex( e ); 
    } while ( e != first );
    if ( in_edges.size() != out_edges.size() )
      DGtal::trace.warning() << "[AVT::localIsocontour] point=" << v->point() 
			     << " #in=" << in_edges.size()
			     << " #out=" << out_edges.size()
			     << std::endl;
    return ( in_edges.size() != 0 ) && ( out_edges.size() == in_edges.size() );
  }

  static Component radius( Component L0, Component L1, Component L01 )
  {
    if ( ( L1 + L01 <= L0 ) || ( L0 + L01 <= L1 ) || ( L0 + L1 <= L01 ) ) 
      return DBL_MAX;
    if ( ( L0 <= 0.0 )  || ( L1 <= 0.0 ) || ( L01 <= 0.0 ) ) return 0.0;
    return ( L0 * L1 * L01 ) 
      / sqrt( ( L0 + L1 + L01 ) * ( -L0 + L1 + L01 ) * ( L0 - L1 + L01 ) * ( L0 + L1 - L01 ) );
  }

  /**
     @return 'true' iff the two edges form a corner.
     @pre e1 and e2 should be consecutive edges, i.e. target( e1 ) = source( e2 )
  */
  bool isSmooth( const Edge & e1, const Edge & e2, Component dh ) const
  {
    ASSERT( TH.target( e1 ) == TH.source( e2 ) );
    Point v0 = TH.source( e1 )->point();
    Point v  = TH.target( e1 )->point();
    Point v1 = TH.target( e2 )->point();
    Component LL0 = ( v - v0 ).squared_length();
    // bool aligned0 = ( ( v - v0 ).x() == 0 ) || ( ( v - v0 ).y() == 0 );
    Component LL1 = ( v1 - v ).squared_length();
    // bool aligned0 = ( ( v1 - v ).x() == 0 ) || ( ( v1 - v ).y() == 0 );
    Component LL01 = ( v1 - v0 ).squared_length();
    Component L0 = sqrt( LL0 ), L1 = sqrt( LL1 ), L01 = sqrt( LL01 );
    Component R = radius( L0, L1, L01 );
    Component eps1 = R - sqrt( R*R - LL0 / 4.0 );
    Component eps2 = R - sqrt( R*R - LL1 / 4.0 );
    return ( eps1 <= sqrt(2.0) * dh ) && ( eps2 <= sqrt(2.0) * dh );

    // Component R0_min = ( LL0 / (4.0*dh) + dh ) / ( 2.0 * sqrt( 2.0 ) );
    // Component R0_max = ( LL0 / (4.0*dh) + dh ) * sqrt( 2.0 ) / 2.0;
    // Component R1_min = ( LL1 / (4.0*dh) + dh ) / ( 2.0 * sqrt( 2.0 ) );
    // Component R1_max = ( LL1 / (4.0*dh) + dh ) * sqrt( 2.0 ) / 2.0;

    
    // Component R = radius( L0, L1, L01 ); 
    // Component Rmin = std::min( R, radius( L0 - dh, L1, L01 ) );
    // Rmin = std::min( Rmin, radius( L0, L1 - dh, L01 ) );
    // Rmin = std::min( Rmin, radius( L0, L1, L01 - dh ) );
    // Component Rmax = std::max( R, radius( L0 + dh, L1, L01 ) );
    // Rmax = std::max( Rmin, radius( L0, L1 + dh, L01 ) );
    // Rmax = std::max( Rmin, radius( L0, L1, L01 + dh ) );

    // return ( ( R0_min <= Rmax ) && ( Rmin <= R0_max ) )
    //   || ( ( R1_min <= Rmax ) && ( Rmin <= R1_max ) );
  }
  
  // ---------------------- value services --------------------------------
public:

  /**
     @return the "invalid" value, used for instance for the infinite vertex.
  */
  inline
  const Value invalid() const
  { return _invalid; }

  /**
     @param v any handle to a vertex of this triangulation.
     @return the value at vertex \a v, or invalid() if not defined.
  */
  Value value( VertexHandle v ) const
  {
    typename VertexHandle2ValueMap::const_iterator it = _vFct.find( v );
    if ( it == _vFct.end() ) return _invalid;
    return it->second;
  }

  /**
     @param e any edge of this triangulation.

     @return the value at edge \a e, which is computed as the minimum
     of the value at the extremeties if it is not yet defined.
  */
  Value value( Edge e )
  {
    // Pick "smallest" edge to find value.
    Edge mirror_e = T().mirror_edge( e );
    if ( TH.source( mirror_e ) < TH.source( e ) ) e = mirror_e;
    std::pair<VertexHandle, VertexHandle> edge = std::make_pair( TH.source( e ), TH.target( e ) );
    typename Edge2ValueMap::const_iterator it = _eFct.find( edge );
    if ( it == _eFct.end() ) 
      {
        Value v = std::min( value( TH.source( e ) ), 
                            value( TH.target( e ) ) );
        _eFct[ edge ] = v;
        return v;
      }
    return it->second;
  }

  /**
     @param f any handle to a face of this triangulation.

     @return the value at face \a f, which is computed as the minimum
     of the value at its edges if it is not yet defined.
  */
  Value value( FaceHandle f )
  {
    typename FaceHandle2ValueMap::const_iterator it = _fFct.find( f );
    if ( it == _fFct.end() ) 
      {
        Value v = std::min( std::min( value( Edge( f, 0 ) ), value( Edge( f, 1 ) ) ), 
                            value( Edge( f, 2 ) ) );
        _fFct[ f ] = v;
        return v;
      }
    return it->second;
  }

  /**
     Sets the value for the given vertex.
     @param vh any vertex handle.
     @param val any value.
  */
  inline void setValue( VertexHandle vh, Value val )
  {
    _vFct[ vh ] = val;
  }

  /**
     Sets the value for the given edge
     @param e any edge.
     @param val any value.
  */
  inline void setValue( Edge e, Value val )
  {
    Edge mirror_e = T().mirror_edge( e );
    if ( TH.source( mirror_e ) < TH.source( e ) ) e = mirror_e;
    std::pair<VertexHandle, VertexHandle> edge = std::make_pair( TH.source( e ), TH.target( e ) );
    _eFct[ edge ] = val;
  }

  /**
     Sets the value for the given face
     @param fh any face handle.
     @param val any value.
  */
  inline void setValue( FaceHandle fh, Value val )
  {
    _fFct[ fh ] = val;
  }

  /**
     Sets the value for the infinite vertex.
  */
  inline void setInfiniteValue( Value v )
  {
    _vFct[ _T.infinite_vertex() ] = v;
  }

  /**
     Erases the value associated with the vertex handle vh.
  */
  void eraseValue( VertexHandle vh )
  {
    _vFct.erase( vh );
  }

  /**
     Erases the value associated with the edge e.
  */
  void eraseValue( Edge e )
  {
    // Pick "smallest" edge to find value.
    Edge mirror_e = T().mirror_edge( e );
    if ( TH.source( mirror_e ) < TH.source( e ) ) e = mirror_e;
    std::pair<VertexHandle, VertexHandle> edge = std::make_pair( TH.source( e ), TH.target( e ) );
    _eFct.erase( edge );
  }

  /**
     Erases the value associated with the face handle fh.
  */
  void eraseValue( FaceHandle fh )
  {
    _fFct.erase( fh );
  }


  // ---------------------- gradient services --------------------------------
public:

  /**
     Computes the gradient within a triangular face (p_0,p_1,p_2) with
     the standard formula:

     Grad u = (1/ (2*Af)) * ( \sum_i=0^2 u_i N x e_i ) 

     where u is the fonction, u_i the value at point p_i , e_i =
     p_{i+1} - p_i, Af is the area of the triangle, N its normal
     vector.

     @param fh a handle to a face of the triangulation.
     @return the discrete gradient at \a fh.

     @note For computations, the triangles belongs to the xy-plane and
     u_i is the value in the AVT.
  */
  Vector2 gradient( FaceHandle fh ) const
  {
    ASSERT( ! T().is_infinite( fh ) );
    Point p0 = fh->vertex( 0 )->point();
    Point p1 = fh->vertex( 1 )->point();
    Point p2 = fh->vertex( 2 )->point();
    double f0 = (double) value( fh->vertex( 0 ) );
    double f1 = (double) value( fh->vertex( 1 ) );
    double f2 = (double) value( fh->vertex( 2 ) );
    double Af = twiceAreaTriangle( p0, p1, p2 );
    Vector2 e2 = p1 - p0;
    Vector2 e0 = p2 - p1;
    Vector2 e1 = p0 - p2;
    return Vector2( -(f0*e0.y() + f1*e1.y() + f2*e2.y()) / Af,
		     (f0*e0.x() + f1*e1.x() + f2*e2.x()) / Af );
  }

  Vector2 gradient( VertexHandle vh ) const
  {
    Edge start = *( T().incident_edges( vh ) );
    if ( TH.source( start ) != vh ) start = T().mirror_edge( start );
    Edge e = start;
    Vector2 grad( 0.0, 0.0 );
    unsigned int nb = 0;
    do {
      if ( ! T().is_infinite( e.first ) ) {
	grad = grad + gradient( e.first );
	++nb;
      }
      e = TH.nextCCWAroundSourceVertex( e );
    } while ( e != start );    
    return grad / nb;
  }

  Vector2 gradient( Edge e ) const
  {
    Vector2 grad( 0.0, 0.0 );
    unsigned int nb = 0;
    if ( ! T().is_infinite( e.first ) ) {
      grad = grad + gradient( e.first );
      ++nb;
    }
    e = T().mirror_edge( e );
    if ( ! T().is_infinite( e.first ) ) {
      grad = grad + gradient( e.first );
      ++nb;
    }
    return grad / nb;
  }

  Vector2 gradient( double x, double y, FaceHandle hint = FaceHandle() ) const
  {
    typename Triangulation2::Locate_type lt;
    int li;
    FaceHandle	fh = T().locate( Point2( x,y ), lt, li, hint );
    Vector2 result( 0.0, 0.0 );
    switch ( lt ) {
    case Triangulation2::VERTEX:
      result = gradient( fh->vertex( li ) ); break;
    case Triangulation2::EDGE:
      result = gradient( Edge( fh, li ) ); break;
    case Triangulation2::FACE:
      result = gradient( fh ); break;
    }
    return result;
  }

  Vector2 gradientPixel( double x0, double y0, double s, FaceHandle hint = FaceHandle() ) const
  {
    const int n = 8;
    const double N = (double) n;
    const double pas = 1.0 / N;
    const double shift = pas / 2.0;
    hint = T().locate( Point2( x0, y0 ), hint );
    double spas = s * pas;
    double ymin = y0 - s * shift * (N+1.0); 
    double xmin = x0 - s * shift * (N+1.0);
    double y = ymin;
    Vector2 grad( 0.0, 0.0 );
    for ( int i = 0; i < n; ++i, y += spas )
      {
        double x = xmin;
        for ( int j = 0; j < n; ++j, x += spas )
          {
            grad = grad + gradient( x, y, hint );
          }
      }
    return grad * ( s*s * pas * pas );
  }

  // ---------------------- various services --------------------------------
public:
  
  static double linearInterpolation( Point p, Point p1, Point p2, 
				     double v1, double v2 )
  {
    double ux = p2.x() - p1.x();
    double uy = p2.y() - p1.y();
    double mu = ( abs( ux ) >= abs( uy ) )
      ? ( p.x() - p1.x() ) / ux
      : ( p.y() - p1.y() ) / uy;
    return ( 1.0 - mu ) * v1 + mu * v2;
  }

  static double twiceAreaTriangle( Point p1, Point p2, Point p3 )
  {
    double ux = p2.x() - p1.x();
    double uy = p2.y() - p1.y();
    double vx = p3.x() - p1.x();
    double vy = p3.y() - p1.y();
    return abs( ux * vy - uy * vx );
  }

  static double midValue( double v1, double v2, double v3 )
  {
    return ( v1 + v2 + v3 ) / 3.0;
  }
  static double midValue( double v1, double v2 )
  {
    return ( v1 + v2 ) / 2.0;
  }
  static Point midPoint( Point p1, Point p2, Point p3 )
  {
    return Point( ( p1.x() + p2.x() + p3.x() ) / 3.0, ( p1.y() + p2.y() + p3.y() ) / 3.0 );
  }
  static Point midPoint( Point p1, Point p2 )
  {
    return Point( ( p1.x() + p2.x() ) / 2.0, ( p1.y() + p2.y() ) / 2.0 );
  }

  static double linearInterpolation( Point p, Point p1, Point p2, Point p3, 
				     double v1, double v2, double v3 )
  {
    double area = twiceAreaTriangle( p1, p2, p3 );
    double a1 = twiceAreaTriangle( p, p2, p3 ) / area;
    double a2 = twiceAreaTriangle( p, p3, p1 ) / area;
    double a3 = twiceAreaTriangle( p, p1, p2 ) / area;
    return a1 * v1 + a2 * v2 + a3 * v3;
  }

  /**
     @param[in] p any point.
     @param[in] gouraud when 'true', performs Gouraud interpolation.
     @return the value at position \a p. 
  */
  double preciseValue( Point p, bool gouraud, FaceHandle hint = FaceHandle() )
  {
    typename Triangulation2::Locate_type lt;
    int li;
    FaceHandle	fh = T().locate( p, lt, li, hint );
    double result = 0.0;
    switch ( lt ) {
    case Triangulation2::VERTEX:
      { // If lt==VERTEX the variable li is set to the index of the vertex
	result = (double) value( fh->vertex( li ) );
      } break;
    case Triangulation2::EDGE:
      { // if lt==EDGE li is set to the index of the vertex opposite to the edge
	VertexHandle v1 = fh->vertex( (li+1)%3 );
	VertexHandle v2 = fh->vertex( (li+2)%3 );
	result = ( gouraud )
	  ? linearInterpolation( p, v1->point(), v2->point(), value( v1 ), value( v2 ) )
	  : (double) value( Edge( fh, li ) );
      } break;
    case Triangulation2::FACE:
      {
	VertexHandle v0 = fh->vertex( 0 );
	VertexHandle v1 = fh->vertex( 1 );
	VertexHandle v2 = fh->vertex( 2 );
	result = ( gouraud )
	  ? linearInterpolation( p, v0->point(), v1->point(), v2->point(), 
				 value( v0 ), value( v1 ), value( v2 ) )
	  : (double) value( fh );
      } break;
    case Triangulation2::OUTSIDE_CONVEX_HULL:
      std::cerr << "OUTSIDE_CONVEX_HULL" << std::endl;
      break;
    case Triangulation2::OUTSIDE_AFFINE_HULL:
      std::cerr << "OUTSIDE_AFFINE_HULL" << std::endl;
      break;
    }
    return result;
  }

  double errorTVOnTriangle( Point p1, Point p2, Point p3, 
			    double v1, double v2, double v3,
			    bool gouraud, int subdivision = 2 )
  {
    if ( subdivision <= 0 ) 
      {
	double val_one = midValue( v1, v2, v3 );
	double val_two = preciseValue( midPoint( p1, p2, p3 ), gouraud );
	return ( val_one - val_two ) * ( val_one - val_two ) 
	  * twiceAreaTriangle( p1, p2, p3 );
      }
    else
      return errorTVOnTriangle( p1, midPoint( p1, p2 ), midPoint( p1, p3 ),
				  v1, midValue( v1, v2 ), midValue( v1, v3 ),
				  gouraud, subdivision - 1 )
	+ errorTVOnTriangle( midPoint( p1, p2 ), p2, midPoint( p2, p3 ),
			     midValue( v1, v2 ), v2, midValue( v2, v3 ),
			     gouraud, subdivision - 1 )
	+ errorTVOnTriangle( midPoint( p1, p2 ), midPoint( p1, p3 ), midPoint( p2, p3 ),
			     midValue( v1, v2 ), midValue( v1, v3 ), midValue( v2, v3 ),
			     gouraud, subdivision - 1 )
	+ errorTVOnTriangle( midPoint( p1, p3 ), midPoint( p2, p3 ), p3,
			     midValue( v1, v3 ), midValue( v2, v3 ), v3,
			     gouraud, subdivision - 1 );
  }

  /**
     Estimates the resulting error in the reconstruction when removing vertex \a vh.

     @param[in] vh the vertex whose removal is checked.
  */
  double errorTVWhenRemoved( VertexHandle vh, bool gouraud, int sub = 2 )
  {
    Point p = vh->point();
    Value val = value( vh );
    std::vector<Point> neighbors;
    std::vector<Value> values;

    Self miniAVT( _invalid );
    typename Triangulation2::Vertex_circulator ci_start = T().incident_vertices( vh );
    typename Triangulation2::Vertex_circulator ci = ci_start;
    do {
      VertexHandle n = ci;
      neighbors.push_back( n->point() );
      values.push_back( value( n ) );
      miniAVT.add( neighbors.back(), values.back() );
      ++ci;
    } while ( ci != ci_start );
    double error = 0.0;
    for ( unsigned int i = 0; i < neighbors.size(); ++i )
      {
	unsigned int j = ( i+1 ) % neighbors.size();
	error += miniAVT.errorTVOnTriangle( p, neighbors[ i ], neighbors[ j ],
					    val, values[ i ], values[ j ], gouraud, sub );
      }
    return error;
  }

  struct VertexError {
    VertexHandle _vh;
    double _error;
    VertexError( VertexHandle vh, double error )
      : _vh( vh ), _error( error ) {}
  };

  struct VertexErrorComparator {
    bool operator()( const VertexError & ve1,
		     const VertexError & ve2 ) const
    {
      return ( ve1._error < ve2._error )
	|| ( ( ve1._error == ve2._error )
	     && ( ve1._vh < ve2._vh ) );
    }
  };

  bool isOnConvexHull( VertexHandle vh ) const
  {
    typename Triangulation2::Vertex_circulator ci_start 
      = T().incident_vertices( vh );
    typename Triangulation2::Vertex_circulator ci = ci_start;
    do {
      if ( T().is_infinite( ci ) ) return true;
      ++ci;
    } while ( ci != ci_start );    
    return false;
  }


  void compress( double ratio, bool gouraud, int sub )
  {
    // First constrain the whole triangulation
    DGtal::trace.info() << "[AVT::compress] Setting constraints." << std::endl;
    for ( FiniteEdgesIterator it = T().finite_edges_begin(), ite = T().finite_edges_end();
	  it != ite; ++it )
      {
	_T.insert_constraint( TH.source( *it ), TH.target( *it ) );
      }

    DGtal::trace.info() << "[AVT::compress] Removing vertices." << std::endl;
    std::set< VertexError, VertexErrorComparator > candidates;
    std::set< VertexHandle > modified_vertices;
    unsigned int nb = T().number_of_vertices();
    unsigned int to_remove = (unsigned int) ceil( ((double)nb) * ratio );
    unsigned int nb_candidates = 0;
    unsigned int nb_removed = 0;
    unsigned int nb_rescan = 0;
    double ratio_rescan = 0.25;
    while ( to_remove != 0 )
      {
	if ( candidates.empty() 
	     || ( nb_removed >= nb_rescan ) )
	  {
	    modified_vertices.clear();
	    candidates.clear();
	    nb_removed = 0;
	    nb_candidates = 0;
	    for ( FiniteVerticesIterator it = T().finite_vertices_begin(), 
		    ite = T().finite_vertices_end(); it != ite; ++it )
	      {
		if ( ! isOnConvexHull( it ) )
		  {
		    VertexError ve( it, errorTVWhenRemoved( it, gouraud, sub ) );
		    candidates.insert( ve );
		    ++nb_candidates;
		  }
	      }
	    nb_rescan = (unsigned int) ceil( ratio_rescan * (double) nb_candidates );
	    DGtal::trace.info() << "[AVT::compress]"
				<< " to_remove=" << to_remove
				<< " candidates=" << nb_candidates
				<< " rescan=" << nb_rescan << std::endl;
	  }
	VertexError ve = *( candidates.begin() );
	candidates.erase( candidates.begin() );
	// DGtal::trace.info() << "[" << to_remove << "]"
	// 		    << "Compressing pt=" << ve._vh->point() 
	// 		    << " error=" << ve._error << std::endl;
	typename std::set< VertexHandle >::iterator itm 
	  = modified_vertices.find( ve._vh );
	if ( itm != modified_vertices.end() )
	  { // this vertex was closed to a removed vertex.
	    continue;
	  }
	// Mark nearby vertices
	Edge start = *( T().incident_edges( ve._vh ) );
	if ( TH.source( start ) != ve._vh ) start = T().mirror_edge( start );
	Edge e = start;
	do {
	  modified_vertices.insert( TH.target( e ) );
	  eraseValue( e );
	  eraseValue( e.first );
	  e = TH.nextCCWAroundSourceVertex( e );
	} while ( e != start );
	// Remove vertex.
	eraseValue( ve._vh );
	modified_vertices.erase( ve._vh );
	_T.remove_incident_constraints( ve._vh );
	_T.remove( ve._vh );
	--to_remove;
	++nb_removed;
      }
    for ( FiniteVerticesIterator it = T().finite_vertices_begin(), 
	    ite = T().finite_vertices_end(); it != ite; ++it )
      _T.remove_incident_constraints( it );
  }
  
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



/**
   This class is intended for visualizing Affine Valued triangulation with CAIRO.
*/
template <typename TAVT, typename TValue>
class CairoViewerAVT
{
public:
  typedef TAVT                                 AVT;
  typedef TValue                               Value;
  typedef typename AVT::Triangulation2         Triangulation2;
  typedef typename AVT::VertexHandle           VertexHandle;
  typedef typename AVT::FiniteVerticesIterator FiniteVerticesIterator;
  typedef typename AVT::Edge                   Edge;
  typedef typename AVT::EdgeIterator           EdgeIterator;
  typedef typename AVT::FiniteEdgesIterator    FiniteEdgesIterator;
  typedef typename AVT::FiniteFacesIterator    FiniteFacesIterator;
  typedef typename AVT::FaceHandle             FaceHandle;
  typedef typename AVT::Point2                 Point2;
  typedef typename AVT::Vector2                Vector2;
  typedef DGtal::Board2D                       Board;
  typedef DGtal::Color                         Color;
  typedef DGtal::Z2i::Point                    PointZ2;

private:
  int _x0, _y0;
  int _width, _height;
  double _xf, _yf;
  bool _gouraud;
  cairo_surface_t* _surface;
  cairo_t* _cr;
  double _redf, _greenf, _bluef;

public:

  enum Mode { Gray, Red, Green, Blue };

  /**
     Constructor. Requires a board \a board for display and an affine
     valued triangulation \a avt.
  */
  CairoViewerAVT( int x0, int y0, int width, int height, 
                  double xfactor = 1.0, double yfactor = 1.0, bool gouraud = false )
    : _x0( x0 ), _y0( y0 ), _width( width ), _height( height ),
      _xf( xfactor ), _yf( yfactor ), _gouraud( gouraud )
  {
    _surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, width, height );
    _cr = cairo_create ( _surface );
    // Fill the background with black
    cairo_set_source_rgba( _cr, 0.0, 0.0, 0.0, 1.0 );
    cairo_rectangle ( _cr, 0, 0, _width, _height );
    cairo_fill( _cr );
  }

  ~CairoViewerAVT()
  {
    cairo_destroy( _cr );
    cairo_surface_destroy( _surface );
  }
  
  void save( const char* file_name ) const
  {
    cairo_surface_write_to_png( _surface, file_name );
  }

  void viewGouraudTriangle( Point2 a, Point2 b, Point2 c, 
			    double val_a, double val_b, double val_c ) 
  {
    cairo_pattern_t * pattern = cairo_pattern_create_mesh();
    /* Add a Gouraud-shaded triangle */
    cairo_mesh_pattern_begin_patch (pattern);
    cairo_mesh_pattern_move_to (pattern, i( a.x() ), j( a.y() ) );
    cairo_mesh_pattern_line_to (pattern, i( b.x() ), j( b.y() ) );
    cairo_mesh_pattern_line_to (pattern, i( c.x() ), j( c.y() ) );
    cairo_mesh_pattern_set_corner_color_rgb (pattern, 0, val_a * _redf, val_a * _greenf, val_a * _bluef );
    cairo_mesh_pattern_set_corner_color_rgb (pattern, 1, val_b * _redf, val_b * _greenf, val_b * _bluef );
    cairo_mesh_pattern_set_corner_color_rgb (pattern, 2, val_c * _redf, val_c * _greenf, val_c * _bluef );
    cairo_mesh_pattern_end_patch (pattern);
    cairo_set_source( _cr, pattern );
    cairo_move_to( _cr, i( a.x() ), j( a.y() ) );
    cairo_line_to( _cr, i( b.x() ), j( b.y() ) );
    cairo_line_to( _cr, i( c.x() ), j( c.y() ) );
    cairo_close_path( _cr );
    cairo_fill( _cr );
    cairo_pattern_destroy( pattern );
  }

  void viewFlatTriangle( Point2 a, Point2 b, Point2 c, double val )
  {
    cairo_set_source_rgb( _cr, val * _redf, val * _greenf, val * _bluef );
    cairo_move_to( _cr, i( a.x() ), j( a.y() ) );
    cairo_line_to( _cr, i( b.x() ), j( b.y() ) );
    cairo_line_to( _cr, i( c.x() ), j( c.y() ) );
    cairo_close_path( _cr );
    cairo_fill( _cr );
  }

  void viewAVTGouraudTriangle( AVT & avt, FaceHandle fh )
  {
    Point2 a = fh->vertex( 0 )->point();
    Point2 b = fh->vertex( 1 )->point();
    Point2 c = fh->vertex( 2 )->point();
    viewGouraudTriangle( a, b, c, 
			 (double) avt.value( fh->vertex( 0 ) ),
			 (double) avt.value( fh->vertex( 1 ) ),
			 (double) avt.value( fh->vertex( 2 ) ) );
  }

  template <typename FFunctor, typename FxFunctor, typename FyFunctor>
  void viewFunctorTriangle( const FFunctor & f, 
			    const FxFunctor & fx, const FyFunctor & fy, 
			    Point2 a, Point2 b, Point2 c )
  {
    double val_a = f( a.x(), a.y() );
    double val_b = f( b.x(), b.y() );
    double val_c = f( c.x(), c.y() );
    Vector2 gradf_a( fx( a.x(), a.y() ), fy( a.x(), a.y() ) );
    Vector2 gradf_b( fx( b.x(), b.y() ), fy( b.x(), b.y() ) );
    Vector2 gradf_c( fx( c.x(), c.y() ), fy( c.x(), c.y() ) );
    Vector2 ab = b - a;
    Vector2 bc = c - b;
    Vector2 ca = a - c;
    Point2 mid_ab = a + (ab/2.0);
    Point2 mid_bc = b + (bc/2.0);
    Point2 mid_ca = c + (ca/2.0);
    ab = ab / sqrt( ab.squared_length() );
    bc = bc / sqrt( bc.squared_length() );
    ca = ca / sqrt( ca.squared_length() );
    double val_ab = ( gradf_a * ab - gradf_b * ab + 4.0 * ( val_a + val_b ) ) / 8.0;
    double val_bc = ( gradf_b * bc - gradf_c * bc + 4.0 * ( val_b + val_c ) ) / 8.0;
    double val_ca = ( gradf_c * ca - gradf_a * ca + 4.0 * ( val_c + val_a ) ) / 8.0;
    viewGouraudTriangle( a, mid_ab, mid_ca, val_a, val_ab, val_ca );
    viewGouraudTriangle( b, mid_bc, mid_ab, val_b, val_bc, val_ab );
    viewGouraudTriangle( c, mid_ca, mid_bc, val_c, val_ca, val_bc );
    viewGouraudTriangle( mid_ab, mid_bc, mid_ca, val_ab, val_bc, val_ca );
  }

  void viewAVTFlatTriangle( AVT & avt, FaceHandle fh )
  {
    Point2 a = fh->vertex( 0 )->point();
    Point2 b = fh->vertex( 1 )->point();
    Point2 c = fh->vertex( 2 )->point();
    double val = (double) avt.value( fh );
    viewFlatTriangle( a, b, c, val );
  }

  void viewAVT( AVT & avt, Mode mode )
  {
    cairo_set_operator( _cr,  CAIRO_OPERATOR_ADD );
    _redf   = ( mode == Red )   || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
    _greenf = ( mode == Green ) || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
    _bluef  = ( mode == Blue )  || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
    cairo_set_line_width( _cr, 0.0 ); 
    cairo_set_line_cap( _cr, CAIRO_LINE_CAP_BUTT );
    cairo_set_line_join( _cr, CAIRO_LINE_JOIN_BEVEL );
    for ( FiniteFacesIterator it = avt.T().finite_faces_begin(), itend = avt.T().finite_faces_end();
	  it != itend; ++it )
      {
	if ( _gouraud ) viewAVTGouraudTriangle( avt, it );
	else            viewAVTFlatTriangle( avt, it );
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

  template <typename ImageDerivatives>
  void viewAVTWithDerivatives( AVT & avt, const ImageDerivatives & dImg, Mode mode )
  {
    cairo_set_operator( _cr,  CAIRO_OPERATOR_ADD );
    _redf   = ( mode == Red )   || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
    _greenf = ( mode == Green ) || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
    _bluef  = ( mode == Blue )  || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
    cairo_set_line_width( _cr, 0.0 ); 
    cairo_set_line_cap( _cr, CAIRO_LINE_CAP_BUTT );
    cairo_set_line_join( _cr, CAIRO_LINE_JOIN_BEVEL );
    typename ImageDerivatives::FFunctor   f( dImg );
    typename ImageDerivatives::FxFunctor fx( dImg );
    typename ImageDerivatives::FyFunctor fy( dImg );
    for ( FiniteFacesIterator it = avt.T().finite_faces_begin(), itend = avt.T().finite_faces_end();
	  it != itend; ++it )
      {
	Point2 a = it->vertex( 0 )->point();
	Point2 b = it->vertex( 1 )->point();
	Point2 c = it->vertex( 2 )->point();
	viewFunctorTriangle( f, fx, fy, a, b, c );
      }
  }  

  inline int i( double x ) const
  {
    return (int)round( (x+0.5) * _xf ) - _x0;
  }

  inline int j( double y ) const
  {
    return _height - (int)(round( (y+0.5) * _yf ) - _y0) - 1;
  }

};


///////////////////////////////////////////////////////////////////////////////
double randomUniform()
{
  return (double) random() / (double) RAND_MAX;
}
///////////////////////////////////////////////////////////////////////////////
template <typename Int>
struct Reverse {
  typedef Reverse<Int> Self;
  Int _i;
  inline Reverse() {}
  inline Reverse( Int i ) : _i( i ) {}
  inline Reverse( const Self & other ) : _i( other._i ) {}
  inline Self & operator=( const Self & other ) 
  { _i = other._i; return *this; }
  inline operator Int() const { return _i; }
  inline bool operator==( const Self & other ) const
  { return _i == other._i; }
  inline bool operator<( const Self & other ) const
  { return _i > other._i; }
  inline bool operator<=( const Self & other ) const
  { return _i >= other._i; }
  inline bool operator>( const Self & other ) const
  { return _i < other._i; }
  inline bool operator>=( const Self & other ) const
  { return _i <= other._i; }
  inline bool operator!=( const Self & other ) const
  { return _i != other._i; }
};

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;
///////////////////////////////////////////////////////////////////////////////

template <typename AVT>
void computeFullRelativeHull( AVT & avt, int l, std::string txt )
{
  using namespace DGtal;
  txt = "Compute full relative hull -- " + txt;
  trace.beginBlock( txt.c_str() );
  unsigned int pass = 0;
  bool changes = true;
  do {
    trace.info() << "- Pass " << pass << std::endl;
    changes = avt.fullRelativeHull();
    ++pass;
    if ( pass >= l ) break;
  } while ( changes );
  trace.endBlock();
}

template <typename AVT>
void viewAffineValuedTriangulation( AVT & avt, double b, 
				    double x0, double y0,
				    double x1, double y1,
				    bool gouraud,
				    std::string fname )
{
  using namespace DGtal;
  trace.beginBlock("View affine valued triangulation");
  CairoViewerAVT< AVT, int > cviewer
    ( (int) round( x0 ), (int) round( y0 ), 
      (int) round( (x1+1 - x0) * b ), (int) round( (y1+1 - y0) * b ), 
      b, b, gouraud );
  cviewer.viewAVT( avt, CairoViewerAVT< AVT, int>::Gray );
  cviewer.save( fname.c_str() );
  trace.endBlock();
}

template <typename AVT, typename ImageDerivatives>
void viewAVTWithDerivatives( AVT & avt, const ImageDerivatives & dImg,
			     double b, 
			     double x0, double y0,
			     double x1, double y1,
			     std::string fname )
{
  using namespace DGtal;
  trace.beginBlock("View affine valued triangulation with derivatives");
  CairoViewerAVT< AVT, int > cviewer
    ( (int) round( x0 ), (int) round( y0 ), 
      (int) round( (x1+1 - x0) * b ), (int) round( (y1+1 - y0) * b ), 
      b, b, true );
  cviewer.viewAVTWithDerivatives( avt, dImg, CairoViewerAVT< AVT, int>::Gray );
  cviewer.save( fname.c_str() );
  trace.endBlock();
}

template <typename AVT, typename ImageDerivatives>
void viewAVTAll( AVT & avt, const ImageDerivatives & dImg,
		 double b, 
		 double x0, double y0,
		 double x1, double y1,
		 std::string fname )
{
  viewAffineValuedTriangulation( avt, b, x0, y0, x1, y1, false, fname + ".png" );
  viewAffineValuedTriangulation( avt, b, x0, y0, x1, y1, true, fname + "-g.png" );
  viewAVTWithDerivatives( avt, dImg, b, x0, y0, x1, y1, fname + "-g2.png" );
}

template <typename AVT>
void viewAffineValuedTriangulationColor
( AVT & avt_red, AVT & avt_green, AVT & avt_blue,
  double b, 
  double x0, double y0,
  double x1, double y1,
  bool gouraud,
  std::string fname )
{
  using namespace DGtal;
  trace.beginBlock("View affine valued triangulation");
  CairoViewerAVT< AVT, int > cviewer
    ( (int) round( x0 ), (int) round( y0 ), 
      (int) round( (x1+1 - x0) * b ), (int) round( (y1+1 - y0) * b ), 
      b, b, gouraud );
  cviewer.viewAVT( avt_red, CairoViewerAVT< AVT, int>::Red );
  cviewer.viewAVT( avt_green, CairoViewerAVT< AVT, int>::Green );
  cviewer.viewAVT( avt_blue, CairoViewerAVT< AVT, int>::Blue );
  cviewer.save( fname.c_str() );
  trace.endBlock();
}

template <typename AVT, 
	  typename RedDerivatives, typename GreenDerivatives, typename BlueDerivatives>
void viewAVTColorWithDerivatives
( AVT & avt_red, AVT & avt_green, AVT & avt_blue,
  const RedDerivatives & dRed,
  const GreenDerivatives & dGreen,
  const BlueDerivatives & dBlue,
  double b, 
  double x0, double y0,
  double x1, double y1,
  std::string fname )
{
  using namespace DGtal;
  trace.beginBlock("View affine valued triangulation with derivatives");
  CairoViewerAVT< AVT, int > cviewer
    ( (int) round( x0 ), (int) round( y0 ), 
      (int) round( (x1+1 - x0) * b ), (int) round( (y1+1 - y0) * b ), 
      b, b, true );
  cviewer.viewAVTWithDerivatives( avt_red, dRed, CairoViewerAVT< AVT, int>::Red );
  cviewer.viewAVTWithDerivatives( avt_green, dGreen, CairoViewerAVT< AVT, int>::Green );
  cviewer.viewAVTWithDerivatives( avt_blue, dBlue, CairoViewerAVT< AVT, int>::Blue );
  cviewer.save( fname.c_str() );
  trace.endBlock();
}

template <typename AVT, 
	  typename RedDerivatives, typename GreenDerivatives, typename BlueDerivatives>
void viewAVTColorAll
( AVT & avt_red, AVT & avt_green, AVT & avt_blue,
  const RedDerivatives & dRed,
  const GreenDerivatives & dGreen,
  const BlueDerivatives & dBlue,
  double b, 
  double x0, double y0,
  double x1, double y1,
  std::string fname )
{
  viewAffineValuedTriangulationColor( avt_red, avt_green, avt_blue,
				      b, x0, y0, x1, y1, false, fname + ".png" );
  viewAffineValuedTriangulationColor( avt_red, avt_green, avt_blue,
				      b, x0, y0, x1, y1, true, fname + "-g.png" );
  viewAVTColorWithDerivatives( avt_red, avt_green, avt_blue, dRed, dGreen, dBlue,
			       b, x0, y0, x1, y1, fname + "-g2.png" );
}

template <typename OutImage, typename Derivatives>
void viewDerivatives( OutImage & dImage, const Derivatives & der, const std::string & s )
{
  typedef typename OutImage::Domain Domain;
  for ( typename Domain::ConstIterator it = dImage.domain().begin(), ite = dImage.domain().end();
        it != ite; ++it )
    {
      double val = der.dx( (double) (*it)[ 0 ] / 1.0, (double) (*it)[ 1 ] / 1.0 );
      val = std::min( 127.0, std::max( -128.0, val ) ) + 128.0;
      dImage.setValue( *it, (int) round(val) );
    }
  DGtal::PPMWriter<OutImage, GrayToRedGreen>::exportPPM( s + "-dx.ppm", dImage );
  for ( typename Domain::ConstIterator it = dImage.domain().begin(), ite = dImage.domain().end();
        it != ite; ++it )
    {
      double val = der.dy( (double) (*it)[ 0 ] / 1.0, (double) (*it)[ 1 ] / 1.0 );
      val = std::min( 127.0, std::max( -128.0, val ) ) + 128.0;
      dImage.setValue( *it, (int) round(val) );
    }
  DGtal::PPMWriter<OutImage, GrayToRedGreen>::exportPPM( s + "-dy.ppm", dImage );
}

template <typename OutImage, typename AVT>
void viewAVTGradient( OutImage & dImage, AVT & avt, const std::string & s )
{
  typedef typename OutImage::Domain Domain;
  for ( typename Domain::ConstIterator it = dImage.domain().begin(), ite = dImage.domain().end();
        it != ite; ++it )
    {
      double val = avt.gradientPixel( (double) (*it)[ 0 ], 
                                      (double) (*it)[ 1 ], 1.0 ).x();
      val = std::min( 127.0, std::max( -128.0, val ) ) + 128.0;
      dImage.setValue( *it, (int) round(val) );
    }
  DGtal::PPMWriter<OutImage, GrayToRedGreen>::exportPPM( s + "-adx.ppm", dImage );
  for ( typename Domain::ConstIterator it = dImage.domain().begin(), ite = dImage.domain().end();
        it != ite; ++it )
    {
      double val = avt.gradientPixel( (double) (*it)[ 0 ], 
                                      (double) (*it)[ 1 ], 1.0 ).y();
      val = std::min( 127.0, std::max( -128.0, val ) ) + 128.0;
      dImage.setValue( *it, (int) round(val) );
    }
  DGtal::PPMWriter<OutImage, GrayToRedGreen>::exportPPM( s + "-ady.ppm", dImage );
}

template <typename Value>
int affineValuedTriangulation( po::variables_map & vm )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel2;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel2> Triangulation2;
  typedef Triangulation2::Point                               Point2;
  typedef DGtal::Z2i::Space Space;
  typedef DGtal::Z2i::Domain Domain;
  typedef DGtal::Z2i::Point DPoint;
  typedef DGtal::Z2i::KSpace KSpace;
  typedef DGtal::Z2i::Integer Integer;
  typedef KSpace::SCell SCell;
  typedef DGtal::Z3i::Point                                   PointZ3;
  typedef AVT<Triangulation2, Kernel2, Value>                 AVTriangulation2;

  using namespace DGtal;

  double x0, y0, x1, y1;
  trace.beginBlock("Construction of the triangulation");
  AVTriangulation2 avt_red( -1 );
  AVTriangulation2 avt_green( -1 );
  AVTriangulation2 avt_blue( -1 );
  bool gouraud = vm.count( "gouraud" );
  double b = vm[ "bitmap" ].as<double>();
  double ratio = vm[ "random" ].as<double>();
  bool color = false;
  if ( ! vm.count( "image" ) ) return 1;

  typedef ImageSelector < Z2i::Domain, unsigned int>::Type Image;
  typedef Image::Domain Domain;
  ColorToRedFunctor c2r;
  ColorToGreenFunctor c2g;
  ColorToBlueFunctor c2b;
  GrayToGrayFunctor g2g;

  std::string imageFileName = vm[ "image" ].as<std::string>();
  Image image = GenericReader<Image>::import( imageFileName ); 
  std::string extension = imageFileName.substr(imageFileName.find_last_of(".") + 1);
  if ( extension == "ppm" ) color = true;
  x0 = 0.0; y0 = 0.0;
  x1 = (double) image.domain().upperBound()[ 0 ];
  y1 = (double) image.domain().upperBound()[ 1 ];
  if ( color )
    for ( Domain::ConstIterator it = image.domain().begin(), ite = image.domain().end();
          it != ite; ++it )
      {
        Point2 pt2( (*it)[ 0 ], (*it)[ 1 ] );
        unsigned int val = image( *it );
        if ( randomUniform() <= ratio )
          {
            avt_red.add  ( pt2, c2r( val ) );
            avt_green.add( pt2, c2g( val ) );
            avt_blue.add ( pt2, c2b( val ) );
          }
      }
  else
    for ( Domain::ConstIterator it = image.domain().begin(), ite = image.domain().end();
          it != ite; ++it )
      {
        Point2 pt2( (*it)[ 0 ], (*it)[ 1 ] );
        unsigned int val = image( *it );
        if ( randomUniform() <= ratio )
          avt_blue.add ( pt2, g2g( val ) );
      }
  trace.endBlock();

  trace.beginBlock("Computes derivatives.");
  typedef Derivatives<Image,ColorToRedFunctor>   RedDerivatives;
  typedef Derivatives<Image,ColorToGreenFunctor> GreenDerivatives;
  typedef Derivatives<Image,ColorToBlueFunctor>  BlueDerivatives;
  typedef Derivatives<Image,GrayToGrayFunctor>   GrayDerivatives;
  RedDerivatives redDer( image, c2r );
  GreenDerivatives greenDer( image, c2g );
  BlueDerivatives blueDer( image, c2b );
  GrayDerivatives grayDer( image, g2g );
  DPoint size( image.domain().upperBound() - image.domain().lowerBound() );
  // size *= 2;
  typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> GrayImage2D;
  GrayImage2D dImage( Domain( DPoint( 0, 0 ), size ) );
  
  if ( color ) {
    viewDerivatives( dImage, redDer, "image-red" );
    viewDerivatives( dImage, greenDer, "image-green" );
    viewDerivatives( dImage, blueDer, "image-blue" );
  } else {
    viewDerivatives( dImage, grayDer, "image-gray" );
  }
  trace.endBlock();

  if ( vm.count( "saddle" ) )
    {
      trace.beginBlock("Process saddle points");
      if ( color ) {
        avt_red.processSaddlePoints();
        avt_green.processSaddlePoints();
      }
      avt_blue.processSaddlePoints();
      trace.endBlock();
    }

  if ( color ) viewAVTColorAll( avt_red, avt_green, avt_blue, 
				redDer, greenDer, blueDer,
				b, x0, y0, x1, y1, "avt-before" );
  else         viewAVTAll( avt_blue, grayDer, 
			   b, x0, y0, x1, y1, "avt-before" );
  if ( color )
    {
      computeFullRelativeHull( avt_red, vm[ "limit" ].as<int>(), "RED" );
      computeFullRelativeHull( avt_green, vm[ "limit" ].as<int>(), "GREEN" );
    }
  computeFullRelativeHull( avt_blue, vm[ "limit" ].as<int>(), "BLUE" );

  if ( color ) viewAVTColorAll( avt_red, avt_green, avt_blue, 
				redDer, greenDer, blueDer,
				b, x0, y0, x1, y1, "avt-after" );
  else         viewAVTAll( avt_blue, grayDer,
			   b, x0, y0, x1, y1, "avt-after" );
  if ( color ) {
    viewAVTGradient( dImage, avt_red, "avt-after-red" );
    viewAVTGradient( dImage, avt_green, "avt-after-green" );
    viewAVTGradient( dImage, avt_blue, "avt-after-blue" );
  } else {
    viewAVTGradient( dImage, avt_blue, "avt-after-gray" );
  }
  
  if ( vm.count( "compress" ) )
    {
      double ratio = vm[ "compress" ].as<double>();
      int sub = 2;
      if ( color ) {
        trace.beginBlock("Compressing triangulation -- RED");
	avt_red.compress( ratio, gouraud, sub );
	trace.endBlock();
	trace.beginBlock("Compressing triangulation -- GREEN");
	avt_green.compress( ratio, gouraud, sub );
	trace.endBlock();
      }
      trace.beginBlock("Compressing triangulation -- BLUE");
      avt_blue.compress( ratio, gouraud, sub );
      trace.endBlock();
      if ( color ) viewAVTColorAll( avt_red, avt_green, avt_blue, 
				    redDer, greenDer, blueDer,
				    b, x0, y0, x1, y1, "avt-compressed" );
      else         viewAVTAll( avt_blue, grayDer,
			       b, x0, y0, x1, y1, "avt-compressed" );
      if ( color )
	{
	  computeFullRelativeHull( avt_red, vm[ "limit" ].as<int>(), "RED" );
	  computeFullRelativeHull( avt_green, vm[ "limit" ].as<int>(), "GREEN" );
	}
      computeFullRelativeHull( avt_blue, vm[ "limit" ].as<int>(), "BLUE" );
      if ( color ) viewAVTColorAll( avt_red, avt_green, avt_blue, 
				    redDer, greenDer, blueDer,
				    b, x0, y0, x1, y1, "avt-compressed-rhull" );
      else         viewAVTAll( avt_blue, grayDer,
			       b, x0, y0, x1, y1, "avt-compressed-rhull" );

    }

  return 0;
}

int main( int argc, char** argv )
{

  using namespace DGtal;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("image,i", po::value<std::string>(), "Specifies the input shape as a 2D image PPM filename.")
    ("random,r", po::value<double>()->default_value(1.0), "Keep only a proportion [arg] (between 0 and 1) of the input data point.")  
    ("limit,L", po::value<int>()->default_value(1000), "Gives the maximum number of passes (default is 10000).")  
    ("reverse,R", "Reverses the topology of the image (black has more priority than white).") 
    ("bitmap,b", po::value<double>()->default_value( 2.0 ), "Rasterization magnification factor [arg] for PNG export." )
    ("gouraud,g", "Displays faces with Gouraud-like shading.")  
    ("2ndorder,2", "Displays faces with 2nd order Gouraud-like shading.")  
    ("saddle,s", "Process saddle points.")  
    ("compress,z", po::value<double>(), "Compress image with compress ratio [arg], e.g. arg=0.9 means keeping only 10% of vertices.")  
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
  if( ! parseOK || vm.count("help") || argc <= 1 )
    {
      trace.info()<< "Generate a 2D triangulation from an arbitrary set of points. The triangulation provides a kind of piecewise linear approximation of any arbitrary color image." <<std::endl << "Basic usage: " << std::endl
		  << "\t2d-triangulation-color [options] -i <image.ppm> -b 4"<<std::endl
		  << general_opt << "\n";
      return 0;
    }

  int result = vm.count( "reverse" ) 
    ? affineValuedTriangulation< Reverse<int> >( vm )
    : affineValuedTriangulation<int>( vm );

  return result;
}
