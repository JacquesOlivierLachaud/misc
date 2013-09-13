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
#include <DGtal/geometry/curves/ArithmeticalDSS.h>
#include <DGtal/geometry/curves/SaturatedSegmentation.h>
#include "DGtal/io/readers/PointListReader.h"

#include <CGAL/Delaunay_triangulation_2.h>
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

};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/**
   This class is intended for visualizing Affine Valued triangulation.
*/
template <typename TAVT, typename TValue>
class ViewerAVT
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
  typedef DGtal::Board2D                       Board;
  typedef DGtal::Color                         Color;
  typedef DGtal::Z2i::Point                    PointZ2;

private:
  Board & _board;
  AVT & _avt;
  Value _min, _max;
  bool _gouraud;
  std::vector<Color> _label_colors;
  std::vector<Color> _other_colors;

public:
  /**
     Constructor. Requires a board \a board for display and an affine
     valued triangulation \a avt.
  */
  ViewerAVT( DGtal::Alias<Board> board, 
             DGtal::Alias<AVT> avt,
             Value min, Value max, 
	     bool gouraud )
    : _board( board ), _avt( avt ),
      _min( min ), _max( max ), _gouraud( gouraud )
  {
    DGtal::GrayscaleColorMap<double> grayShade( (double) min, (double) max );
    for ( Value i = min; i <= max; i++ )
      _label_colors.push_back( grayShade( (double) i ) );
    _other_colors.push_back( Color::Blue );
  }

  inline Color color( Value val ) const
  {
    if ( val == (Value) _avt.invalid() )
      return _other_colors[ 0 ];
    else
      return _label_colors[ std::min( _max, std::max( _min, val ) ) - _min ];
  }

  /**
     View vertices of the triangulation. 
  */
  void viewVertices()
  {
    PointZ2 dummy;
    std::string specificStyle =  dummy.className() + "/Grid";
    _board << DGtal::SetMode( dummy.className(), "Grid" );
    for ( FiniteVerticesIterator it = _avt.T().finite_vertices_begin(), itend = _avt.T().finite_vertices_end();
          it != itend; ++it )
      {
        VertexHandle vh = it;
        Value l1 = (Value) _avt.value( vh );
        Color c = color( l1 );
        PointZ2 a = toDGtal( it->point() );
        _board << DGtal::CustomStyle( specificStyle, new DGtal::CustomColors( c, c ) );
        _board << a;
      }
  }

  /**
     View edges of the triangulation. 
  */
  void viewEdges()
  {
    for ( FiniteEdgesIterator it = _avt.T().finite_edges_begin(), itend = _avt.T().finite_edges_end();
          it != itend; ++it )
      {
        Edge e = *it;
        Value l1 = (Value) _avt.value( e );
        Color col = color( l1 );
        _board.setPenColor( col );
        _board.setFillColor( col );
        _board.setLineWidth( 2.0f );
        PointZ2 a = toDGtal( _avt.TH.source( e )->point() );
        PointZ2 b = toDGtal( _avt.TH.target( e )->point() );
        _board.drawLine( a[ 0 ], a[ 1 ], b[ 0 ], b[ 1 ] );
      }
  }

  /**
     View edges of the triangulation. 
  */
  void viewTriangulation( Color col )
  {
    _board.setPenColor( col );
    _board.setFillColor( col );
    _board.setLineWidth( 1.0f );
    for ( FiniteEdgesIterator it = _avt.T().finite_edges_begin(), itend = _avt.T().finite_edges_end();
          it != itend; ++it )
      {
        Edge e = *it;
        Value l1 = (Value) _avt.value( e );
        PointZ2 a = toDGtal( _avt.TH.source( e )->point() );
        PointZ2 b = toDGtal( _avt.TH.target( e )->point() );
        _board.drawLine( a[ 0 ], a[ 1 ], b[ 0 ], b[ 1 ] );
      }
  }

  /**
     View faces of the triangulation. 
  */
  void viewFaces()
  {
    for ( FiniteFacesIterator it = _avt.T().finite_faces_begin(), itend = _avt.T().finite_faces_end();
          it != itend; ++it )
      {
	_board.setLineWidth( 0.0f );
        FaceHandle fh = it;
	if ( _gouraud )
	  {
	    Value lf = (Value) _avt.value( fh );
	    Value l0 = (Value) _avt.value( fh->vertex( 0 ) );
	    Value l1 = (Value) _avt.value( fh->vertex( 1 ) );
	    Value l2 = (Value) _avt.value( fh->vertex( 2 ) );
	    PointZ2 a = toDGtal( fh->vertex( 0 )->point() );
	    PointZ2 b = toDGtal( fh->vertex( 1 )->point() );
	    PointZ2 c = toDGtal( fh->vertex( 2 )->point() );
	    _board.setPenColor( color( lf ) );
	    _board.setFillColor( color( lf ) );
	    _board.drawTriangle( a[ 0 ], a[ 1 ], b[ 0 ], b[ 1 ], c[ 0 ], c[ 1 ] );
	    _board.setPenColor( DGtal::Color( 0, 0, 0, 255 ) );
	    _board.fillGouraudTriangle( a[ 0 ], a[ 1 ], lf <= l0 ? color( l0 ) : color( lf ),
	     				b[ 0 ], b[ 1 ], lf <= l1 ? color( l1 ) : color( lf ),
	     				c[ 0 ], c[ 1 ], lf <= l2 ? color( l2 ) : color( lf ),
					2 );
	  }
	else
	  {
	    Value l1 = (Value) _avt.value( fh );
	    Color col = color( l1 );
	    _board.setPenColor( col );
	    _board.setFillColor( col );
	    PointZ2 a = toDGtal( fh->vertex( 0 )->point() );
	    PointZ2 b = toDGtal( fh->vertex( 1 )->point() );
	    PointZ2 c = toDGtal( fh->vertex( 2 )->point() );
	    _board.drawTriangle( a[ 0 ], a[ 1 ], b[ 0 ], b[ 1 ], c[ 0 ], c[ 1 ] );
	  }
      }
  }

  /**
     View edges of the triangulation. 
  */
  void viewContour( Value val )
  {
    PointZ2 dummy;
    std::string specificStyle =  dummy.className() + "/Grid";
    _board << DGtal::SetMode( dummy.className(), "Grid" );
    std::vector<Edge> in_edges, out_edges;
    for ( FiniteVerticesIterator it = _avt.T().finite_vertices_begin(), itend = _avt.T().finite_vertices_end();
          it != itend; ++it )
      {
        VertexHandle vh = it;
        Value l1 = (Value) _avt.value( vh );
	if ( _avt.localIsocontour( in_edges, out_edges, vh, val ) )
	  {
	    Color col = Color::Green;
	    PointZ2 b = toDGtal( it->point() );
	    _board << DGtal::CustomStyle( specificStyle, new DGtal::CustomColors( col, col ) );
	    _board << b;
	    _board.setPenColor( col );
	    _board.setFillColor( col );
	    _board.setLineWidth( 2.0f );
	    for ( unsigned int i = 0; i < in_edges.size(); ++i )
	      {
		PointZ2 a = toDGtal( _avt.TH.source( in_edges[ i ] )->point() );
		PointZ2 c = toDGtal( _avt.TH.target( out_edges[ i ] )->point() );
		_board.drawLine( a[ 0 ], a[ 1 ], b[ 0 ], b[ 1 ] );
		_board.drawLine( b[ 0 ], b[ 1 ], c[ 0 ], c[ 1 ] );
	      }
	  }
      }
  }

  /**
     View edges of the triangulation. 
  */
  void viewCorners( Value val )
  {
    PointZ2 dummy;
    std::string specificStyle =  dummy.className() + "/Grid";
    _board << DGtal::SetMode( dummy.className(), "Grid" );
    std::vector<Edge> in_edges, out_edges;
    for ( FiniteVerticesIterator it = _avt.T().finite_vertices_begin(), itend = _avt.T().finite_vertices_end();
          it != itend; ++it )
      {
        VertexHandle vh = it;
        Value l1 = (Value) _avt.value( vh );
	if ( _avt.localIsocontour( in_edges, out_edges, vh, val ) )
	  {
	    Color col = Color::Magenta;
	    for ( unsigned int i = 0; i < in_edges.size(); ++i )
	      {
		if ( ! _avt.isSmooth( in_edges[ i ], out_edges[ i ], 1.0 ) )
		  {
		    PointZ2 b = toDGtal( it->point() );
		    _board << DGtal::CustomStyle( specificStyle, new DGtal::CustomColors( col, col ) );
		    _board << b;
		    _board.setPenColor( col );
		    _board.setFillColor( col );
		    _board.setLineWidth( 2.0f );
		    PointZ2 a = toDGtal( _avt.TH.source( in_edges[ i ] )->point() );
		    PointZ2 c = toDGtal( _avt.TH.target( out_edges[ i ] )->point() );
		    _board.drawLine( (a[ 0 ] + b[ 0 ])/2.0, (a[ 1 ] + b[ 1 ])/2.0, b[ 0 ], b[ 1 ] );
		    _board.drawLine( b[ 0 ], b[ 1 ], (b[ 0 ] + c[ 0 ])/2.0, (b[ 1 ] + c[ 1 ])/2.0 );
		  }
	      }
	  }
      }
  }
  

};



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

  void viewAVT( AVT & avt, Mode mode )
  {
    cairo_operator_t op;
    switch ( mode ) {
    case Red: 
    case Green:
    case Blue: 
      op = CAIRO_OPERATOR_ADD; break;
    case Gray:
    default: op = CAIRO_OPERATOR_SOURCE; break;
    };
    op = CAIRO_OPERATOR_ADD;
    cairo_set_operator( _cr, op );
    double redf   = ( mode == Red )   || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
    double greenf = ( mode == Green ) || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
    double bluef  = ( mode == Blue )  || ( mode == Gray ) ? 1.0 / 255.0 : 0.0;
    cairo_set_line_width( _cr, 0.0 ); 
    cairo_set_line_cap( _cr, CAIRO_LINE_CAP_BUTT );
    cairo_set_line_join( _cr, CAIRO_LINE_JOIN_BEVEL );
    if ( ! _gouraud )
      for ( FiniteFacesIterator it = avt.T().finite_faces_begin(), itend = avt.T().finite_faces_end();
	    it != itend; ++it )
	{
	  FaceHandle fh = it;
	  double val = (double) avt.value( fh );
	  cairo_set_source_rgb( _cr, val * redf, val * greenf, val * bluef );
	  Point2 a = fh->vertex( 0 )->point();
	  Point2 b = fh->vertex( 1 )->point();
	  Point2 c = fh->vertex( 2 )->point();
	  cairo_move_to( _cr, i( a.x() ), j( a.y() ) );
	  cairo_line_to( _cr, i( b.x() ), j( b.y() ) );
	  cairo_line_to( _cr, i( c.x() ), j( c.y() ) );
	  cairo_close_path( _cr );
	  cairo_fill( _cr );
	  //cairo_stroke( _cr );
	}
    else
      for ( FiniteFacesIterator it = avt.T().finite_faces_begin(), itend = avt.T().finite_faces_end();
	    it != itend; ++it )
	{
	  FaceHandle fh = it;
	  Point2 a = fh->vertex( 0 )->point();
	  Point2 b = fh->vertex( 1 )->point();
	  Point2 c = fh->vertex( 2 )->point();
	  cairo_pattern_t * pattern = cairo_pattern_create_mesh();
	  /* Add a Gouraud-shaded triangle */
	  cairo_mesh_pattern_begin_patch (pattern);
	  cairo_mesh_pattern_move_to (pattern, i( a.x() ), j( a.y() ) );
	  cairo_mesh_pattern_line_to (pattern, i( b.x() ), j( b.y() ) );
	  cairo_mesh_pattern_line_to (pattern, i( c.x() ), j( c.y() ) );
	  for ( int idx = 0; idx < 3; ++idx )
	    {
	      double val = (double) avt.value( fh->vertex( idx ) );
	      cairo_mesh_pattern_set_corner_color_rgb (pattern, idx, val * redf, val * greenf, val * bluef );
	    }
	  cairo_mesh_pattern_end_patch (pattern);
	  cairo_set_source( _cr, pattern );
	  cairo_move_to( _cr, i( a.x() ), j( a.y() ) );
	  cairo_line_to( _cr, i( b.x() ), j( b.y() ) );
	  cairo_line_to( _cr, i( c.x() ), j( c.y() ) );
	  cairo_close_path( _cr );
	  cairo_fill( _cr );
	  cairo_pattern_destroy( pattern );
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

template <typename Value>
int affineValuedTriangulation( po::variables_map & vm )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel2;
  typedef CGAL::Delaunay_triangulation_2<Kernel2>             Triangulation2;
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
  double ratio = vm[ "random" ].as<double>();
  if ( vm.count( "image" ) )
    {
      typedef ImageSelector < Z2i::Domain, unsigned int>::Type Image;
      typedef Image::Domain Domain;
      std::string imageFileName = vm[ "image" ].as<std::string>();
      Image image = GenericReader<Image>::import( imageFileName ); 
      x0 = 0.0; y0 = 0.0;
      x1 = (double) image.domain().upperBound()[ 0 ];
      y1 = (double) image.domain().upperBound()[ 1 ];
      for ( Domain::ConstIterator it = image.domain().begin(), ite = image.domain().end();
            it != ite; ++it )
        {
          Point2 pt2( (*it)[ 0 ], (*it)[ 1 ] );
          unsigned int val = image( *it );
          if ( randomUniform() <= ratio )
	    {
	      avt_red.add  ( pt2, (val >> 16) & 0xff );
	      avt_green.add( pt2, (val >> 8) & 0xff );
	      avt_blue.add ( pt2, val & 0xff );
	    }
        }
    }
  trace.endBlock();
  
  bool gouraud = vm.count( "gouraud" );

  {
    trace.beginBlock("View affine valued triangulation");
    double b = vm[ "bitmap" ].as<double>();
    CairoViewerAVT< AVTriangulation2, int > cviewer
      ( (int) round( x0 ), (int) round( y0 ), 
        (int) round( (x1+1 - x0) * b ), (int) round( (y1+1 - y0) * b ), 
        b, b, gouraud );
    cviewer.viewAVT( avt_red, CairoViewerAVT< AVTriangulation2, int>::Red );
    cviewer.viewAVT( avt_green, CairoViewerAVT< AVTriangulation2, int>::Green );
    cviewer.viewAVT( avt_blue, CairoViewerAVT< AVTriangulation2, int>::Blue );
    cviewer.save( "avt-before.png" );
    trace.endBlock();
  }

  
  {
    trace.beginBlock("Compute full relative hull -- RED");
    unsigned int pass = 0;
    bool changes = true;
    do {
      trace.info() << "- Pass " << pass << std::endl;
      changes = avt_red.fullRelativeHull();
      ++pass;
      if ( pass >= vm[ "limit" ].as<int>() ) break;
    } while ( changes );
    trace.endBlock();
  }
  {
    trace.beginBlock("Compute full relative hull -- GREEN");
    unsigned int pass = 0;
    bool changes = true;
    do {
      trace.info() << "- Pass " << pass << std::endl;
      changes = avt_green.fullRelativeHull();
      ++pass;
      if ( pass >= vm[ "limit" ].as<int>() ) break;
    } while ( changes );
    trace.endBlock();
  }
  {
    trace.beginBlock("Compute full relative hull -- BLUE");
    unsigned int pass = 0;
    bool changes = true;
    do {
      trace.info() << "- Pass " << pass << std::endl;
      changes = avt_blue.fullRelativeHull();
      ++pass;
      if ( pass >= vm[ "limit" ].as<int>() ) break;
    } while ( changes );
    trace.endBlock();
  }

  {
    trace.beginBlock("View affine valued triangulation");
    double b = vm[ "bitmap" ].as<double>();
    CairoViewerAVT< AVTriangulation2, int > cviewer
      ( (int) round( x0 ), (int) round( y0 ), 
        (int) round( (x1+1 - x0) * b ), (int) round( (y1+1 - y0) * b ), 
        b, b, gouraud );
    cviewer.viewAVT( avt_red, CairoViewerAVT< AVTriangulation2, int>::Red );
    cviewer.viewAVT( avt_green, CairoViewerAVT< AVTriangulation2, int>::Green );
    cviewer.viewAVT( avt_blue, CairoViewerAVT< AVTriangulation2, int>::Blue );
    cviewer.save( "avt-after.png" );
    trace.endBlock();
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
