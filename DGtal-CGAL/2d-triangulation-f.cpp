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
#include "DGtal/geometry/curves/FreemanChain.h"

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

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
  typedef typename std::map<Edge,Value>                     Edge2ValueMap;
  typedef typename std::map<FaceHandle,Value>               FaceHandle2ValueMap;
  typedef DGtal::UmbrellaPart2D<Triangulation2,Kernel2>     Strip;
  typedef typename std::set<Edge>                           EdgeSet;

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
    ASSERT( std::min( val_v0, val_v1 ) <= val_e );
    ASSERT( std::max( val_v0, val_v1 ) >= val_e );
    if ( val_v1 < val_v0 )
      {
	Strip strip( T() );
        NoGreaterThanValuePredicate predNoGreaterThanV1( *this, val_e );
        strip.init( predNoGreaterThanV1, predNoGreaterThanV1, e );
        if ( strip.isConcave() ) 
          inQueue.insert( std::make_pair( v0, v1 ) );
      }
    else if ( val_v0 < val_v1 )
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
        NoGreaterThanValuePredicate predNoGreaterThanV1( *this, val_V1 );
        if ( ! predNoGreaterThanV1( e ) ) 
	  {
	    inQueue.erase( inQueue.begin() );
	    continue; // Edge (v0,v1) has changed value.
	  }
	strip.init( predNoGreaterThanV1, predNoGreaterThanV1, e );
	ASSERT( strip.isValid() );
        if ( strip.isConcave() ) 
          {
            // We may already update the value of border edges.
            setValue( strip.e0(), std::max( value( strip.e0() ), val_V1+1 ) );
            setValue( strip.en(), std::max( value( strip.en() ), val_V1+1 ) );
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
                DGtal::trace.info() << "  - Flipping " << TH.source( fedge )->point() 
                                    << " -> " <<  TH.target( fedge )->point() << std::endl;
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
		    inQueue.erase( inQueue.begin() );
		  }
		changes = true; ++nb_flipped;
	      }
	    else // none are flippable. This is a concave/flat piece.
	      { // All faces and edges are set to value V1+1
                setValue( strip.f( 0 ), val_V1+1 );
                for ( unsigned int i = 1; i < strip.size() - 1; ++i )
                  {
                    setValue( strip.e( i ), val_V1+1 );
                    setValue( strip.f( i ), val_V1+1 );
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
    typename Edge2ValueMap::const_iterator it = _eFct.find( e );
    if ( it == _eFct.end() ) 
      {
        Value v = std::min( value( TH.source( e ) ), 
                            value( TH.target( e ) ) );
        _eFct[ e ] = v;
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
    _eFct[ e ] = val;
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
    _eFct.erase( e );
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
template <typename TAVT>
class ViewerAVT
{
public:
  typedef TAVT                                 AVT;
  typedef typename AVT::Triangulation2         Triangulation2;
  typedef typename AVT::Value                  Value;
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
  std::vector<Color> _label_colors;
  std::vector<Color> _other_colors;

public:
  /**
     Constructor. Requires a board \a board for display and an affine
     valued triangulation \a avt.
  */
  ViewerAVT( DGtal::Alias<Board> board, 
             DGtal::Alias<AVT> avt,
             Value min, Value max )
    : _board( board ), _avt( avt ),
      _min( min ), _max( max )
  {
    DGtal::GrayscaleColorMap<double> grayShade( (double) min, (double) max );
    for ( Value i = min; i <= max; i++ )
      _label_colors.push_back( grayShade( (double) i ) );
    _other_colors.push_back( Color::Blue );
  }

  inline Color color( Value val ) const
  {
    if ( val == _avt.invalid() )
      return _other_colors[ 0 ];
    else
      return _label_colors[ std::min( _max, std::max( _min, val ) ) ];
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
        Value l1 = _avt.value( vh );
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
        Value l1 = _avt.value( e );
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
        Value l1 = _avt.value( e );
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
        FaceHandle fh = it;
        Value l1 = _avt.value( fh );
        Color col = color( l1 );
        _board.setPenColor( col );
        _board.setFillColor( col );
        _board.setLineWidth( 0.0f );
        PointZ2 a = toDGtal( fh->vertex( 0 )->point() );
        PointZ2 b = toDGtal( fh->vertex( 1 )->point() );
        PointZ2 c = toDGtal( fh->vertex( 2 )->point() );
        _board.drawTriangle( a[ 0 ], a[ 1 ], b[ 0 ], b[ 1 ], c[ 0 ], c[ 1 ] );
      }
  }

};

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
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
  typedef AVT<Triangulation2, Kernel2, int>                   AVTriangulation2;

  using namespace DGtal;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("point-fct-list,l", po::value<std::string>(), "Specifies the input shape as a list of 2d integer points + value, 'x y f(x,y)' per line.");

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
      trace.info()<< "Generate a 2D triangulation from an arbitrary set of points. The triangulation provides a kind of piecewise linear approximation of the digitized set. The digitized set has label 0, the exterior points have label 1." <<std::endl << "Basic usage: " << std::endl
		  << "\t2d-triangulation-f [options] -l <point-list> -v -1"<<std::endl
		  << general_opt << "\n";
      return 0;
    }

  trace.beginBlock("Construction of the triangulation");
  AVTriangulation2 avt( -1 );
  int min_value = -1;
  int max_value = -1;
  if ( vm.count( "point-fct-list" ) )
    {
      PointZ3 lo, hi;
      std::vector<PointZ3> pts;
      if ( readPointList<PointZ3>( pts, lo, hi, vm["point-fct-list"].as<std::string>() ) != 0 )
        {
          trace.error() << "Error reading file <" << vm["point-fct-list"].as<std::string>() << ">." 
                        << std::endl;
          return 1;
        }
      // ks.init( lo, hi, true );
      for ( std::vector<PointZ3>::const_iterator it = pts.begin(), ite = pts.end();
            it != ite; ++it )
        {
          Point2 pt2( (*it)[ 0 ], (*it)[ 1 ] );
          int val = (*it)[ 2 ];
          if ( ( min_value == -1 ) || ( min_value > val ) ) min_value = val;
          if ( ( max_value == -1 ) || ( max_value < val ) ) max_value = val;
          avt.add( pt2, val );
        }
    }
  trace.endBlock();

  trace.beginBlock("View initial triangulation");
  Board2D board;
  ViewerAVT< AVTriangulation2 > viewer( board, avt, min_value, max_value );
  viewer.viewFaces();
  viewer.viewEdges();
  viewer.viewVertices();
  board.saveSVG("avt-before.svg");
  board.saveEPS("avt-before.eps");
  board.clear();
  trace.endBlock();
  
  trace.beginBlock("Compute full relative hull");
  unsigned int pass = 0;
  bool changes = true;
  do {
    trace.info() << "- Pass " << pass << std::endl;
    changes = avt.fullRelativeHull();
    viewer.viewFaces();
    viewer.viewEdges();
    viewer.viewVertices();
    std::ostringstream ostr;
    ostr << "tmp/avt-" << pass << ".eps";
    board.saveEPS( ostr.str().c_str() );
    ++pass;
  } while ( changes );
  trace.endBlock();

  trace.beginBlock("View affine valued triangulation");
  viewer.viewFaces();
  viewer.viewEdges();
  viewer.viewVertices();
  viewer.viewTriangulation( DGtal::Color::Red );
  board.saveSVG("avt-after.svg");
  board.saveEPS("avt-after.eps");
  board.clear();
  trace.endBlock();

  return 0;
}
