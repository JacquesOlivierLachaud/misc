#include <iostream>
#include <vector>
#include <string>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/Shapes.h>
#include <DGtal/shapes/ShapeFactory.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/geometry/curves/GridCurve.h>
#include <DGtal/io/boards/Board2D.h>
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include <DGtal/geometry/curves/ArithmeticalDSS.h>
#include <DGtal/geometry/curves/SaturatedSegmentation.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

#include "Auxiliary.h"

static const double EPSILON = 0.00001;
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
   A structure that gives a few methods for moving in a Triangulation.
*/
template <typename Triangulation, typename Kernel>
struct TriangulationHelper {
  typedef typename Triangulation::Vertex_circulator       VertexCirculator;
  typedef typename Triangulation::Vertex_handle           VertexHandle;
  typedef typename Triangulation::Edge_iterator           EdgeIterator;
  typedef typename Triangulation::Edge_circulator         EdgeCirculator;
  typedef typename Triangulation::Edge                    Edge;
  typedef typename Triangulation::Finite_faces_iterator   FacesIterator;
  typedef typename Triangulation::Face_handle             FaceHandle;
  typedef typename Triangulation::Point                   Point;
  typedef typename Kernel::Vector_2                       Vector;
  typedef typename Kernel::FT                             Coordinate;
  typedef typename Kernel::FT                             Component;

  const Triangulation & _T;

public:
  TriangulationHelper( const Triangulation & T )
    : _T( T ) {}

  inline const Triangulation & T() const
  { return _T; }

  inline
  Component det( const Vector & u, const Vector & v ) const
  {
    return u.x() * v.y() - u.y() * v.x();
  }

  /**
     @param e an oriented edge [st] within a face [stu].
     @return the (positive) angle of [st] and [su].
  */
  inline 
  Component innerAngle( const Edge & e ) const
  {
    Vector v1 = target( e )->point() - source( e )->point();
    Vector v2 = ( e.first->vertex( e.second ) )->point() - source( e )->point();
    return acos( v1 * v2 / sqrt( v1.squared_length() ) / sqrt( v2.squared_length() ) );
  }

  /// The incident face is considered ccw.
  VertexHandle source( const Edge & e ) const
  {
    return e.first->vertex( T().ccw( e.second ) );
  }
  
  /// The incident face is considered ccw.
  VertexHandle target( const Edge & e ) const
  {
    return e.first->vertex( T().cw( e.second ) );
  }

  Edge nextCCWAroundFace( const Edge & e ) const
  {
    return Edge( e.first, T().ccw( e.second ) );
  }
  
  
  Edge nextCWAroundFace( const Edge & e ) const
  {
    return Edge( e.first, T().cw( e.second ) );
  }
  
  /// @return the next edge around the source vertex, such that this
  /// edge has the same source and is CCW around it.
  Edge nextCCWAroundSourceVertex( const Edge & e ) const
  {
    Edge e1 = nextCWAroundFace( e );
    Edge n = T().mirror_edge( e1 );
    if ( source( n ) != source( e ) )
      DGtal::trace.error() << "[DigitalCore::nextCCWAroundSourceVertex] sources are not consistent." << std::endl;
    if ( target( n ) == target( e ) )
      DGtal::trace.error() << "[DigitalCore::nextCCWAroundSourceVertex] targets are equal." << std::endl;
    return n;
  }
  
  /// @return the next edge around the source vertex, such that this
  /// edge has the same source and is CW around it.
  Edge nextCWAroundSourceVertex( const Edge & e ) const
  {
    return nextCCWAroundFace( T().mirror_edge( e ) );
  }

  /// The order is ccw around some fictive point from v1 to v4.
  /// @return 'true' if the convex hull of [v1,v2,v3,v4] is a convex quadrilateral.
  bool isConvexQuadrilateral( const VertexHandle & v1, 
                              const VertexHandle & v2,
                              const VertexHandle & v3,
                              const VertexHandle & v4 ) const
  {
    Point a( v1->point() );
    Point b( v2->point() );
    Point c( v3->point() );
    Point d( v4->point() );
    return ( det( ( b-a ), ( c-b ) ) > 0 )
      && ( det( ( c-b ), ( d-c ) ) > 0 )
      && ( det( ( d-c ), ( a-d ) ) > 0 )
      && ( det( ( a-d ), ( b-a ) ) > 0 );
  }

  /**
     @return 'true' iff the edge e is stabbed by the segment [v1v2],
     assuming the following point locations.
     .....x......
     .....|e.....
     .v1..|......
     .x...|...v2.
     .....|...x..
     .....|......
     .....v......
  */
  bool isStabbedCCW( const Edge & e, VertexHandle v1, VertexHandle v2 ) const
  {
    return isConvexQuadrilateral( source( e ), v1, target( e ), v2 );
  }

  /**
     @return 'true' iff the edge e is stabbed by the segment [v1v2],
     without assumption on point locations.
     .....x......
     .....|e.....
     .v1..|......
     .x...|...v2.
     .....|...x..
     .....|......
     .....x......
  */
  bool isStabbed( const Edge & e, VertexHandle v1, VertexHandle v2 ) const
  {
    return isConvexQuadrilateral( source( e ), v1, target( e ), v2 )
      || isConvexQuadrilateral( source( e ), v2, target( e ), v1 );
  }

  /**
     @return 'true' iff [v1,v2] is an edge of the triangulation. In
     this case, \a e is this edge as an out parameter.
     @pre v1 and v2 are vertex of this triangulation.
  */
  bool findEdge( Edge & e, VertexHandle v1, VertexHandle v2 ) const
  {
    EdgeCirculator ci = T().incident_edges( v1 );
    if ( ci != 0 ) {
      e = *ci;
      if ( source( e ) != v1 ) e = T().mirror_edge( e );
      if ( source( e ) != v1 )
        DGtal::trace.error() << "[TriangulationHelper::isEdge] Invalid incident edge." << std::endl;
      Edge start = e;
      bool found = false;
      do {
        if ( target( e ) == v2 ) return true;
        e = nextCCWAroundSourceVertex( e );
      } while ( e != start );
      return false;
    }
    else return false;
  }
};

template <typename Triangulation, typename Kernel>
struct SimplicialStrip
{
public:
  typedef SimplicialStrip<Triangulation, Kernel>          Self;
  typedef typename Triangulation::Vertex_circulator       VertexCirculator;
  typedef typename Triangulation::Vertex_handle           VertexHandle;
  typedef typename Triangulation::Edge_iterator           EdgeIterator;
  typedef typename Triangulation::Edge                    Edge;
  typedef typename Triangulation::Finite_faces_iterator   FacesIterator;
  typedef typename Triangulation::Face_handle             FaceHandle;
  typedef typename Triangulation::Point                   Point;
  typedef typename Kernel::Vector_2                       Vector;
  typedef typename Kernel::FT                             Coordinate;
  typedef typename Kernel::FT                             Component;

public:
  /// The pivot v
  VertexHandle pivot; 
  /// The sequence of edges (v,v_i)
  std::deque<Edge> umbrella;
  /// True iff the umbrella contains all incident edges to the pivot.
  bool _loop;
  /// Used for computations.
  TriangulationHelper<Triangulation, Kernel> TH;

  /// Constructor from triangulation.
  SimplicialStrip( const Triangulation & T ) : TH( T ) {};

  /// The only way to define a strip arount \a pivot, starting at frontier.
  /// The \a predicate should return true for all vertices inside the strip.
  template <typename VertexHandlePredicate>
  void initStrip( const VertexHandlePredicate & predicate, 
                  VertexHandle v, const Edge & frontier )
  {
    ASSERT( TH.source( frontier ) == v );
    umbrella.clear();
    pivot = v;
    _loop = false;
    Edge e = frontier;
    while ( (! _loop) && predicate( TH.target( e ) ) )
      {
        // DGtal::trace.info() << "[initStrip] " << TH.source( e )->point()
        //                     << " -> " << TH.target( e )->point() << std::endl;
        umbrella.push_back( e );
        e = TH.nextCCWAroundSourceVertex( e );
        if ( e == frontier )
          _loop = true; 
        ASSERT( TH.source( e ) == pivot );
      }
    if ( ! _loop )
      {
        umbrella.push_back( e ); // last
        Edge e = TH.nextCWAroundSourceVertex( frontier );
        while ( predicate( TH.target( e ) ) )
          {
            umbrella.push_front( e );
            e = TH.nextCWAroundSourceVertex( e );
            ASSERT( TH.source( e ) == pivot );
          }
        umbrella.push_front( e );
      }
  }
  
  /// @return 'true' iff the strip contains all incident edges to the
  /// pivot (ie the whole umbrella).
  inline bool isLoop() const 
  { return _loop; }

  inline bool isTrivial() const
  {
    ASSERT( umbrella.size() >= 2 );
    return umbrella.size() == 2;
  }
  /// \a i-th vertex of umbrella.
  inline VertexHandle v( unsigned int i ) const
  {
    return TH.target( umbrella.at( i ) );
  }
  /// \a i-th edge of umbrella.
  inline Edge e( unsigned int i ) const
  {
    return umbrella.at( i );
  }
  /// first vertex of umbrella, equivalent to v( 0 )
  inline VertexHandle v0() const
  {
    return TH.target( umbrella.front() );
  }
  /// last vertex of umbrella, equivalent to v( size()-1 )
  inline VertexHandle vn() const
  {
    return TH.target( umbrella.back() );
  }
  /// Number of edges/vertices in umbrella (at least 2, which are then equal).
  inline unsigned int size() const
  { 
    return umbrella.size();
  }
  /// A priority value for the strip, the smaller, the highest.
  inline double priority() const
  {
    return angle(); // ( vn()->point() - v0()->point() ).squared_length();
  }

  /// Angle of umbrella.
  inline Component angle() const
  {
    if ( isTrivial() || isLoop() ) return 2.0*M_PI;
    Component totalAngle = 0.0;
    for ( typename std::deque<Edge>::const_iterator it = umbrella.begin(), ite = umbrella.end() - 1;
          it != ite; ++it )
      {
        totalAngle += TH.innerAngle( *it );
      }
    return totalAngle;
  }
  /// Concavity test.
  bool isConcave() const
  {
    if ( angle() >= M_PI - EPSILON ) return false;
    bool quadrilaterals = true;
    for ( unsigned int i = 1; i < size() - 1; ++i )
      {
        if ( TH.T().is_constrained( e( i ) )
             || ( ! TH.isStabbedCCW( e( i ), v0(), vn() ) ) )
          {
            quadrilaterals = false;
            break;
          }
      }
    return quadrilaterals;
  }

  // @return the index of an edge that is flippable in the strip.
  unsigned int getFlippableEdgeIndex() const
  {
    ASSERT( size() > 2 );
    for ( unsigned int i = 1; i < size() - 1; ++i )
      if ( ( TH.isStabbedCCW( e( i ), 
                              v( i-1 ), 
                              v( i+1 ) ) ) 
           && ( ! TH.T().is_constrained( e( i ) ) ) )
        return i;
    DGtal::trace.error() << "[SimplicialStrip::getFlippableEdge] Unable to find a valid edge." << std::endl;
    return size();
  } 
};

enum GeometryTag 
  { Unknown = -2, Convex = 1, Flat = 0, Concave = -1, Multiple = 2 
  };

template <typename Kernel>
class DAC
{
public:
  typedef typename CGAL::Triangulation_vertex_base_2<Kernel>       Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel>      Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
  typedef CGAL::Exact_predicates_tag                               Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> Triangulation;
  // typedef typename CGAL::Delaunay_triangulation_2<Kernel> Delaunay;
  typedef typename Triangulation::Vertex_circulator       VertexCirculator;
  typedef typename Triangulation::Vertex_handle           VertexHandle;
  typedef typename Triangulation::Edge_iterator           EdgeIterator;
  typedef typename Triangulation::Finite_edges_iterator   FiniteEdgesIterator;
  typedef typename Triangulation::Edge                    Edge;
  typedef typename Triangulation::Finite_faces_iterator   FacesIterator;
  typedef typename Triangulation::Face_handle             FaceHandle;
  typedef typename Triangulation::Point                   Point;
  typedef typename Kernel::Vector_2                  Vector;
  typedef int                                        Label;
  typedef typename std::map<VertexHandle,Label>      VertexLabeling;
  typedef SimplicialStrip<Triangulation,Kernel>      Strip;
  typedef typename std::map<VertexHandle,GeometryTag> VertexGeometryTagging;
  typedef typename std::map<Edge,GeometryTag>         EdgeGeometryTagging;
  typedef typename std::set<Edge>                     EdgeSet;

private:
  /// The current triangulation.
  Triangulation _T;
  /// A mapping VertexHandle -> Label that stores for each vertex to which set it belongs.
  VertexLabeling _vLabeling;
  /// The set of all points of the triangulation (Debug).
  std::set<Point> _points;

public:
  /// Provides useful methods on triangulation
  TriangulationHelper<Triangulation,Kernel> TH;
  static const Label INVALID = -1;

  struct Simplex : public std::vector<VertexHandle>
  {
    Simplex()
    {}
    Simplex( VertexHandle v0 )
    {
      this->push_back( v0 );
    }
    Simplex( VertexHandle v0, VertexHandle v1 )
    {
      this->push_back( v0 );
      this->push_back( v1 );
    }
    Simplex( VertexHandle v0, VertexHandle v1, VertexHandle v2 )
    {
      this->push_back( v0 );
      this->push_back( v1 );
      this->push_back( v2 );
    }
    unsigned int dim() const
    {
      return this->size();
    }
  };

  /// A predicate that returns 'true' whenever the labeling is not the
  /// one given at instanciation.
  struct CheckVertexLabelingInequality {
    const VertexLabeling & _vl;
    Label _l;
    
    inline
    CheckVertexLabelingInequality( const VertexLabeling & vl, Label l )
      : _vl( vl ), _l( l ) 
    {}

    inline
    CheckVertexLabelingInequality( const CheckVertexLabelingInequality & other )
      : _vl( other._vl ), _l( other._l ) 
    {}

    inline
    bool operator()( const VertexHandle & v ) const
    {
      typename VertexLabeling::const_iterator it = _vl.find( v );
      return ( it == _vl.end() ) ? false
        : ( it->second != _l );
    }
  };

  /// A potential concavity.
  struct Concavity {
    VertexHandle _v0;
    VertexHandle _v1;
    const TriangulationHelper<Triangulation,Kernel>* _TH;
    const CheckVertexLabelingInequality* _pred;
    double _n01;

    Concavity( VertexHandle v0, VertexHandle v1, 
               DGtal::ConstAlias< TriangulationHelper<Triangulation,Kernel> > TH,
               DGtal::ConstAlias< CheckVertexLabelingInequality > pred )
      : _v0( v0 ), _v1( v1 ), _TH( TH ), _pred( pred )
    {
      //_n01 = ( _v1->point() - _v0->point() ).squared_length();
      Edge e;
      Strip strip( _TH->T() );
      _n01 = computePriority( e, strip );
    }

    Concavity( VertexHandle v0, VertexHandle v1, 
               DGtal::ConstAlias< TriangulationHelper<Triangulation,Kernel> > TH,
               DGtal::ConstAlias< CheckVertexLabelingInequality > pred,
               double prior )
      : _v0( v0 ), _v1( v1 ), _TH( TH ), _pred( pred ), _n01( prior )
    {}

    Concavity( const Concavity & other )
      : _v0( other._v0 ), _v1( other._v1 ), _n01( other._n01 ),
        _TH( other._TH ), _pred( other._pred )
    {}

    Concavity & operator=( const Concavity & other )
    {
      if ( this != &other )
        {
          _v0 = other._v0;
          _v1 = other._v1;
          _TH = other._TH;
          _pred = other._pred;
          _n01 = other._n01;
        }
      return *this;
    }

    double priority() const
    { 
      return _n01;
    }

    double computePriority( Edge & e, Strip & strip ) const
    {
      bool exist = _TH->findEdge( e, _v0, _v1 );
      strip.initStrip( *_pred, _TH->source( e ), e );
      return strip.priority();
    }

    bool operator<( const Concavity & other ) const
    {
      return _n01 > other._n01;
    }
  };

  //---------------------------------------------------------------------------
public:
  /// Constructor.
  DAC() : TH( _T ) {}

  /// The object is reseted.
  void clear()
  {
    _T.clear();
    _vLabeling.clear();
    _points().clear();
  }

  inline const Triangulation & T() const
  { return _T; }

  inline const VertexLabeling & labeling() const
  { return _vLabeling; }

  inline Label label( VertexHandle v ) const
  {
    typename VertexLabeling::const_iterator it = labeling().find( v );
    if ( it == labeling().end() ) return INVALID;
    return it->second;
  }
  
  /** 
      Adds (digital) points to the triangulation.
      
      @param k the digital topology: 0 is 4-adjacency, 1 is 8.
      @param l the label for the set of points (0 stands for P, 1 for others).
  */
  template <typename PointIterator>
  void add( PointIterator itb, PointIterator ite, int k, Label l )
  {
    std::set<Point> current;
    for ( ; itb != ite; ++itb )
      {
        Point p = *itb;
        if ( _points.find( p ) != _points.end() )
          DGtal::trace.error() << "[DAC::add] Points " << p << " has already been inserted. Ignored." << std::endl;
        else
          {
            VertexHandle v = _T.insert( p );
            _points.insert( p );
            current.insert( p );
            _vLabeling[ v ] = l;
          }
      }
    addConstraints0( current );
    if ( k == 1 ) addConstraints1( current );
  }
  
  void addConstraints0( const std::set<Point> & pts )
  {
    for ( typename std::set<Point>::const_iterator it = pts.begin(), ite = pts.end();
          it != ite; ++it )
      {
        Point p = *it;
        if ( pts.find( p + Vector( 1, 0 ) ) != pts.end() ) _T.insert_constraint( p, p + Vector( 1, 0 ) );
        if ( pts.find( p + Vector( 0, 1 ) ) != pts.end() ) _T.insert_constraint( p, p + Vector( 0, 1 ) );
      }
  }
  
  void addConstraints1( const std::set<Point> & pts )
  {
    for ( typename std::set<Point>::const_iterator it = pts.begin(), ite = pts.end();
          it != ite; ++it )
      {
        Point p = *it;
        if ( ( pts.find( p + Vector( 1, 1 ) ) != pts.end() ) 
             && ( ( pts.find( p + Vector( 1, 0 ) ) == pts.end() ) 
                  || ( pts.find( p + Vector( 0, 1 ) ) == pts.end() ) ) )
          _T.insert_constraint( p, p + Vector( 1, 1 ) );
        if ( ( pts.find( p + Vector( -1, 1 ) ) != pts.end() ) 
             && ( ( pts.find( p + Vector( -1, 0 ) ) == pts.end() ) 
                  || ( pts.find( p + Vector( 0, 1 ) ) == pts.end() ) ) )
          _T.insert_constraint( p, p + Vector( -1, 1 ) );
      }
  }

  bool removeConcavities( Label l ) 
  {
    bool changes = false;
    std::priority_queue<Concavity> Q;
    std::priority_queue<Concavity> Q1;
    CheckVertexLabelingInequality predNotL( labeling(), l );
    DGtal::trace.beginBlock( "Searching concavities" );
    for ( FiniteEdgesIterator it = T().finite_edges_begin(), itend = T().finite_edges_end();
          it != itend; ++it )
      {
        Edge e = *it;
        if ( _T.is_constrained( e ) ) continue;
        VertexHandle v0 = TH.source( e );
        VertexHandle v1 = TH.target( e );
        Label l0 = label( v0 );
        Label l1 = label( v1 );
        if ( ( l0 == l ) && ( l1 != l ) && ( l1 != INVALID ) )
          Q.push( Concavity( v0, v1, TH, predNotL ) );
        else if ( ( l1 == l ) && ( l0 != l ) && ( l0 != INVALID ) )
          Q.push( Concavity( v1, v0, TH, predNotL ) );
      }
    DGtal::trace.info() << "- Found " << Q.size() << " concavities." << std::endl;
    DGtal::trace.endBlock();
    DGtal::trace.beginBlock( "Flipping concavities" );
    while ( ( ! Q.empty() ) || ( ! Q1.empty() ) )
      {
        Concavity concavity = Q1.empty() ? Q.top() : Q1.top();
        if ( Q1.empty() ) Q.pop(); 
        else Q1.pop();
        Edge e;
        Strip strip( _T );
        double p = concavity.computePriority( e, strip );
        if ( p <= 0.0 ) continue;           // Invalid concavity.
        if ( ( p != concavity.priority() )  // priority has changed
             && ( ! Q.empty() )             // but the queue is not empty
             && ( p > Q.top().priority() ) )// and the priority is not the best
          { 
            Q.push( Concavity( concavity._v0, concavity._v1, TH, predNotL, p ) );
            continue;
          }
        // DGtal::trace.info() << "- Edge found " << TH.source( e )->point()
        //                     << "(" << label( TH.source( e ) ) << ")"
        //                     << " -> " << TH.target( e )->point() 
        //                     << "(" << label( TH.target( e ) ) << ")"
        //                     << std::endl;
        // DGtal::trace.info() << "Strip size = " << strip.size() 
        //                     << " angle(deg) = " << ( strip.angle() * 180.0 / M_PI ) << std::endl;
        if ( strip.isConcave() )
          {
            unsigned int idx = strip.getFlippableEdgeIndex();
            Edge fedge = strip.e( idx );
            // DGtal::trace.info() << " => strip is concave. Flipping "
            //                     << TH.source( fedge )->point()
            //                     << "(" << label( TH.source( fedge ) ) << ")"
            //                     << " -> " << TH.target( fedge )->point() 
            //                     << "(" << label( TH.target( fedge ) ) << ")"<< std::endl;
            for ( unsigned int i = 1; i < strip.size() - 1; ++i )
              if ( i != idx ) Q.push( Concavity( strip.pivot, strip.v( i ), TH, predNotL, strip.priority() ) );
            _T.flip( fedge.first, fedge.second );
            changes = true;
          }
      }
    DGtal::trace.endBlock();
    return changes;
  }

  /**
     Returns the boundary with label \a l within the digital affine
     complex.
     
     @param[out] bEdges returns the set of edges of label that lies on
     the boundary. Each edge is a pair <Face,i> such that the \a i-th
     vertex of \a face has a label different from \a l, while the
     other vertices have label \a l.

     @param[in] l the label of the boundary.
  */
  void getBoundarySurface( EdgeSet & bEdges, Label l ) const
  {
    bEdges.clear();
    for ( FacesIterator it = T().finite_faces_begin(), itend = T().finite_faces_end();
          it != itend; ++it )
      {
        FaceHandle f = it;
        int n = 0;
        int j = 0;
        for ( int i = 0; i < 3; ++i )
          {
            if ( label( f->vertex( i ) ) == l )
              ++n;
            else j = i;
          }
        if ( n == 2 )
          { // Boundary edge.
            bEdges.insert( Edge( f, j ) );
          }
      }
  }

  /**
     Return in \a bEdges all edges (s,t) of the triangulation such that label(s)=l0, label(t)!=l0.
     @param[in] l0 any label
     @param[out] the set of such edges.
  */
  void getFrontierEdges( EdgeSet & bEdges, Label l0 ) const
  {
    bEdges.clear();
    for ( FiniteEdgesIterator it = T().finite_edges_begin(), itend = T().finite_edges_end();
          it != itend; ++it )
      {
        Edge e = *it;
        if ( ( label( TH.source( e ) ) == l0 ) 
             && ( label( TH.target( e ) ) != l0 ) )
          bEdges.insert( e );
        else if ( ( label( TH.source( e ) ) != l0 ) 
                  && ( label( TH.target( e ) ) == l0 ) )
          bEdges.insert( T().mirror_edge( e ) );
      }
  }

  /**
     Tags every source of given edges as convex(1),flat(0), concave(-1).
  */
  template <typename EdgeIterator>
  void tagFrontierEdges( VertexGeometryTagging & tag, EdgeIterator it, EdgeIterator ite ) const
  {
    Strip strip( _T );
    for ( ; it != ite; ++it )
      {
        Edge e = *it;
        VertexHandle v = TH.source( e );
        CheckVertexLabelingInequality predNotL( labeling(), label( v ) );
        strip.initStrip( predNotL, v, e );
        double angle = strip.angle();
        if ( angle > M_PI + EPSILON )      tag[ v ] = Convex;
        else if ( angle < M_PI - EPSILON ) tag[ v ] = Concave;
        else                               tag[ v ] = Flat;
      }
  }

  /**
     Tags every edge as convex(1),flat(0), concave(-1).
  */
  template <typename EdgeIterator>
  void tagFrontierEdges( EdgeGeometryTagging & tag, EdgeIterator it, EdgeIterator ite ) const
  {
    for ( ; it != ite; ++it )
      tag[ *it ] = tagFrontierEdge( *it );
  }

  GeometryTag tagFrontierEdge( Edge e ) const
  {
    Strip strip( _T );
    VertexHandle v = TH.source( e );
    CheckVertexLabelingInequality predNotL( labeling(), label( v ) );
    strip.initStrip( predNotL, v, e );
    double angle = strip.angle();
    GeometryTag tag;
    if ( angle > M_PI + EPSILON )      tag = Convex;
    else if ( angle < M_PI - EPSILON ) tag = Concave;
    else                               tag = Flat;
    return tag;
  }
};

/**
   The frontier of a region in a DAC.
*/
template <typename DigitalAffineComplex>
class Frontier
{
public:
  typedef typename DigitalAffineComplex::Triangulation Triangulation;
  typedef typename DigitalAffineComplex::Label         Label;
  typedef typename DigitalAffineComplex::VertexHandle  VertexHandle;
  typedef typename DigitalAffineComplex::Edge          Edge;
  typedef typename DigitalAffineComplex::EdgeIterator  EdgeIterator;
  typedef typename DigitalAffineComplex::Point         Point;
  typedef typename DigitalAffineComplex::VertexGeometryTagging  VertexGeometryTagging;
  typedef typename DigitalAffineComplex::EdgeGeometryTagging    EdgeGeometryTagging;
  typedef typename DigitalAffineComplex::Strip         Strip;
  typedef typename DigitalAffineComplex::Simplex       Simplex;
  typedef typename DigitalAffineComplex::EdgeSet       EdgeSet;
  typedef typename EdgeSet::const_iterator             ConstIterator;


private: 
  const DigitalAffineComplex* _dac;
  Label _l;
  EdgeSet _frontier;
  EdgeGeometryTagging _eTags;

public:
  Frontier( DGtal::ConstAlias<DigitalAffineComplex> dac, Label l )
    : _dac( dac ), _l( l )
  {
    _dac->getFrontierEdges( _frontier, _l );
    _dac->tagFrontierEdges( _eTags, _frontier.begin(), _frontier.end() );
  }

  ConstIterator begin() const
  { return _frontier.begin(); }

  ConstIterator end() const
  { return _frontier.end(); }

  // Circulates to the next frontier edge.
  Edge next( const Edge & e ) const
  {
    // Edge is (s,t) with l(s) = _l, l(st) != _l
    if ( _dac->label( e.first->vertex( e.second ) ) == _l )
      // Turn around t.
      return _dac->T().mirror_edge( _dac->TH.nextCCWAroundFace( e ) );
    else
      // Turn around s
      return _dac->T().mirror_edge( _dac->TH.nextCWAroundFace( e ) );
  }

  void getAllFrontierContours( std::vector< std::vector<VertexHandle> > & contours ) const
  {
    EdgeSet visited;
    for ( ConstIterator it = begin(), ite = end(); it != ite; ++it )
      {
        if ( visited.find( *it ) == visited.end() )
          { // If this edge is not visited yet.
            contours.push_back( std::vector<VertexHandle>() );
            contours.push_back( std::vector<VertexHandle>() );
            getFrontierContour( contours[ contours.size() - 2 ], contours[ contours.size() - 1 ], visited, *it );
          }
      }
  }

  void getAllMLPContours( std::vector< std::vector<VertexHandle> > & contours ) const
  {
    EdgeSet visited;
    for ( ConstIterator it = begin(), ite = end(); it != ite; ++it )
      {
        if ( visited.find( *it ) == visited.end() )
          { // If this edge is not visited yet.
            std::deque<VertexHandle> C;
            getMLPContour( C, visited, *it );
            contours.push_back( std::vector<VertexHandle>( C.size() ) );
            std::copy( C.begin(), C.end(), contours.back().begin() );
          }
      }
  }

  void getAllFrontierSimplices( std::vector< std::vector<Simplex> > & complexes ) const
  {
    EdgeSet visited;
    for ( ConstIterator it = begin(), ite = end(); it != ite; ++it )
      {
        if ( visited.find( *it ) == visited.end() )
          { // If this edge is not visited yet.
            complexes.push_back( std::vector<Simplex>() );
            getFrontierSimplices( complexes.back(), visited, *it );
          }
      }
  }

  void getFrontierSimplices( std::vector<Simplex> & simplices, EdgeSet & fEdges, 
                             const Edge & start ) const
  {
    Edge e = start;
    GeometryTag tag_next;
    VertexHandle v0, v1, v2;
    do {
      fEdges.insert( e );
      Edge next_e = next( e );
      GeometryTag tag_e  = _dac->tagFrontierEdge( e );
      GeometryTag tag_me = _dac->tagFrontierEdge( _dac->T().mirror_edge( e ) );
      v0 = _dac->TH.source( e );
      v1 = _dac->TH.target( e );
      if ( _dac->label( v0 ) != _l )
        DGtal::trace.error() << "[Frontier::getFrontierSimplices] Invalid label for v0." << std::endl;
      if ( _dac->label( v1 ) == _l )
        DGtal::trace.error() << "[Frontier::getFrontierSimplices] Invalid label for v1." << std::endl;
      if ( v0 == _dac->TH.source( next_e ) )
        {
          tag_next = _dac->tagFrontierEdge( _dac->T().mirror_edge( next_e ) );
          v2 = _dac->TH.target( next_e );
        }
      else
        {
          tag_next = _dac->tagFrontierEdge( next_e );
          v2 = _dac->TH.source( next_e );
        }
      if ( ( tag_e == Convex ) || ( tag_e == Flat ) )
        {
          if ( ( tag_me == Convex ) || ( tag_me == Flat ) )
            {
              if ( ( tag_next == Convex ) || ( tag_next == Flat ) )
                simplices.push_back( Simplex( v0, v1, v2 ) );
              else simplices.push_back( Simplex( v0, v1 ) );
            }
          else
            {
              if ( ( tag_next == Convex ) || ( tag_next == Flat ) )
                simplices.push_back( Simplex( v0, v2 ) );
              else simplices.push_back( Simplex( v0 ) );
            }
        }
      else
        {
          if ( ( tag_me == Convex ) || ( tag_me == Flat ) )
            {
              if ( ( tag_next == Convex ) || ( tag_next == Flat ) )
                simplices.push_back( Simplex( v1, v2 ) );
              else simplices.push_back( Simplex( v1 ) );
            }
          else
            {
              if ( ( tag_next == Convex ) || ( tag_next == Flat ) )
                simplices.push_back( Simplex( v2 ) );
              else
                {
                  DGtal::trace.error() << "[Frontier::getFrontierSimplices] Empty simplex." << std::endl;
                  simplices.push_back( Simplex() );
                }
            }
        }
      // Go to next
      e = next_e;
      if ( v0 == _dac->TH.source( next_e ) )
        tag_me = tag_next;
      else
        tag_e = tag_next;
    } while ( e != start );
  }

  /**
     Adds a point with label _l to the current MLP.
     @pre Q0.front() == Q1.front() (it is called a corner).
     @post Q0.front() == Q1.front() (it is called a corner).
  */
  bool addLPoint( std::deque<VertexHandle> & C,
                  std::deque<VertexHandle> & Q0,
                  std::deque<VertexHandle> & Q1,
                  VertexHandle v0 ) const
  {
    bool looped = false;
    ASSERT( ! Q0.empty() && ! Q1.empty() && ( Q0.front() == Q1.front() ) );
    ASSERT( _dac->label( v0 ) == _l );
    if ( v0 == Q0.back() ) return false; // the vertex must be a new vertex.
    while ( ( Q0.size() >= 2 )
            && ( _dac->TH.det( Q0.back()->point() - Q0.front()->point(),
                               v0->point() - Q0.front()->point() ) < 0 ) ) // non convex
      Q0.pop_back();
    Q0.push_back( v0 );
    if ( Q0.size() == 2 )
      {
        while ( ( Q1.size() >= 2 )
                && ( _dac->TH.det( Q1[ 1 ]->point() - Q1.front()->point(), 
                                   Q0.back()->point() - Q1.front()->point() ) < 0 ) ) // non convex
          {
            Q1.pop_front();
            C.push_back( Q1.front() );
            if ( ( C.size() >= 2 ) && ( C.front() == C.back() ) )
              looped = true;
          }
        Q0.front() = Q1.front();
      }
    return looped;
  }
  
  /**
     Adds a point with label _l to the current MLP.
     @pre Q0.front() == Q1.front() (it is called a corner).
     @post Q0.front() == Q1.front() (it is called a corner).
  */
  bool addNotLPoint( std::deque<VertexHandle> & C,
                     std::deque<VertexHandle> & Q0,
                     std::deque<VertexHandle> & Q1,
                     VertexHandle v1 ) const
  {
    bool looped = false;
    ASSERT( ! Q0.empty() && ! Q1.empty() && ( Q0.front() == Q1.front() ) );
    ASSERT( _dac->label( v1 ) != _l );
    if ( v1 == Q1.back() ) return false; // the vertex must be a new vertex.
    while ( ( Q1.size() >= 2 )
            && ( _dac->TH.det( Q1.back()->point() - Q1.front()->point(),
                               v1->point() - Q1.front()->point() ) > 0 ) ) // non concave
      Q1.pop_back();
    Q1.push_back( v1 );
    if ( Q1.size() == 2 )
      {
        while ( ( Q0.size() >= 2 )
                && ( _dac->TH.det( Q0[ 1 ]->point() - Q0.front()->point(),
                                   Q1.back()->point() - Q0.front()->point() ) > 0 ) ) // non concave
          {
            Q0.pop_front();
            C.push_back( Q0.front() );
            if ( ( C.size() >= 2 ) && ( C.front() == C.back() ) )
              looped = true;
          }
        Q1.front() = Q0.front();
      }
    return looped;
  }
  
  void getMLPContour( std::deque<VertexHandle> & C,
                      EdgeSet & fEdges, 
                      const Edge & start ) const
  {
    std::deque<VertexHandle> Q0, Q1;
    // Find first corner.
    VertexHandle v0, v1;
    GeometryTag tag_v0, tag_v1, tag_v2;
    Edge e, next_e;

    DGtal::trace.beginBlock( "Computing MLP of frontier." );
    // Finding a first reasonnable corner.
    e = start;
    do {
      fEdges.insert( e );
      next_e = next( e );
      tag_v0 = _dac->tagFrontierEdge( e );
      tag_v1 = _dac->tagFrontierEdge( _dac->T().mirror_edge( e ) );
      v0 = _dac->TH.source( e );
      v1 = _dac->TH.target( e );
      if ( _dac->TH.source( next_e ) == v0 )
        { // v2 is the target of next_e
          tag_v2 = _dac->tagFrontierEdge( _dac->T().mirror_edge( next_e ) );
          if ( ( tag_v1 == Convex ) && ( tag_v2 == Convex ) )
            {
              Q0.push_back( v1 );
              Q1.push_back( v1 );
              break;
            }
        }
      else
        { // v2 is the source of next_e
          tag_v2 = _dac->tagFrontierEdge( next_e );
          if ( ( tag_v0 == Convex ) && ( tag_v2 == Convex ) )
            {
              Q0.push_back( v0 );
              Q1.push_back( v0 );
              break;
            }
        }
      e = next_e;
    } while ( e != start );

    // If none was found, the MLP is reduced to one point.
    if ( Q0.empty() )
      {
        if ( tag_v0 == Convex ) C.push_back( v0 );
        else if ( tag_v1 == Convex ) C.push_back( v1 );
        else 
          DGtal::trace.error() << "[Frontier::getMLPContour] Unable to find first corner." << std::endl;
        DGtal::trace.endBlock();
        return;
      }

    // Extracting the whole MLP with the invariant that Q0.front() == Q1.front() is a corner.
    e = next_e;
    Edge new_start = e;
    int nbloop = 0;
    do {
      fEdges.insert( e );
      Edge next_e = next( e );
      tag_v0  = _dac->tagFrontierEdge( e );
      tag_v1 = _dac->tagFrontierEdge( _dac->T().mirror_edge( e ) );
      v0 = _dac->TH.source( e );
      ASSERT( _dac->label( v0 ) == _l && "[Frontier::getMLPContour] Invalid label for v0." );
      if ( addLPoint( C, Q0, Q1, v0 ) ) break; // the contour has looped
      v1 = _dac->TH.target( e );
      ASSERT( _dac->label( v1 ) != _l && "[Frontier::getMLPContour] Invalid label for v1." );
      if ( addNotLPoint( C, Q0, Q1, v1 ) ) break; // the contour has looped
      // Go to next
      e = next_e;
      // Both following lines are for debug.
      if ( e == new_start ) ++nbloop; 
      if ( nbloop == 2 ) break;
    } while ( true ); //( e != new_start );
    if ( nbloop == 2 )
      {
        DGtal::trace.error() << "[getMLPContour] nbloop=2 ";
        for ( unsigned int i = 0; i < C.size(); ++i )
          DGtal::trace.info() << C[ i ]->point();
        DGtal::trace.info() << std::endl;
      }
    else {
      while ( C.front() != C.back() )
        C.pop_front();
      C.pop_front();
    }
    DGtal::trace.info() << "[getMLPContour] #contour=" << C.size() << std::endl;
    DGtal::trace.endBlock();
  }
  
  void getFrontierContour( std::vector<VertexHandle> & Q0, std::vector<VertexHandle> & Q1,
                           EdgeSet & fEdges, 
                           const Edge & start ) const
  {
    //std::vector<VertexHandle> Q0, Q1;
    Edge e = start;
    GeometryTag tag_next;
    VertexHandle v0, v1, v2;
    DGtal::trace.beginBlock( "getFrontierContour" );
    do {
      fEdges.insert( e );
      Edge next_e = next( e );
      GeometryTag tag_e  = _dac->tagFrontierEdge( e );
      GeometryTag tag_me = _dac->tagFrontierEdge( _dac->T().mirror_edge( e ) );
      // only convex vertices belong to the MLP.
      v0 = _dac->TH.source( e );
      //DGtal::trace.info() << "[getFrontierContour] Convex  v0=" << v0->point() << std::endl;
      if ( tag_e == Convex )
        {
          ASSERT( _dac->label( v0 ) == _l && "[Frontier::getFrontierContour] Invalid label for v0." );
          if ( Q0.empty() ) 
            Q0.push_back( v0 );
          else if ( v0 != Q0.back() )
            {
              while ( ( Q0.size() >= 2 )                     // at least 2 vertices
                      && ( ( ( _dac->label( Q0.back() ) != _l )  // previous vertex has not the same label
                             && ( _dac->TH.det( Q0.back()->point() - Q0[ Q0.size() - 2 ]->point(),
                                                v0->point() - Q0.back()->point() ) > 0 ) ) // non convex
                           || ( ( _dac->label( Q0.back() ) == _l )  // previous vertex has not the same label
                                && ( _dac->TH.det( Q0.back()->point() - Q0[ Q0.size() - 2 ]->point(),
                                                   v0->point() - Q0.back()->point() ) < 0 ) ) // non convex
                           || ( ( _dac->label( Q0.back() ) != _l )  // previous vertex has not the same label
                                && ( _dac->TH.det( Q0.back()->point() - Q0[ Q0.size() - 2 ]->point(),
                                                   v0->point() - Q0.back()->point() ) == 0 ) // flat
                                && ( ( Q0.back()->point() - Q0[ Q0.size() - 2 ]->point() )
                                     * ( v0->point() - Q0.back()->point() ) < 0.0 ) )
                           ) )
                Q0.pop_back();
              Q0.push_back( v0 );
              while ( ( Q1.size() >= 2 )                     // at least 2 vertices
                      && ( _dac->label( Q1.back() ) == _l )  // previous vertex has not the same label
                      && ( _dac->TH.det( Q1.back()->point() - Q1[ Q1.size() - 2 ]->point(),
                                         v0->point() - Q1.back()->point() ) < 0 ) ) // non concave
                Q1.pop_back();
              Q1.push_back( v0 );
            }
        }
      v1 = _dac->TH.target( e );
      if ( tag_me == Convex )
        {
          //DGtal::trace.info() << "[getFrontierContour] Concave v1=" << v1->point() << std::endl;
          ASSERT( _dac->label( v1 ) != _l && "[Frontier::getFrontierContour] Invalid label for v1." );
          if ( Q1.empty() )
            Q1.push_back( v1 );
          else if ( v1 != Q1.back() )
            {
              while ( ( Q0.size() >= 2 )                     // at least 2 vertices
                      && ( ( ( _dac->label( Q0.back() ) != _l )  // previous vertex has not the same label
                             && ( _dac->TH.det( Q0.back()->point() - Q0[ Q0.size() - 2 ]->point(),
                                                v1->point() - Q0.back()->point() ) > 0 ) ) // non convex
                           || ( ( _dac->label( Q0.back() ) == _l )  // previous vertex has not the same label
                                && ( _dac->TH.det( Q0.back()->point() - Q0[ Q0.size() - 2 ]->point(),
                                                   v1->point() - Q0.back()->point() ) < 0 ) ) // non convex
                           || ( ( _dac->label( Q0.back() ) != _l )  // previous vertex has not the same label
                                && ( _dac->TH.det( Q0.back()->point() - Q0[ Q0.size() - 2 ]->point(),
                                                   v1->point() - Q0.back()->point() ) == 0 ) // flat
                                && ( ( Q0.back()->point() - Q0[ Q0.size() - 2 ]->point() )
                                     * ( v1->point() - Q0.back()->point() ) < 0.0 ) )
                           ) )
                Q0.pop_back();
              Q0.push_back( v1 );
              while ( ( Q1.size() >= 2 )                     // at least 2 vertices
                      && ( _dac->label( Q1.back() ) == _l )  // previous vertex has not the same label
                      && ( _dac->TH.det( Q1.back()->point() - Q1[ Q1.size() - 2 ]->point(),
                                         v1->point() - Q1.back()->point() ) <= 0 ) ) // non concave
                Q1.pop_back();
              Q1.push_back( v1 );
            }
        }
      
      // Go to next
      e = next_e;
    } while ( e != start );
    DGtal::trace.info() << "[getFrontierContour] #Q0=" << Q0.size() << std::endl;
    DGtal::trace.info() << "[getFrontierContour] #Q1=" << Q1.size() << std::endl;
    DGtal::trace.endBlock();
  }
  
  void getContour( std::vector<VertexHandle> & contour, const std::vector<Simplex> & complex ) const
  {
    if ( complex.empty() ) return;
    unsigned int j = complex.size()-1;
    unsigned int k = 0; 
    bool loop = false;
    while ( true )
      {
        if ( complex[ k ].dim() <= 2 )
          {
            if ( complex[ j ].dim() != 3 )
              contour.push_back( complex[ k ].front() );
          }
        j = k++;
        if ( k == complex.size() )
          {
            loop = true;
            k = 0;
          }
        if ( loop && ( contour.front() == contour.back() ) )
          break;
      }
  }
  
};
  
  

template <typename DigitalAffineComplex>
class ViewerDAC
{
public:
  typedef typename DigitalAffineComplex::Triangulation Triangulation;
  typedef typename DigitalAffineComplex::Label         Label;
  typedef typename DigitalAffineComplex::VertexHandle  VertexHandle;
  typedef typename DigitalAffineComplex::Edge          Edge;
  typedef typename DigitalAffineComplex::EdgeIterator  EdgeIterator;
  typedef typename DigitalAffineComplex::Point         CPoint;
  typedef typename DigitalAffineComplex::VertexGeometryTagging  VertexGeometryTagging;
  typedef typename DigitalAffineComplex::EdgeSet       EdgeSet;
  typedef typename DigitalAffineComplex::Simplex       Simplex;
  typedef DGtal::Board2D Board;
  typedef DGtal::Color Color;
  typedef DGtal::Z2i::Point DPoint;

  DPoint toDGtal( const CPoint & p )
  {
    return DPoint( p.x(),p.y() );
  }

private:
  Board & _board;
  const DigitalAffineComplex & _dac;

public:
  ViewerDAC( Board & board, const DigitalAffineComplex & dac )
    : _board( board ), _dac( dac ) 
  {}

  void viewAll()
  {
    Color colors[ 4 ];
    colors[ 0 ] = Color::Red;
    colors[ 1 ] = Color::Green;
    colors[ 2 ] = Color::Blue;
    _board.setLineStyle( Board::Shape::SolidStyle );
    for ( EdgeIterator it = _dac.T().edges_begin(), itend = _dac.T().edges_end();
          it != itend; ++it )
      {
        Edge edge = *it;
        DPoint a = toDGtal( _dac.TH.source( edge )->point() );
        DPoint b = toDGtal( _dac.TH.target( edge )->point() );
        int l1 = _dac.labeling().find( _dac.TH.source( edge ) )->second;
        int l2 = _dac.labeling().find( _dac.TH.target( edge ) )->second;
        int i = ( l1 == 0 ) && ( l2 == 0 ) ? 0
          : ( l1 == 1 ) && ( l2 == 1 ) ? 1
          : 2;
        double w = _dac.T().is_constrained( edge ) ? 2.0 : 1.0;
        _board.setPenColor( colors[ i ] );
        _board.setFillColor( DGtal::Color::None );
        _board.setLineWidth( w );
        _board.drawLine(a[0],a[1],b[0],b[1]);
      }
  }
  
  template <typename EdgeIterator>
  void viewEdges( EdgeIterator itb, EdgeIterator ite, Color c, double w )
  {
    for ( ; itb != ite; ++itb )
      {
        Edge edge = *itb;
        DPoint a = toDGtal( _dac.TH.source( edge )->point() );
        DPoint b = toDGtal( _dac.TH.target( edge )->point() );
        _board.setPenColor( c );
        _board.setFillColor( DGtal::Color::None );
        _board.setLineWidth( w );
        _board.drawLine(a[0],a[1],b[0],b[1]);
      }
  }
  
  void viewBoundarySurface( Label l, Color c, double w )
  {
    EdgeSet edges;
    _dac.getBoundarySurface( edges, l );
    viewEdges( edges.begin(), edges.end(), c, w );
  }

  void viewConvexSurface( Label l, Color c, double w )
  {
    EdgeSet edges;
    VertexGeometryTagging vTags;
    _dac.getFrontierEdges( edges, l );
    _dac.tagFrontierEdges( vTags, edges.begin(), edges.end() );
    _dac.getBoundarySurface( edges, l );
    std::vector<Edge> nonConcaveEdges;
    for ( typename EdgeSet::const_iterator it = edges.begin(), ite = edges.end();
          it != ite; ++it )
      {
        if ( ( vTags[ _dac.TH.source( *it ) ] >= 0 )
             && ( vTags[ _dac.TH.target( *it ) ] >= 0 ) )
          nonConcaveEdges.push_back( *it );
      }
    viewEdges( nonConcaveEdges.begin(), nonConcaveEdges.end(), c, w );
  }


  void viewFrontierContour( Label l, Color c, double w )
  {
    typedef typename Frontier<DigitalAffineComplex>::ConstIterator FrontierConstIterator;
    typedef typename Frontier<DigitalAffineComplex>::ConstIterator FrontierConstIterator;
    Frontier<DigitalAffineComplex> frontier( _dac, l );
    std::vector< std::vector<VertexHandle> > contours;
    frontier.getAllFrontierContours( contours );
    for ( unsigned int i = 0; i < contours.size(); ++i )
      {
        if ( contours[ i ].empty() ) continue;
        std::vector<VertexHandle> & contour = contours[ i ];
        DPoint a = toDGtal( contour[ 0 ]->point() );
        _board.setPenColor( ( i % 2 == 0 ) ? c : DGtal::Color::Magenta  );
        _board.setFillColor( DGtal::Color::None );
        _board.setLineWidth( ( i % 2 == 0 ) ? w : w/2.0 );
        for ( typename std::vector<VertexHandle>::const_iterator its = contour.begin(), itsend = contour.end();
              its != itsend; ++its )
          {
            DPoint b = toDGtal( (*its)->point() );
            _board.drawDot(a[0],a[1]);
            if ( a != b ) _board.drawArrow(a[0],a[1],b[0],b[1]);
            a = b;
          }
      }
  }

  void viewMLPContour( Label l, Color c, double w )
  {
    typedef typename Frontier<DigitalAffineComplex>::ConstIterator FrontierConstIterator;
    typedef typename Frontier<DigitalAffineComplex>::ConstIterator FrontierConstIterator;
    Frontier<DigitalAffineComplex> frontier( _dac, l );
    std::vector< std::vector<VertexHandle> > contours;
    frontier.getAllMLPContours( contours );
    for ( unsigned int i = 0; i < contours.size(); ++i )
      {
        if ( contours[ i ].empty() ) continue;
        std::vector<VertexHandle> & contour = contours[ i ];
        DGtal::HueShadeColorMap<unsigned int,1> hueShade( 0, contour.size() );
        DPoint a = toDGtal( contour[ contour.size()-1 ]->point() );
        _board.setFillColor( DGtal::Color::None );
        _board.setLineWidth( w );
        unsigned int j = 0;
        for ( typename std::vector<VertexHandle>::const_iterator its = contour.begin(), itsend = contour.end();
              its != itsend; ++its, ++j )
          {
            DPoint b = toDGtal( (*its)->point() );
            _board.setPenColor( c ); //hueShade( j ) );
            _board.drawDot(a[0],a[1]);
            if ( a != b ) _board.drawArrow(a[0],a[1],b[0],b[1]);
            a = b;
          }
      }
  }

  void viewFrontier( Label l, Color c, double w )
  {
    typedef typename Frontier<DigitalAffineComplex>::ConstIterator FrontierConstIterator;
    typedef typename Frontier<DigitalAffineComplex>::ConstIterator FrontierConstIterator;
    Frontier<DigitalAffineComplex> frontier( _dac, l );
    std::vector< std::vector<Simplex> > complexes;
    frontier.getAllFrontierSimplices( complexes );
    for ( unsigned int i = 0; i < complexes.size(); ++i )
      {
        std::vector<Simplex> & simplices = complexes[ i ];
        for ( typename std::vector<Simplex>::const_iterator its = simplices.begin(), itsend = simplices.end();
              its != itsend; ++its )
          {
            Simplex s = *its;
            _board.setPenColor( c );
            _board.setFillColor( c );
            _board.setLineWidth( w );
            if ( s.dim() == 1 )
              {
                DPoint a = toDGtal( s[ 0 ]->point() );
                _board.drawDot(a[0],a[1]);
              }
            else if ( s.dim() == 2 )
              {
                DPoint a = toDGtal( s[ 0 ]->point() );
                DPoint b = toDGtal( s[ 1 ]->point() );
                _board.drawLine(a[0],a[1],b[0],b[1]);
              }
            else if ( s.dim() == 3 )
              {
                DPoint a = toDGtal( s[ 0 ]->point() );
                DPoint b = toDGtal( s[ 1 ]->point() );
                DPoint c = toDGtal( s[ 2 ]->point() );
                _board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
              }
          }
      }
  }

  
};


int main( int argc, char** argv )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef DAC<K> DigitalAffineComplex;
  typedef DigitalAffineComplex::Point Point;
  typedef DGtal::Z2i::Space Space;
  typedef DGtal::Z2i::Domain Domain;
  typedef DGtal::Z2i::Point DPoint;
  using namespace DGtal;

  std::set<DPoint> p,q;
  DigitalAffineComplex dac;
  // p.insert( DPoint( 0,0 ) );
  // p.insert( DPoint( 1,1 ) );
  // p.insert( DPoint( 2,1 ) );
  // p.insert( DPoint( 3,1 ) );
  // p.insert( DPoint( 4,0 ) );
  // p.insert( DPoint( 3,-1 ) );
  // p.insert( DPoint( 2,-1 ) );
  // p.insert( DPoint( 1,-1 ) );
  // p.insert( DPoint( 0,-1 ) );
  // p.insert( DPoint( 10,10 ) );
  // p.insert( DPoint( 11,11 ) );
  // p.insert( DPoint( 10,12 ) );

  trace.beginBlock("Construction the shape");
  // typedef Ellipse2D<Z2i::Space> Ellipse; 
  // int a = 5, b = 1;
  // Ellipse2D<Z2i::Space> ellipse(Z2i::Point(0,0), a, b, 0.3 );
  // double h = 0.06125; // 06125; 
  // int N = 16;
  // GaussDigitizer<Z2i::Space,Ellipse> dig;  
  // dig.attach( ellipse );
  // dig.init( ellipse.getLowerBound()+Z2i::Vector(-1,-1),
  //           ellipse.getUpperBound()+Z2i::Vector(1,1), h ); 
  typedef Flower2D<Z2i::Space> Flower; 
  int a = 19, b = 9;
  Flower2D<Z2i::Space> flower(Z2i::Point(0,0), 15, 2, 5, 0);
  double h = 0.5; 
  int N = 5;
  GaussDigitizer<Z2i::Space,Flower> dig;  
  dig.attach( flower );
  dig.init( flower.getLowerBound()+Z2i::Vector(-1,-1),
            flower.getUpperBound()+Z2i::Vector(1,1), h ); 
  Z2i::KSpace ks;
  ks.init( dig.getLowerBound(), dig.getUpperBound(), true );
  SurfelAdjacency<2> sAdj( true );
  Z2i::SCell bel = Surfaces<Z2i::KSpace>::findABel( ks, dig, 1000 );
  std::vector<Z2i::Point> boundaryPoints;
  Surfaces<Z2i::KSpace>
    ::track2DBoundaryPoints( boundaryPoints, ks, sAdj, dig, bel );
  Z2i::Curve c;
  c.initFromVector( boundaryPoints );  
  typedef Z2i::Curve::PointsRange Range; 
  Range r = c.getPointsRange(); 
  for ( Range::ConstIterator it = r.begin(), itE = r.end(); it != itE; ++it )
    {
      Z2i::Point P( *it );
      p.insert( P );
      p.insert( P + Z2i::Point( 2*N+2, -2*N ) );
      p.insert( Z2i::Point( P[ 1 ], P[ 0 ] ) );
    }
  for ( int x = 0; x < 10*N; ++x )
    {
      p.insert( Z2i::Point( x, (12*x)/13 - 4*N - 2 ) );
    }
  trace.endBlock();



  trace.beginBlock( "Computing border." );
  computeBorder<Space>( q, p );
  trace.endBlock();

  trace.beginBlock( "Creating constrained Delaunay triangulation." );
  std::vector<Point> pts;
  for ( std::set<DPoint>::const_iterator it = p.begin(), ite = p.end();
        it != ite; ++it )
    pts.push_back( toCGAL<Point>( *it ) );
  dac.add( pts.begin(), pts.end(), 1, 0 );
  pts.clear();
  for ( std::set<DPoint>::const_iterator it = q.begin(), ite = q.end();
        it != ite; ++it )
    pts.push_back( toCGAL<Point>( *it ) );
  dac.add( pts.begin(), pts.end(), 0, 1 );
  trace.endBlock();

  trace.beginBlock( "Visualizing Delaunay complex." );
  Board2D board;
  ViewerDAC<DigitalAffineComplex> viewer( board, dac );
  for ( std::set<DPoint>::const_iterator it = p.begin(), ite = p.end();
        it != ite; ++it )
    board << *it;
  // board << Domain( DPoint( -2, -2 ), DPoint( 13, 13 ) );
  viewer.viewAll();
  board.saveSVG("delaunay.svg");
  board.saveEPS("delaunay.eps");
  board.clear();
  trace.endBlock();

  bool rc;
  do {
    rc = dac.removeConcavities( 0 );
  } while ( rc );
  do {
    rc = dac.removeConcavities( 1 );
  } while ( rc );
  //   rc = dac.removeConcavities( 0 );
  // rc = dac.removeConcavities( 1 );
  // rc = dac.removeConcavities( 1 );
  //rc = dac.removeConcavities( 0 );
    

  trace.beginBlock( "Visualizing Digital Affine Complex." );
  // Board2D board;
  // ViewerDAC<DigitalAffineComplex> viewer( board, dac );
  for ( std::set<DPoint>::const_iterator it = p.begin(), ite = p.end();
        it != ite; ++it )
    board << *it;
  // board << Domain( DPoint( -2, -2 ), DPoint( 13, 13 ) );
  viewer.viewAll();
  board.saveSVG("dac.svg");
  board.saveEPS("dac.eps");
  board.clear();
  for ( std::set<DPoint>::const_iterator it = p.begin(), ite = p.end();
        it != ite; ++it )
    board << *it;
  // board << Domain( DPoint( -2, -2 ), DPoint( 13, 13 ) );
  std::set<DigitalAffineComplex::Edge> edges;
  viewer.viewConvexSurface( 0, DGtal::Color( 255, 0, 0 ), 4.0 );
  viewer.viewConvexSurface( 1, DGtal::Color( 0, 255, 0 ), 4.0 );
  viewer.viewBoundarySurface( 0, DGtal::Color( 200, 0, 0 ), 1.0 );
  viewer.viewBoundarySurface( 1, DGtal::Color( 0, 200, 0 ), 1.0 );
  // Frontier 0 and Frontier 1 are the same.
  //viewer.viewFrontier( 0, DGtal::Color( 255, 255, 0, 0 ), 3.0 );
  //viewer.viewFrontierContour( 0, DGtal::Color( 0, 255, 0, 0 ), 3.0 );
  viewer.viewMLPContour( 0, DGtal::Color( 0, 0, 255, 0 ), 2.0 );
  board.saveSVG("dac-bdy.svg");
  board.saveEPS("dac-bdy.eps");
  trace.endBlock();
  return 0;
}
