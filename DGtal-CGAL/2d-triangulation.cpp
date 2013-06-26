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
#include "DGtal/io/readers/GenericReader.h"
#include <DGtal/geometry/curves/ArithmeticalDSS.h>
#include <DGtal/geometry/curves/SaturatedSegmentation.h>
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/geometry/curves/FreemanChain.h"



#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

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

/**
   This class represents a sequence of edge-adjacent faces around a
   given edge (\a s, \a t). The user must also provide a predicate \a
   P. Each face of the strip must be incident to vertex \a
   s. Furthermore, inside edges (\a s, \a t_i) of the strip satisfy
   \f$ P(t_i) \f$.

   @code
   Triangulation T;
   ...
   Edge e;
   SimplicialStrip<Triangulation,K> strip( T );
   strip.init( P, e ); // P is a predicate VertexHandle -> bool.
   std::cout << strip.angle() << std::endl; // angle is the sum of the inner angle of each face.
   @endcode

   @tparam Triangulation any kind of CGAL 2D triangulation.
   @tparam Kernel any kind of CGAL 2D Kernel.
*/
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
  /// The sequence of edges (v,v_i)
  std::deque<Edge> umbrella;
  /// True iff the umbrella contains all incident edges to the pivot.
  bool _loop;
  /// Used for computations.
  TriangulationHelper<Triangulation, Kernel> TH;

  /// Constructor from triangulation.
  SimplicialStrip( const Triangulation & T ) : TH( T ) {};

  /// Define a strip around the source edge of \a border, starting around \a border.
  /// The \a predicate should return true for all vertices inside the strip.
  template <typename VertexHandlePredicate>
  void init( const VertexHandlePredicate & predicate, 
             const Edge & border )
  {
    umbrella.clear();
    VertexHandle pivot = TH.source( border );
    _loop = false;
    Edge e = border;
    while ( (! _loop) && predicate( TH.target( e ) ) )
      {
        umbrella.push_back( e );
        e = TH.nextCCWAroundSourceVertex( e );
        if ( e == border )
          _loop = true; 
        ASSERT( TH.source( e ) == pivot );
      }
    if ( ! _loop )
      {
        umbrella.push_back( e ); // last
        Edge e = TH.nextCWAroundSourceVertex( border );
        while ( predicate( TH.target( e ) ) )
          {
            umbrella.push_front( e );
            e = TH.nextCWAroundSourceVertex( e );
            ASSERT( TH.source( e ) == pivot );
          }
        umbrella.push_front( e );
      }
  }

  /// Define a strip around the source edge of \a border, starting
  /// around \a border.  The \a predicate should return true for all
  /// vertices inside the strip.  It also checks if faces have been
  /// extended in the Digital Affine Complex \a dac.
  template <typename DAC, typename VertexHandlePredicate>
  void init( const DAC & dac, 
	     const VertexHandlePredicate & predicate, 
             const Edge & border )
  {
    umbrella.clear();
    VertexHandle pivot = TH.source( border );
    _loop = false;
    Edge e = border;
    while ( (! _loop) && predicate( TH.target( e ) )
	    && ( ! dac.isFaceExtended( dac.T().mirror_edge( e ).first ) ) )
      {
	e = TH.nextCWAroundSourceVertex( e );
	if ( e == border )
	  _loop = true; 
	ASSERT( TH.source( e ) == pivot );
      }
    if ( _loop ) return; // umbrella is the whole neighborhood.
    while ( ! dac.isFaceExtended( e.first ) ) // face has already been extended.
      {
	umbrella.push_back( e );
	e = TH.nextCCWAroundSourceVertex( e );
	ASSERT( TH.source( e ) == pivot );
	if ( ! predicate( TH.target( e ) ) ) break;
      }
    umbrella.push_back( e );
    // while ( (! _loop) && predicate( TH.target( e ) )
    // 	    && ( ! dac.isFaceExtended( e.first ) ) ) // face has already been extended.
    //   {
    //     umbrella.push_back( e );
    //     e = TH.nextCCWAroundSourceVertex( e );
    //     if ( e == border )
    //       _loop = true; 
    //     ASSERT( TH.source( e ) == pivot );
    //   }
    // if ( ! _loop )
    //   {
    //     umbrella.push_back( e ); // last
    //     Edge e = TH.nextCWAroundSourceVertex( border );
    //     while ( predicate( TH.target( e ) ) 
    // 		&& ( ! dac.isFaceExtended
    // 		     ( dac.T().mirror_edge( e ).first ) ) ) // face has already been extended.
    //       {
    // 	    umbrella.push_front( e );
    //         e = TH.nextCWAroundSourceVertex( e );
    //         ASSERT( TH.source( e ) == pivot );
    //       }
    //     umbrella.push_front( e );
    //   }
    // std::cout << "[Umbrella] v=" << pivot->point() << " s=" << umbrella.size()<< std::endl;
    ASSERT( ! predicate( pivot ) );
    ASSERT( predicate( TH.target( border ) ) );
  }

  /// @return 'true' iff the strip was initialized. It must have an
  /// umbrella of size at least 2.
  inline bool isValid() const
  { return ( umbrella.size() >= 1 ) || isLoop(); }

  /// @return 'true' iff the strip contains all incident edges to the
  /// pivot (ie the whole umbrella).
  inline bool isLoop() const 
  { return _loop; }

  /// @return 'true' iff the umbrella is reduced to the initial edge.
  inline bool isTrivial() const
  {
    ASSERT( isValid() );
    return umbrella.size() == 1;
  }

  /// @return 'true' iff the strip was initialized. It must have an
  /// umbrella of size at least 2.
  inline VertexHandle pivot() const
  { 
    ASSERT( isValid() );
    return TH.source( umbrella.front() );
  }

  /**
     @param[in] an index between 0 and umbrella.size()-1.
     @return the \a i-th vertex of umbrella.
  */
  inline VertexHandle v( unsigned int i ) const
  {
    ASSERT( isValid() && ( i < umbrella.size() ) );
    return TH.target( umbrella.at( i ) );
  }
  /**
     @param[in] an index between 0 and umbrella.size()-1.
     @return the i-th edge of umbrella.
  */
  inline Edge e( unsigned int i ) const
  {
    ASSERT( isValid() && ( i < umbrella.size() ) );
    return umbrella.at( i );
  }
  /// first vertex of umbrella, equivalent to v( 0 )
  inline VertexHandle v0() const
  {
    ASSERT( isValid() );
    return TH.target( umbrella.front() );
  }
  /// last vertex of umbrella, equivalent to v( size()-1 )
  inline VertexHandle vn() const
  {
    ASSERT( isValid() );
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
    ASSERT( isValid() );
    if ( isLoop() ) return 2.0*M_PI;
    if ( isTrivial() ) return 2.0*M_PI;
    // if ( ( size() == 2 ) && ( umbrella.front() == umbrella.back() ) )
    //   return 2.0*M_PI;
    Component totalAngle = 0.0;
    for ( typename std::deque<Edge>::const_iterator it = umbrella.begin(), ite = umbrella.end() - 1;
          it != ite; ++it )
      {
        totalAngle += TH.innerAngle( *it );
      }
    return totalAngle;
  }
  /// Concavity test.
  bool isConcaveAndFlippable() const
  {
    ASSERT( isValid() );
    if ( isTrivial() ) return false;
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
  /// Concavity test.
  bool isConcave() const
  {
    ASSERT( isValid() );
    if ( isTrivial() ) return false;
    return ( angle() < M_PI - EPSILON );
  }
  /// Concavity test.
  bool isConvex() const
  {
    ASSERT( isValid() );
    if ( isTrivial() ) return false;
    return ( angle() > M_PI + EPSILON );
  }
  /// Concavity test.
  bool isFlat() const
  {
    ASSERT( isValid() );
    if ( isTrivial() ) return false;
    double a = angle();
    return ( a <= M_PI + EPSILON ) && ( a >= M_PI - EPSILON );
  }
  /// Concavity test.
  bool isFlippable() const
  {
    ASSERT( isValid() );
    if ( isTrivial() ) return false;
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
    for ( unsigned int i = 1; i < size() - 1; ++i )
      if ( ( TH.isStabbedCCW( e( i ), 
                              v( i-1 ), 
                              v( i+1 ) ) ) 
           && ( ! TH.T().is_constrained( e( i ) ) ) )
        return i;
    // DGtal::trace.error() << "[SimplicialStrip::getFlippableEdge] Unable to find a valid edge." << std::endl;
    return size();
  } 

  // @return the index of the concavity in the strip.
  unsigned int getUnflippableEdgeIndex() const
  {
    ASSERT( size() > 2 );
    for ( unsigned int i = 1; i < size() - 1; ++i )
      if ( ! TH.isStabbedCCW( e( i ), 
                              v( i-1 ), 
                              v( i+1 ) ) )
        return i;
    DGtal::trace.error() << "[SimplicialStrip::getUnflippableEdge] Unable to find a valid edge." << std::endl;
    return size();
  } 
};

enum GeometryTag 
  { Unknown = -2, Convex = 1, Flat = 0, Concave = -1, Multiple = 2, Infinite = 3  
  };

/**
   A DAC, short word for Digital Affine Complex is a triangulation of
   digital points with specific properties.  Digital points are
   labelled. Each set of points with the same label is given digital
   topology (0: 4-adjacency, 1: 8-adjacency). The triangulation
   contains the points, edges, faces, induced by each
   adjacency. Furthermore, the triangulation satisfies the relative
   convex hull property, i.e. forall P,Q points with same label l,
   [PQ] does not intersect a point, edge, face of another region
   implies that [PQ] belongs to simplices whose vertices have label l.

   @tparam Kernel it is the chosen Kernel for the CGAL triangulation,
   for instance CGAL::Exact_predicates_inexact_constructions_kernel
*/
template <typename Kernel>
class DAC
{
public:
  typedef typename CGAL::Triangulation_vertex_base_2<Kernel>       Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel>      Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
  typedef CGAL::No_intersection_tag                                Itag;
  // Exact_predicates_tag                               Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> Triangulation;
  // typedef typename CGAL::Delaunay_triangulation_2<Kernel> Delaunay;
  typedef typename Triangulation::Finite_vertices_iterator FiniteVerticesIterator;
  typedef typename Triangulation::Vertex_circulator       VertexCirculator;
  typedef typename Triangulation::Vertex_handle           VertexHandle;
  typedef typename Triangulation::Edge_iterator           EdgeIterator;
  typedef typename Triangulation::Finite_edges_iterator   FiniteEdgesIterator;
  typedef typename Triangulation::All_edges_iterator      AllEdgesIterator;
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
  // maps inside edges to outside edges. Each edge has either no or one outside edge.
  typedef typename std::map<Edge,Edge>                In2OutEdgeMapping;
  // maps outside edges to inside edges. Each edge has either no or two inside edges.
  typedef typename std::map<Edge,std::pair<Edge,Edge> > Out2inEdgeMapping;
  

  /**
     This structure is associated with a face. When a face is added to
     the relative hull, some edges of the face are inside the hull,
     some edges are outside the hull. We can therefore build an
     extension from the inside edges to the outside edges, as well as
     a retraction from the outside edges and the inner face to the
     inside edge(s).
   */
  struct Extension {
    FaceHandle face;
    int inside_edges[ 2 ]; //< stores the inside edges (indices in the face)
    int outside_edges[ 2 ];//< stores the outside edges (indices in the face)
    bool type;             //< 'true' 1 outside, 2 inside edges,
			   //'false' 2 outside, 1 inside.
    Extension() {}

    Extension( const Extension & other )
      : face( other.face ), type( other.type )
    {
      inside_edges[ 0 ] = inside_edges[ 1 ] = -1;
      outside_edges[ 0 ] = outside_edges[ 1 ] = -1;
    }

    Extension( FaceHandle f, int outside )
    {
      face = f;
      outside_edges[ 0 ] = outside;
      outside_edges[ 1 ] = -1;
      inside_edges[ 0 ] = ( outside + 1 ) % 3;
      inside_edges[ 1 ] = ( outside + 2 ) % 3;
      type = true;
    }

    Extension & operator=( const Extension & other )
    {
      if ( this != &other )
	{
	  face = other.face;
	  inside_edges[ 0 ] = other.inside_edges[ 0 ];
	  inside_edges[ 1 ] = other.inside_edges[ 1 ];
	  outside_edges[ 0 ] = other.outside_edges[ 0 ];
	  outside_edges[ 1 ] = other.outside_edges[ 1 ];
	}
      return *this;
    }

    Extension( FaceHandle f, int outside1, int outside2 )
    {
      face = f;
      outside_edges[ 0 ] = ( outside1+1 ) % 3 == outside2 ? outside1 : outside2;
      outside_edges[ 1 ] = ( outside1+1 ) % 3 == outside2 ? outside2 : outside1;
      inside_edges[ 0 ] = ( outside_edges[ 1 ] + 1 ) % 3;
      inside_edges[ 1 ] = -1;
      type = false;
    }
    
    Point retract( const Triangulation & T, const Point & p ) const
    {
      Point p1 = T.segment( Edge( face, inside_edges[ 0 ] ) )
	.supporting_line().projection( p );
      if ( ! type ) return p1;
      
    }
  };

  /**
     A simplex of dimension n is the convex hull of n+1 points.
  */ 
  struct Simplex : public std::vector<VertexHandle>
  {
    typedef typename std::vector<VertexHandle>::iterator Iterator;
    typedef typename std::vector<VertexHandle>::const_iterator ConstIterator;

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
      std::sort( this->begin(), this->end() );
    }
    Simplex( VertexHandle v0, VertexHandle v1, VertexHandle v2 )
    {
      this->push_back( v0 );
      this->push_back( v1 );
      this->push_back( v2 );
      std::sort( this->begin(), this->end() );
    }
    Simplex( VertexHandle v0, VertexHandle v1, VertexHandle v2, VertexHandle v3 )
    {
      this->push_back( v0 );
      this->push_back( v1 );
      this->push_back( v2 );
      this->push_back( v3 );
      std::sort( this->begin(), this->end() );
    }
    int dim() const
    {
      return ((int)this->size())-1;
    }
    bool operator<( const Simplex & other ) const
    {
      ASSERT( this->dim() == other.dim() );
      for ( ConstIterator it1 = this->begin(), it2 = other.begin(), it1e = this->end();
            ( it1 != it1e ); ++it1, ++it2 )
        if ( *it1 < *it2 ) return true;
        else if ( *it2 < *it1 ) return false;
      return false;
    }
    bool operator==( const Simplex & other ) const
    {
      ASSERT( this->dim() == other.dim() );
      ConstIterator it1 = this->begin(), it2 = other.begin(), it1e = this->end();
      while ( ( it1 != it1e ) && ( *it1 == *it2 ) )
        { ++it1, ++it2; }
      return it1 == it1e;
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

  /// A predicate that returns 'false' when the vertex is marked, or
  /// if not, returns 'true' whenever the labeling is not the one
  /// given at instanciation, .
  struct CheckVertexLabelingInequalityWhenUnmarked {
    const VertexLabeling & _vl;
    const VertexLabeling & _vmark;
    Label _l;
    Label _mark;
    
    inline
    CheckVertexLabelingInequalityWhenUnmarked
    ( const VertexLabeling & vl, const VertexLabeling & vmark,
      Label l, Label mark )
      : _vl( vl ), _vmark( vmark ), _l( l ), _mark( mark )
    {}

    inline
    CheckVertexLabelingInequalityWhenUnmarked
    ( const CheckVertexLabelingInequalityWhenUnmarked & other )
      : _vl( other._vl ), _vmark( other._vmark ), _l( other._l ), _mark( other._mark )
    {}

    inline
    bool operator()( const VertexHandle & v ) const
    {
      typename VertexLabeling::const_iterator it = _vmark.find( v );
      if ( it != _vmark.end() ) return it->second != _mark;
      it = _vl.find( v );
      return ( it == _vl.end() ) ? false
        : ( it->second != _l );
    }
  };


  /// A potential concavity.
  template <typename Predicate>
  struct GenericConcavity {
    VertexHandle _v0;
    VertexHandle _v1;
    const TriangulationHelper<Triangulation,Kernel>* _TH;
    const Predicate* _pred;
    double _n01;

    GenericConcavity( VertexHandle v0, VertexHandle v1, 
                      DGtal::ConstAlias< TriangulationHelper<Triangulation,Kernel> > TH,
                      DGtal::ConstAlias< Predicate > pred )
      : _v0( v0 ), _v1( v1 ), _TH( TH ), _pred( pred )
    {
      //_n01 = ( _v1->point() - _v0->point() ).squared_length();
      Edge e;
      Strip strip( _TH->T() );
      _n01 = computePriority( e, strip );
      ASSERT( ! (*_pred)( _v0 ) );
      ASSERT( (*_pred)( _v1 ) );
    }

    GenericConcavity( VertexHandle v0, VertexHandle v1, 
                      DGtal::ConstAlias< TriangulationHelper<Triangulation,Kernel> > TH,
                      DGtal::ConstAlias< Predicate > pred,
                      double prior )
      : _v0( v0 ), _v1( v1 ), _TH( TH ), _pred( pred ), _n01( prior )
    {}

    GenericConcavity( const GenericConcavity & other )
      : _v0( other._v0 ), _v1( other._v1 ), _n01( other._n01 ),
        _TH( other._TH ), _pred( other._pred )
    {}

    GenericConcavity & operator=( const GenericConcavity & other )
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
      strip.init( *_pred, e );
      return strip.priority();
    }

    bool operator<( const GenericConcavity & other ) const
    {
      return _n01 > other._n01;
    }
  };


private:
  /// The current triangulation.
  Triangulation _T;
  /// A mapping VertexHandle -> Label that stores for each vertex to which set it belongs.
  VertexLabeling _vLabeling;
  /// A mapping VertexHandle -> Label that stores for each vertex if it is marked.
  VertexLabeling _vMarked;

  typedef typename std::map<FaceHandle, Extension> FaceMapping;
  typedef typename FaceMapping::iterator FaceMappingIterator;
  typedef typename FaceMapping::const_iterator FaceMappingConstIterator;
  /// A mapping Face -> Info that stores faces that have been added to
  /// the relative hull for each vertex if it is marked.
  FaceMapping _vFaces;
  /// The set of all points of the triangulation (Debug).
  std::set<Point> _points;
  /// The set of quads that have already been flipped.
  std::set<Simplex> _flippedQuads;

  DGtal::Board2D _debugBoard;

public:
  /// Provides useful methods on triangulation
  TriangulationHelper<Triangulation,Kernel> TH;
  static const Label INVALID = -1;


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
    _vMarked.clear();
    _vFaces.clear();
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

  inline Label mark( VertexHandle v ) const
  {
    typename VertexLabeling::const_iterator it = _vMarked.find( v );
    if ( it == _vMarked.end() ) return -1;
    return it->second;
  }

  /**
     @param f any face of the triangulation
     @return 'true' iff the face has been extended during the relative
     hull process and is part of global retraction of the relative
     hull onto the canonic complex.
  */
  inline bool isFaceExtended( FaceHandle f ) const
  {
    FaceMappingConstIterator it = _vFaces.find( f );
    return it != _vFaces.end();
  }
  
  FaceMappingConstIterator itFaceExtension( FaceHandle f ) const
  {
    return _vFaces.find( f );
  }
  FaceMappingConstIterator endFaceExtension() const
  {
    return _vFaces.end();
  }
  /**
     @param f a face of the triangulation that has been extended during the relative
     hull process.
     @return a const_reference to the associated extension.
     @pre 'isFaceExtended( f ) == true'
  */
  const Extension & extension( FaceHandle f ) const
  {
    ASSERT( isFaceExtended( f ) );
    return *( itFaceExtension( f ) );
  }

  inline void setInfiniteLabel( Label l )
  {
    _vLabeling[ _T.infinite_vertex() ] = l;
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

  void addConstraint( const Point & p, const Point & q )
  {
    VertexHandle vp = _T.insert( p );
    VertexHandle vq = _T.insert( q );
    _T.insert_constraint( vp, vq );
  }
  // For speed-up
  void addBand( const Point & p1, const Point & p2, const Point & midp, 
                const Point & q1, const Point & q2, const Point & midq,
                bool p1cvx, bool p2cvx )
  {
    DGtal::trace.info() << "[DAC.addBand] ";
    int val = ( p1cvx ? 2 : 0 ) + ( p2cvx ? 1 : 0 );
    switch ( val ) {
    case 0: // ccv1, ccv2
      DGtal::trace.info() << " C=" << q1 << " -> " << midp;
      DGtal::trace.info() << " C=" << midp << " -> " << q2 << std::endl;
      _T.insert_constraint( q1, midp );
      _T.insert_constraint( midp, q2 );
      break;
    case 1: // ccv1, cvx2
      // These constraints may be false afterwards for the relative
      // hull. We remove them.
      // DGtal::trace.info() << " C=" << q1 << " -> " << p2 << std::endl;
      // _T.insert_constraint( q1, p2 );
      break;
    case 2: // cvx1, ccv2
      // These constraints may be false afterwards for the relative
      // hull. We remove them.
      // DGtal::trace.info() << " C=" << p1 << " -> " << q2 << std::endl;
      // _T.insert_constraint( p1, q2 );
      break;
    case 3: // cvx1, cvx2
      DGtal::trace.info() << " C=" << p1 << " -> " << midq;
      DGtal::trace.info() << " C=" << midq << " -> " << p2 << std::endl;
      _T.insert_constraint( p1, midq );
      _T.insert_constraint( midq, p2 );
      break;
    }
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

  /**
     This procedure removes all concavities around vertices with label \a
     l. Complexity is O(e x n), where e is the number of edges.

     @param[in] any label, as given with method \ref add.
     @return 'true' if some concavity was flipped, 'false' when no concavity was found.
  */
  void removeAllConcavities( Label l ) 
  {
    _flippedQuads.clear();
    bool rc;
    unsigned int i = 0;
    do {
      _debugBoard.clear();
      _debugBoard.setPenColor( DGtal::Color::Black );
      _debugBoard.setFillColor( DGtal::Color::None );
      _debugBoard.setLineWidth( 2.0 );
      rc = this->removeConcavities( l );
      std::ostringstream sname;
      sname << "tmp/dac-rc-" << l << "-" << i++ << ".eps";
      _debugBoard.saveEPS( sname.str().c_str() );
    } while ( rc );
  }

  /**
     This procedure removes concavities around vertices with label \a
     l. Complexity is O(e), where e is the number of edges.

     @param[in] any label, as given with method \ref add.
     @return 'true' if some concavity was flipped, 'false' when no concavity was found.
  */
  bool removeConcavities( Label l ) 
  {
    bool changes = false;
    typedef GenericConcavity< CheckVertexLabelingInequality > Concavity;
    std::priority_queue<Concavity> Q;
    std::set< std::pair< VertexHandle, VertexHandle > > inQueue;

    CheckVertexLabelingInequality predNotL( labeling(), l );
    // DGtal::trace.beginBlock( "Searching concavities" );
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
          {
            Concavity c( v0, v1, TH, predNotL );
            if ( c.priority() < M_PI - EPSILON ) 
              { 
                Q.push( c ); // only concave vertices are flippable.
                inQueue.insert( std::make_pair( v0, v1 ) );
              }
          }
        else if ( ( l1 == l ) && ( l0 != l ) && ( l0 != INVALID ) )
          {
            Concavity c( v1, v0, TH, predNotL );
            if ( c.priority() <  M_PI - EPSILON ) 
              {
                Q.push( c ); // only concave vertices are flippable.
                inQueue.insert( std::make_pair( v1, v0 ) );
              }
          }
      }
    DGtal::trace.info() << "- Found " << Q.size() << " potential concavities." << std::endl;
    // DGtal::trace.endBlock();
    // DGtal::trace.beginBlock( "Flipping concavities" );
    unsigned int nb_checked = 0;
    unsigned int nb_flipped = 0;
    while ( ! Q.empty() )
      {
        Concavity concavity = Q.top();
        Q.pop(); ++nb_checked;
        if ( nb_checked % 1000 == 0 ) 
          DGtal::trace.info() << "- Queue=" << Q.size()
                              << ", flipped " << nb_flipped << "/" << nb_checked 
                              << " edges in a concavity." << std::endl;
        inQueue.erase( std::make_pair( concavity._v0, concavity._v1 ) );
        Edge e;
        Strip strip( _T );
        double p = concavity.computePriority( e, strip );
        if ( p <= 0.0 ) continue;           // Invalid concavity.
        if ( ( p != concavity.priority() )  // priority has changed
             && ( ! Q.empty() )             // but the queue is not empty
             && ( p > Q.top().priority() ) )// and the priority is not the best
          { 
            Q.push( Concavity( concavity._v0, concavity._v1, TH, predNotL, p ) );
            inQueue.insert( std::make_pair( concavity._v0, concavity._v1 ) );
            continue;
          }
        if ( strip.isConcaveAndFlippable() )
          {
            unsigned int idx = strip.getFlippableEdgeIndex();
            Edge fedge = strip.e( idx );
            Simplex quad( concavity._v0, strip.v( idx - 1 ), concavity._v1, strip.v( idx + 1) );
            if ( _flippedQuads.find( quad ) != _flippedQuads.end() )
              continue; // already flipped.
            _flippedQuads.insert( quad );
            std::vector<Concavity> smallQ;
            for ( unsigned int i = 1; i < strip.size() - 1; ++i )
              if ( i != idx ) 
                smallQ.push_back( Concavity( strip.pivot(), strip.v( i ), TH, predNotL, strip.priority() ) );
            // if ( idx != 1 ) 
            smallQ.push_back( Concavity( strip.v0(), strip.v( 1 ), TH, predNotL ) );
            // if ( idx != strip.size() - 2 ) 
            smallQ.push_back( Concavity( strip.vn(), strip.v( strip.size() - 2 ), TH, predNotL ) );
            for ( typename std::vector<Concavity>::const_iterator it = smallQ.begin(),
                    ite = smallQ.end(); it != ite; ++it )
              { 
                if ( inQueue.find( std::make_pair( it->_v0, it->_v1 ) ) == inQueue.end() )
                  {
                    Q.push( *it );
                    inQueue.insert( std::make_pair( it->_v0, it->_v1 ) );
                  }
              }
            DGtal::Z2i::Point a = toDGtal( TH.source( fedge )->point() );
            DGtal::Z2i::Point b = toDGtal( TH.target( fedge )->point() );
            _debugBoard.drawLine(a[0],a[1],b[0],b[1]);
            _T.flip( fedge.first, fedge.second );
            changes = true; ++nb_flipped;
          }
      }
    ASSERT( inQueue.empty() );
    ASSERT( Q.empty() );
    DGtal::trace.info() << "- Flipped " << nb_flipped << "/" << nb_checked << " edges in a concavity." << std::endl;
    // DGtal::trace.endBlock();
    return changes;
  }

  void insertQueue( std::set< std::pair< VertexHandle, VertexHandle > > & inQueue,
		    VertexHandle v0, VertexHandle v1,
		    const CheckVertexLabelingInequality & predicate ) const
  {
    if ( ( ! predicate( v0 ) ) && predicate( v1 ) )
      {
	Edge e;
	if ( T().is_edge( v0, v1, e.first, e.second ) )
	  {
	    Strip strip( T() );
	    strip.init( *this, predicate, e );
	    if ( strip.isConcave() ) 
	      inQueue.insert( std::make_pair( v0, v1 ) );
          }
      }
  }

  /**
     This procedure computes the relative hull with the same principle
     as removeConcavities. around vertices with label \a l. Complexity
     is O(e), where e is the number of edges.

     @param[in] any label, as given with method \ref add.
     @return 'true' if some concavity was flipped, 'false' when no concavity was found.
  */
  bool relativeHull2( Label l ) 
  {
    bool changes = false;
    static const Label rmark = 1;
    std::set< std::pair< VertexHandle, VertexHandle > > inQueue;
    CheckVertexLabelingInequality predNotL( labeling(), l );
    Strip strip( T() );
    // DGtal::trace.beginBlock( "Searching concavities" );
    for ( FiniteEdgesIterator it = T().finite_edges_begin(), itend = T().finite_edges_end();
          it != itend; ++it )
      {
        Edge e = *it;
        if ( T().is_constrained( e ) ) continue;
        VertexHandle v0 = TH.source( e );
        VertexHandle v1 = TH.target( e );
        Label l0 = label( v0 );
        Label l1 = label( v1 );
        if ( ( l0 == l ) && predNotL( v1 ) && ( l1 != INVALID ) )
	  // && ( _vMarked.find( v1 ) == _vMarked.end() ) )
          {
	    strip.init( *this, predNotL, e );
	    if ( strip.isConcave() ) 
	      inQueue.insert( std::make_pair( v0, v1 ) );
          }
        else if ( ( l1 == l ) && predNotL( v0 ) && ( l0 != INVALID ) )
	  // && ( _vMarked.find( v0 ) == _vMarked.end() ) )
          {
	    strip.init( *this, predNotL, T().mirror_edge( e ) );
	    if ( strip.isConcave() ) 
	      inQueue.insert( std::make_pair( v1, v0 ) );
          }
      }
    DGtal::trace.info() << "- Found " << inQueue.size() << " potential concavities." 
			<< std::endl;
    // DGtal::trace.endBlock();
    // DGtal::trace.beginBlock( "Flipping concavities" );
    unsigned int nb_checked = 0;
    unsigned int nb_flipped = 0;
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
	if ( T().is_constrained( e ) )
	  {
	    inQueue.erase( inQueue.begin() );
	    continue; // (v0,v1) is constrained
	  }
	strip.init( *this, predNotL, e );
	ASSERT( strip.isValid() );
        if ( strip.isConcave() ) 
          {
	    unsigned int i = strip.getFlippableEdgeIndex();
	    if ( i != strip.size() ) // at least one edge is flippable.
	      {
                Edge fedge = strip.e( i );
		ASSERT( ! T().is_constrained( fedge ) );
		ASSERT( ! isFaceExtended( fedge.first ) );
		Edge mirror_first = T().mirror_edge( strip.e( 0 ) );
		// std::cout << "[relativeHull2]" 
		// 	  << " a=" << strip.pivot()->point()
		// 	  << " b=" << strip.v(i-1)->point()
		// 	  << " c=" << strip.v(i)->point()
		// 	  << " d=" << strip.v(i+1)->point() << std::endl;
		VertexHandle _v0 = strip.v0();
		VertexHandle _v1 = strip.v( 1 );
		VertexHandle _vn = strip.vn();
		VertexHandle _vn1 = strip.v( strip.size() - 2 );
		_T.flip( fedge.first, fedge.second );
		if ( strip.size() == 3 ) 
		  { // last flip made a triangle
		    Edge nedge = T().mirror_edge( mirror_first );
		    VertexHandle sv0 = TH.target( nedge );
		    VertexHandle sv1 = TH.target( Edge( nedge.first, 
							T().ccw( nedge.second ) ) );
		    if ( predNotL( sv0 ) ) _vMarked[ sv0 ] = rmark;
		    if ( predNotL( sv1 ) ) _vMarked[ sv1 ] = rmark;
		    _vFaces[ nedge.first ] = Extension( nedge.first, 
							T().ccw( nedge.second ) );
		    // std::cout << "[relativeHull2] Last flip is a triangle." << std::endl;
		  }
		else
		  {
		    insertQueue( inQueue, _v0, _v1, predNotL );
		    insertQueue( inQueue, _vn, _vn1, predNotL );
		  }
		// else 
		//   std::cout << "[relativeHull2] flip, strip size =" << strip.size() << std::endl;
		changes = true; ++nb_flipped;
	      }
	    else // none are flippable. This is a concave/flat piece.
	      {
		Edge fedge = strip.e( 0 );
		if ( ! isFaceExtended( fedge.first ) )
		  {
		    if ( predNotL( strip.v0() ) ) _vMarked[ strip.v0() ] = rmark;
		    _vFaces[ fedge.first ] = Extension( fedge.first, 
							T().ccw( fedge.second ) );
		    for ( unsigned int i = 1; i < strip.size() - 1; ++i ) 
		      {
			_vMarked[ strip.v( i ) ] = rmark;
			fedge = strip.e( i );
			ASSERT( ! isFaceExtended( fedge.first ) );
			_vFaces[ fedge.first ] = Extension( fedge.first, 
							    fedge.second,
							    T().ccw( fedge.second ) );
			// DGtal::trace.info() << "- mark " << TH.target( fedge )->point()
			// 		    << " source=" << TH.source( fedge )->point() << std::endl;
		      }
		    if ( predNotL( strip.vn() ) ) _vMarked[ strip.vn() ] = rmark;
		    changes = true;
		  }
		inQueue.erase( inQueue.begin() );
	      }
          }
	else
	  inQueue.erase( inQueue.begin() );
      }
    ASSERT( inQueue.empty() );
    DGtal::trace.info() << "- Flipped " << nb_flipped << "/" << nb_checked << " edges in a concavity." << std::endl;
    //if ( nb_flipped == 1 && nb_checked == 2 ) changes = false;
    // DGtal::trace.endBlock();
    return changes;
  }


  /**
     This procedure computes the relative hull with the same principle
     as removeConcavities. around vertices with label \a l. Complexity
     is O(e), where e is the number of edges.

     @param[in] any label, as given with method \ref add.
     @return 'true' if some concavity was flipped, 'false' when no concavity was found.
  */
  bool relativeHull( Label l ) 
  {
    bool changes = false;
    typedef GenericConcavity< CheckVertexLabelingInequalityWhenUnmarked > Concavity;
    std::priority_queue<Concavity> Q;
    std::set< std::pair< VertexHandle, VertexHandle > > inQueue;
    const Label mark = 1;
    CheckVertexLabelingInequalityWhenUnmarked predNotLWU( labeling(), _vMarked, l, mark );
    // DGtal::trace.beginBlock( "Searching concavities" );
    for ( FiniteEdgesIterator it = T().finite_edges_begin(), itend = T().finite_edges_end();
          it != itend; ++it )
      {
        Edge e = *it;
        if ( _T.is_constrained( e ) ) continue;
        VertexHandle v0 = TH.source( e );
        VertexHandle v1 = TH.target( e );
        Label l0 = label( v0 );
        Label l1 = label( v1 );
        if ( ( l0 == l ) && predNotLWU( v1 ) && ( l1 != INVALID ) )
          {
            Concavity c( v0, v1, TH, predNotLWU );
            if ( c.priority() < M_PI - EPSILON ) 
              { 
                Q.push( c ); // only concave vertices are flippable.
                inQueue.insert( std::make_pair( v0, v1 ) );
              }
          }
        else if ( ( l1 == l ) && predNotLWU( v0 ) && ( l0 != INVALID ) )
          {
            Concavity c( v1, v0, TH, predNotLWU );
            if ( c.priority() <  M_PI - EPSILON ) 
              {
                Q.push( c ); // only concave vertices are flippable.
                inQueue.insert( std::make_pair( v1, v0 ) );
              }
          }
      }
    DGtal::trace.info() << "- Found " << Q.size() << " potential concavities." << std::endl;
    // DGtal::trace.endBlock();
    // DGtal::trace.beginBlock( "Flipping concavities" );
    unsigned int nb_checked = 0;
    unsigned int nb_flipped = 0;
    while ( ! Q.empty() )
      {
        Concavity concavity = Q.top();
        Q.pop(); ++nb_checked;
        if ( nb_checked % 1000 == 0 ) 
          DGtal::trace.info() << "- Queue=" << Q.size()
                              << ", flipped " << nb_flipped << "/" << nb_checked 
                              << " edges in a concavity." << std::endl;
        inQueue.erase( std::make_pair( concavity._v0, concavity._v1 ) );
        Edge e;
        Strip strip( _T );
        double p = concavity.computePriority( e, strip );
        if ( p <= 0.0 ) continue;           // Invalid concavity.
        if ( ( p != concavity.priority() )  // priority has changed
             && ( ! Q.empty() )             // but the queue is not empty
             && ( p > Q.top().priority() ) )// and the priority is not the best
          { 
            Q.push( Concavity( concavity._v0, concavity._v1, TH, predNotLWU, p ) );
            inQueue.insert( std::make_pair( concavity._v0, concavity._v1 ) );
            continue;
          }
        if ( strip.isConcave() ) 
          {
            if ( strip.isFlippable() )
              {
                unsigned int idx = strip.getFlippableEdgeIndex();
                Edge fedge = strip.e( idx );
                Simplex quad( concavity._v0, strip.v( idx - 1 ), concavity._v1, strip.v( idx + 1) );
                if ( _flippedQuads.find( quad ) != _flippedQuads.end() )
                  continue; // already flipped.
                _flippedQuads.insert( quad );
                std::vector<Concavity> smallQ;
                for ( unsigned int i = 1; i < strip.size() - 1; ++i )
                  if ( i != idx ) 
                    smallQ.push_back( Concavity( strip.pivot(), strip.v( i ), 
                                                 TH, predNotLWU, strip.priority() ) );
                // if ( idx != 1 ) 
                if ( label( strip.v0() ) == l )
                  smallQ.push_back( Concavity( strip.v0(), strip.v( 1 ), TH, predNotLWU ) );
                // if ( idx != strip.size() - 2 ) 
                if ( label( strip.vn() ) == l )
                  smallQ.push_back( Concavity( strip.vn(), strip.v( strip.size() - 2 ), TH, predNotLWU ) );
                for ( typename std::vector<Concavity>::const_iterator it = smallQ.begin(),
                        ite = smallQ.end(); it != ite; ++it )
                  { 
                    if ( inQueue.find( std::make_pair( it->_v0, it->_v1 ) ) == inQueue.end() )
                      {
                        Q.push( *it );
                        inQueue.insert( std::make_pair( it->_v0, it->_v1 ) );
                      }
                  }
                DGtal::Z2i::Point a = toDGtal( TH.source( fedge )->point() );
                DGtal::Z2i::Point b = toDGtal( TH.target( fedge )->point() );
                _debugBoard.drawLine(a[0],a[1],b[0],b[1]);
                _T.flip( fedge.first, fedge.second );
                changes = true; ++nb_flipped;
              }
            else
              {
                unsigned int idx = strip.getUnflippableEdgeIndex();
                if ( idx != strip.size() )
                  {
                    Edge fedge = strip.e( idx );
                    _vMarked[ TH.target( fedge ) ] = mark;
                    // DGtal::trace.info() << "- mark " << TH.target( fedge )->point()
                    //                     << " source=" << TH.source( fedge )->point() << std::endl;
                    changes = true;
                  }
              }
          }
      }
    ASSERT( inQueue.empty() );
    ASSERT( Q.empty() );
    DGtal::trace.info() << "- Flipped " << nb_flipped << "/" << nb_checked << " edges in a concavity." << std::endl;
    // DGtal::trace.endBlock();
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
  void getBorderEdges( EdgeSet & bEdges, Label l0 ) const
  {
    bEdges.clear();
    // for ( FiniteEdgesIterator it = T().finite_edges_begin(), itend = T().finite_edges_end();
    //       it != itend; ++it )
    for ( AllEdgesIterator it = T().all_edges_begin(), itend = T().all_edges_end();
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
  void tagBorderEdges( VertexGeometryTagging & tag, EdgeIterator it, EdgeIterator ite ) const
  {
    Strip strip( _T );
    for ( ; it != ite; ++it )
      {
        Edge e = *it;
        if ( T().is_infinite( e ) )
          tag[ TH.source( e ) ] = T().is_infinite( TH.source( e ) ) ? Concave : Convex;
        else 
          {
            VertexHandle v = TH.source( e );
            CheckVertexLabelingInequality predNotL( labeling(), label( v ) );
            strip.init( predNotL, e );
            double angle = strip.angle();
            if ( angle > M_PI + EPSILON )      tag[ v ] = Convex;
            else if ( angle < M_PI - EPSILON ) tag[ v ] = Concave;
            else                               tag[ v ] = Flat;
          }
      }
  }

  /**
     Tags every edge as convex(1),flat(0), concave(-1).
  */
  template <typename EdgeIterator>
  void tagBorderEdges( EdgeGeometryTagging & tag, EdgeIterator it, EdgeIterator ite ) const
  {
    for ( ; it != ite; ++it )
      tag[ *it ] = tagBorderEdge( *it );
  }

  GeometryTag tagBorderEdge( Edge e ) const
  {
    GeometryTag tag;
    if ( T().is_infinite( e ) )
      tag = T().is_infinite( TH.source( e ) ) ? Concave : Convex;
    else
      {
        Strip strip( _T );
        VertexHandle v = TH.source( e );
        CheckVertexLabelingInequality predNotL( labeling(), label( v ) );
        strip.init( predNotL, e );
        double angle = strip.angle();
        if ( angle > M_PI + EPSILON )      tag = Convex;
        else if ( angle < M_PI - EPSILON ) tag = Concave;
        else                               tag = Flat;
      }
    return tag;
  }


  bool flipBorder( Label l )
  {
    static const int array_flip[ 16 ] = // 0: no flip, 1: flip always; 2: check length.
      { 0 /* dcba */, 0 /* dcbA */, 1 /* dcBa */, 2 /* dcBA */,
        0 /* dCba */, 0 /* dCbA */, 2 /* dCBa */, 0 /* dCBA */,
        1 /* Dcba */, 2 /* DcbA */, 1 /* DcBa */, 1 /* DcBA */,
        2 /* DCba */, 0 /* DCbA */, 1 /* DCBa */, 2 /* DCBA */ };
      // { 0 /* dcba */, 0 /* dcbA */, 1 /* dcBa */, 2 /* dcBA */,
      //   0 /* dCba */, 0 /* dCbA */, 2 /* dCBa */, 2 /* dCBA */,
      //   1 /* Dcba */, 2 /* DcbA */, 1 /* DcBa */, 2 /* DcBA */,
      //   2 /* DCba */, 2 /* DCbA */, 2 /* DCBa */, 2 /* DCBA */ };
    GeometryTag tag_a, tag_b, tag_c, tag_d;
    EdgeSet bEdges;
    this->getBorderEdges( bEdges, l );
    unsigned int nbflips = 0;
    while ( ! bEdges.empty() )
      { // Get edge
        typename EdgeSet::iterator itEdge = bEdges.begin();
        Edge ac = *itEdge;
        bEdges.erase( itEdge );
        Edge ca = T().mirror_edge( ac );
        Edge ab = TH.nextCWAroundSourceVertex( ac );
        Edge ad = TH.nextCCWAroundSourceVertex( ac );
        if ( T().is_infinite( ac ) || T().is_infinite( ab ) || T().is_infinite( ad ) ) continue;
        // Get vertices E = a -> c
        VertexHandle a = TH.source( ac );
        VertexHandle b = TH.target( ab );
        VertexHandle c = TH.target( ac );
        VertexHandle d = TH.target( ad );
        ASSERT( ( label( a ) == l ) && ( label( c ) != l ) );
        if ( label( b ) == label( d ) ) continue;
        if ( ! TH.isConvexQuadrilateral( a, b, c, d ) ) continue;
        if ( label( b ) == l )
          {
            tag_a = tagBorderEdge( ac );                         // ac
            tag_b = tagBorderEdge( TH.nextCCWAroundFace( ab ) ); // bc
            tag_c = tagBorderEdge( T().mirror_edge( ac ) );      // ca
            tag_d = tagBorderEdge( T().mirror_edge( ad ) );      // da
          }
        else // ( label( d ) == l )
          {
            tag_a = tagBorderEdge( ac );                         // ac
            tag_b = tagBorderEdge( T().mirror_edge( ab ) );      // ba
            tag_c = tagBorderEdge( T().mirror_edge( ac ) );      // ca
            tag_d = tagBorderEdge( T().mirror_edge( TH.nextCCWAroundFace( ac ) ) ); // dc
          }
        int flag = ( ( tag_a != Concave ) ? 1 : 0 ) 
          + ( ( tag_b != Concave ) ? 2 : 0 )
          + ( ( tag_c != Concave ) ? 4 : 0 ) 
          + ( ( tag_d != Concave ) ? 8 : 0 );
        if ( array_flip[ flag ] == 0 ) continue; // no flip
        if ( array_flip[ flag ] == 2 )           // check length
          {
            double l_ac = ( c->point() - a->point() ).squared_length();
            double l_bd = ( d->point() - b->point() ).squared_length();
            if ( l_ac >= l_bd ) continue;
          }
        if ( label( b ) == l )
          {
            bEdges.insert( TH.nextCCWAroundFace( ab ) ); // bc
            bEdges.insert( ad );                         // ad
          }
        else  // ( label( d ) == l )
          {
            bEdges.insert( ab );                         // ab
            bEdges.insert( T().mirror_edge( TH.nextCCWAroundFace( ac ) ) ); // dc
          }
        _T.flip( ac.first, ac.second );
        nbflips++;
      }
    DGtal::trace.info() << "[DAC::flipBorder] nbflips=" << nbflips << std::endl;
    return nbflips != 0;
  }


};

/**
   The border of a region in a DAC.
*/
template <typename DigitalAffineComplex>
class Border
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
  EdgeSet _border;
  EdgeGeometryTagging _eTags;

public:
  Border( DGtal::ConstAlias<DigitalAffineComplex> dac, Label l )
    : _dac( dac ), _l( l )
  {
    _dac->getBorderEdges( _border, _l );
    _dac->tagBorderEdges( _eTags, _border.begin(), _border.end() );
  }

  /// @return an iterator on the first border edge for this label.
  ConstIterator begin() const
  { return _border.begin(); }

  /// @return an iterator after the last border edge for this label.
  ConstIterator end() const
  { return _border.end(); }

  /**
     @param e any edge of the border complex.
     @return the geometric tag for this edge (Convex, Flat, Concave, Unknown).
  */
  GeometryTag tag( const Edge & e ) const
  {
    typename EdgeGeometryTagging::const_iterator it = _eTags.find( e );
    return ( it != _eTags.end() ) ? it->second : Unknown;
  }

  // Circulates to the next border edge.
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

  void getAllBorderContours( std::vector< std::vector<VertexHandle> > & contours ) const
  {
    EdgeSet visited;
    for ( ConstIterator it = begin(), ite = end(); it != ite; ++it )
      {
        if ( visited.find( *it ) == visited.end() )
          { // If this edge is not visited yet.
            contours.push_back( std::vector<VertexHandle>() );
            contours.push_back( std::vector<VertexHandle>() );
            getBorderContour( contours[ contours.size() - 2 ], contours[ contours.size() - 1 ], visited, *it );
          }
      }
  }

  /**
     @param[out] returns the list of simplical complexes that defines the border
     complex of the DAC.
   */
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
      GeometryTag tag_e  = _dac->tagBorderEdge( e );
      GeometryTag tag_me = _dac->tagBorderEdge( _dac->T().mirror_edge( e ) );
      v0 = _dac->TH.source( e );
      v1 = _dac->TH.target( e );
      if ( _dac->label( v0 ) != _l )
        DGtal::trace.error() << "[Border::getFrontierSimplices] Invalid label for v0." << std::endl;
      if ( _dac->label( v1 ) == _l )
        DGtal::trace.error() << "[Border::getFrontierSimplices] Invalid label for v1." << std::endl;
      if ( v0 == _dac->TH.source( next_e ) )
        {
          tag_next = _dac->tagBorderEdge( _dac->T().mirror_edge( next_e ) );
          v2 = _dac->TH.target( next_e );
        }
      else
        {
          tag_next = _dac->tagBorderEdge( next_e );
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
                  DGtal::trace.error() << "[Border::getFrontierSimplices] Empty simplex." << std::endl;
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

     @param[inout] C contains the (potentially updated) sequence of MLP points.

     @param[inout] Q0 contains the (potentially updated) sequence of
     inside points: the first is common with the first of Q1 and is a
     corner, the others are points with label _l.

     @param[inout] Q1 contains the (potentially updated) sequence of
     outside points: the first is common with the first of Q0 and is a
     corner, the others are points with label != _l.

     @param[in] v0 a point with label _l.

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
     Adds a point with label != _l to the current MLP.

     @param[inout] C contains the (potentially updated) sequence of MLP points.

     @param[inout] Q0 contains the (potentially updated) sequence of
     inside points: the first is common with the first of Q1 and is a
     corner, the others are points with label _l.

     @param[inout] Q1 contains the (potentially updated) sequence of
     outside points: the first is common with the first of Q0 and is a
     corner, the others are points with label != _l.

     @param[in] v1 a point with label != _l.

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
  
  /**
     Computes the minimum length polygon (MLP) that remains in the
     border complex, given a starting border edge in the border complex.

     @param[out] C contains the sequence of points of the MLP.
     @param[inout] fEdges is updated with the set of visited border
     edges.  @param[in] start the border edge that determines a
     component of the border complex.
  */
  void getMLPContour( std::deque<VertexHandle> & C,
                      EdgeSet & fEdges, 
                      const Edge & start ) const
  {
    std::deque<VertexHandle> Q0, Q1;
    // Find first corner.
    VertexHandle v0, v1;
    GeometryTag tag_v0, tag_v1, tag_v2;
    Edge e, next_e;

    DGtal::trace.beginBlock( "Computing MLP of border." );
    // Finding a first reasonnable corner.
    e = start;
    do {
      fEdges.insert( e );
      next_e = next( e );
      tag_v0 = tag( e ); //_dac->tagBorderEdge( e );
      tag_v1 = _dac->tagBorderEdge( _dac->T().mirror_edge( e ) );
      v0 = _dac->TH.source( e );
      v1 = _dac->TH.target( e );
      if ( _dac->TH.source( next_e ) == v0 )
        { // v2 is the target of next_e
          tag_v2 = _dac->tagBorderEdge( _dac->T().mirror_edge( next_e ) );
          if ( ( tag_v1 == Convex ) && ( tag_v2 == Convex ) )
            {
              Q0.push_back( v1 );
              Q1.push_back( v1 );
              break;
            }
        }
      else
        { // v2 is the source of next_e
          tag_v2 = tag( next_e ); // _dac->tagBorderEdge( next_e );
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
          DGtal::trace.error() << "[Border::getMLPContour] Unable to find first corner." << std::endl;
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
      // tag_v0  = _dac->tagBorderEdge( e );
      // tag_v1 = _dac->tagBorderEdge( _dac->T().mirror_edge( e ) );
      v0 = _dac->TH.source( e );
      ASSERT( _dac->label( v0 ) == _l && "[Border::getMLPContour] Invalid label for v0." );
      if ( addLPoint( C, Q0, Q1, v0 ) ) break; // the contour has looped
      v1 = _dac->TH.target( e );
      ASSERT( _dac->label( v1 ) != _l && "[Border::getMLPContour] Invalid label for v1." );
      if ( addNotLPoint( C, Q0, Q1, v1 ) ) break; // the contour has looped
      // Go to next
      e = next_e;
      // Both following lines are for debug.
      if ( e == new_start ) ++nbloop; 
      if ( nbloop == 2 ) break;
    } while ( true ); //( e != new_start );
    if ( nbloop == 2 )
      {
        DGtal::trace.info() << "[Border::getMLPContour] Bad initialization, nbloop=2 ";
        for ( unsigned int i = 0; i < C.size(); ++i )
          DGtal::trace.info() << C[ i ]->point();
        DGtal::trace.info() << std::endl;
        DGtal::trace.info() << "[Border::getMLPContour] Correcting contour." << std::endl;
      }
    while ( C.front() != C.back() )
      C.pop_front();
    C.pop_front();
    DGtal::trace.info() << "[Border::getMLPContour] #contour=" << C.size() << std::endl;
    DGtal::trace.endBlock();
  }

  /**
     @return a vector containing polygons as sequence of vertices;
     each polygon is the MLP (minimum length/perimeter polygon) of
     this part of the border complex.

     @see getMLPContour
  */
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


  /**
     Experimental. This is a tentative for extracting two frontier
     contours (like an inside MLP and an outside MLP) from the
     frontier simplicial complex.
  */
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
      GeometryTag tag_e  = _dac->tagBorderEdge( e );
      GeometryTag tag_me = _dac->tagBorderEdge( _dac->T().mirror_edge( e ) );
      // only convex vertices belong to the MLP.
      v0 = _dac->TH.source( e );
      //DGtal::trace.info() << "[getFrontierContour] Convex  v0=" << v0->point() << std::endl;
      if ( tag_e == Convex )
        {
          ASSERT( _dac->label( v0 ) == _l && "[Border::getFrontierContour] Invalid label for v0." );
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
          ASSERT( _dac->label( v1 ) != _l && "[Border::getFrontierContour] Invalid label for v1." );
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
  

  
};
  
  
/**
   This class is intended for visualizing digital affine complex and
   their border complexes.
*/
template <typename DigitalAffineComplex>
class ViewerDAC
{
public:
  typedef typename DigitalAffineComplex::Triangulation   Triangulation;
  typedef typename DigitalAffineComplex::Label           Label;
  typedef typename DigitalAffineComplex::VertexHandle    VertexHandle;
  typedef typename DigitalAffineComplex::FiniteVerticesIterator FiniteVerticesIterator;
  typedef typename DigitalAffineComplex::Edge            Edge;
  typedef typename DigitalAffineComplex::EdgeIterator    EdgeIterator;
  typedef typename DigitalAffineComplex::FacesIterator   FacesIterator;
  typedef typename DigitalAffineComplex::FaceHandle      FaceHandle;
  typedef typename DigitalAffineComplex::Point           CPoint;
  typedef typename DigitalAffineComplex::VertexGeometryTagging  VertexGeometryTagging;
  typedef typename DigitalAffineComplex::EdgeSet         EdgeSet;
  typedef typename DigitalAffineComplex::Simplex         Simplex;
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
  std::vector<Color> _label_colors;
  std::vector<Color> _other_colors;

public:
  /**
     Constructor. Requires a board \a board for display and a digital
     affine complex \a dac.
  */
  ViewerDAC( DGtal::Alias<Board> board, 
             DGtal::ConstAlias<DigitalAffineComplex> dac )
    : _board( board ), _dac( dac ) 
  {
    _label_colors.push_back( Color::Red );
    _label_colors.push_back( Color::Green );
    _label_colors.push_back( Color::Blue );
    _label_colors.push_back( Color::Magenta );
    _other_colors.push_back( Color::Black );
    _other_colors.push_back( Color::Yellow );
  }

  /**
     View vertices of the triangulation. 
     @param l when -1, view all vertices, otherwise view the vertices with given label.
  */
  void viewVertices( int l )
  {
    DPoint dummy;
    std::string specificStyle =  dummy.className() + "/Grid";
    _board << DGtal::SetMode( dummy.className(), "Grid" );
    for ( FiniteVerticesIterator it = _dac.T().finite_vertices_begin(), itend = _dac.T().finite_vertices_end();
          it != itend; ++it )
      {
        int l1 = _dac.label( it );
        if ( ( l == -1 ) || ( l1 == l ) )
          {
            Color c = _label_colors[ l1 ];
            DPoint a = toDGtal( it->point() );
            if ( _dac.mark( it ) != -1 )
              _board << DGtal::CustomStyle( specificStyle, new DGtal::CustomColors( Color::Black, Color::Black ) );
            else
              _board << DGtal::CustomStyle( specificStyle, new DGtal::CustomColors( c, c ) );
            _board << a;
          }
      }
  }

  /**
     View vertices of the triangulation. 
     @param l when -1, view all vertices, otherwise view the vertices with given label.
  */
  void viewPavingVertices( int l )
  {
    DPoint dummy;
    std::string specificStyle =  dummy.className() + "/Paving";
    _board << DGtal::SetMode( dummy.className(), "Paving" );
    
    for ( FiniteVerticesIterator it = _dac.T().finite_vertices_begin(), itend = _dac.T().finite_vertices_end();
          it != itend; ++it )
      {
        int l1 = _dac.label( it );
        if ( ( l == -1 ) || ( l1 == l ) )
          {
            Color c = _label_colors[ l1 ];
            DPoint a = toDGtal( it->point() );
            _board << DGtal::CustomStyle( specificStyle, new DGtal::CustomColors( Color::Black, c ) )
                   << a;
          }
      }
  }

  void viewAllTriangulationEdges()
  {
    _board.setLineStyle( Board::Shape::SolidStyle );
    _board.setFillColor( DGtal::Color::None );
    _board.setLineWidth( 1.0 );
    for ( EdgeIterator it = _dac.T().edges_begin(), itend = _dac.T().edges_end();
          it != itend; ++it )
      {
        Edge edge = *it;
        DPoint a = toDGtal( _dac.TH.source( edge )->point() );
        DPoint b = toDGtal( _dac.TH.target( edge )->point() );
        int l1 = _dac.labeling().find( _dac.TH.source( edge ) )->second;
        int l2 = _dac.labeling().find( _dac.TH.target( edge ) )->second;
        Color c = ( l1 == l2 ) ? _label_colors[ l1 ] : _other_colors[ 0 ]; 
        // double w = 1.0; // _dac.T().is_constrained( edge ) ? 2.0 : 1.0;
        _board.setPenColor( c );
        _board.drawLine(a[0],a[1],b[0],b[1]);
      }
  }

  void viewTriangulationEdges( int l )
  {
    _board.setLineStyle( Board::Shape::SolidStyle );
    _board.setFillColor( DGtal::Color::None );
    _board.setLineWidth( 1.0 );
    for ( EdgeIterator it = _dac.T().edges_begin(), itend = _dac.T().edges_end();
          it != itend; ++it )
      {
        Edge edge = *it;
        DPoint a = toDGtal( _dac.TH.source( edge )->point() );
        DPoint b = toDGtal( _dac.TH.target( edge )->point() );
        int l1 = _dac.labeling().find( _dac.TH.source( edge ) )->second;
        int l2 = _dac.labeling().find( _dac.TH.target( edge ) )->second;
        if ( ( l == -1 ) || ( ( l == l1 ) && ( l == l2 ) ) )
          {
            Color c = ( l1 == l2 ) ? _label_colors[ l1 ] : _other_colors[ 0 ]; 
            // double w = 1.0; // _dac.T().is_constrained( edge ) ? 2.0 : 1.0;
            _board.setPenColor( c );
            _board.drawLine(a[0],a[1],b[0],b[1]);
          }
      }
  }

  void viewConstrainedEdges( int l )
  {
    _board.setLineStyle( Board::Shape::SolidStyle );
    _board.setFillColor( DGtal::Color::None );
    _board.setLineWidth( 3.0 );
    for ( EdgeIterator it = _dac.T().edges_begin(), itend = _dac.T().edges_end();
          it != itend; ++it )
      {
        Edge edge = *it;
        DPoint a = toDGtal( _dac.TH.source( edge )->point() );
        DPoint b = toDGtal( _dac.TH.target( edge )->point() );
        int l1 = _dac.labeling().find( _dac.TH.source( edge ) )->second;
        int l2 = _dac.labeling().find( _dac.TH.target( edge ) )->second;
        if ( ( ( l == -1 ) || ( l == l1 ) )
             && _dac.T().is_constrained( edge ) )
          {
            Color c = ( l1 == l2 ) ? _label_colors[ l1 ] : _other_colors[ 1 ];
            _board.setPenColor( c );
            _board.drawLine(a[0],a[1],b[0],b[1]);
          }
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
  
  void viewBoundarySurface( Label l, double w )
  {
    EdgeSet edges;
    _dac.getBoundarySurface( edges, l );
    viewEdges( edges.begin(), edges.end(), _label_colors[ l ], w );
  }

  void viewConvexSurface( Label l, Color c, double w )
  {
    EdgeSet edges;
    VertexGeometryTagging vTags;
    _dac.getBorderEdges( edges, l );
    _dac.tagBorderEdges( vTags, edges.begin(), edges.end() );
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
    typedef typename Border<DigitalAffineComplex>::ConstIterator BorderConstIterator;
    Border<DigitalAffineComplex> border( _dac, l );
    std::vector< std::vector<VertexHandle> > contours;
    border.getAllFrontierContours( contours );
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

  void viewMLPContour( Label l, Color c, double w, bool arrow )
  {
    typedef typename Border<DigitalAffineComplex>::ConstIterator BorderConstIterator;
    Border<DigitalAffineComplex> border( _dac, l );
    std::vector< std::vector<VertexHandle> > contours;
    border.getAllMLPContours( contours );
    for ( unsigned int i = 0; i < contours.size(); ++i )
      {
        if ( contours[ i ].empty() ) continue;
        std::vector<VertexHandle> & contour = contours[ i ];
        DGtal::GradientColorMap<unsigned int> gradient( 0, contour.size() + 1 );
        gradient.addColor( _label_colors[ l ] );
        gradient.addColor( c );
        gradient.addColor( _label_colors[ l ] );
        DPoint a = toDGtal( contour[ contour.size()-1 ]->point() );
        _board.setFillColor( DGtal::Color::None );
        _board.setLineWidth( w );
        unsigned int j = 0;
        for ( typename std::vector<VertexHandle>::const_iterator its = contour.begin(), itsend = contour.end();
              its != itsend; ++its, ++j )
          {
            DPoint b = toDGtal( (*its)->point() );
            _board.setPenColor( arrow ? gradient( j ) : c );
            _board.drawDot(a[0],a[1]);
            if ( a != b ) 
              {
                if ( arrow ) _board.drawArrow(a[0],a[1],b[0],b[1]);
                else         _board.drawLine(a[0],a[1],b[0],b[1]);
              }
            a = b;
          }
      }
  }

  void viewFrontier( Label l, Color pc, Color c, double w )
  {
    typedef typename Border<DigitalAffineComplex>::ConstIterator BorderConstIterator;
    Border<DigitalAffineComplex> border( _dac, l );
    std::vector< std::vector<Simplex> > complexes;
    border.getAllFrontierSimplices( complexes );
    for ( unsigned int i = 0; i < complexes.size(); ++i )
      {
        std::vector<Simplex> & simplices = complexes[ i ];
        for ( typename std::vector<Simplex>::const_iterator its = simplices.begin(), itsend = simplices.end();
              its != itsend; ++its )
          {
            Simplex s = *its;
            _board.setPenColor( pc );
            _board.setFillColor( c );
            _board.setLineWidth( w );
            if ( s.dim() == 0 )
              {
                DPoint a = toDGtal( s[ 0 ]->point() );
                _board.drawDot(a[0],a[1]);
              }
            else if ( s.dim() == 1 )
              {
                DPoint a = toDGtal( s[ 0 ]->point() );
                DPoint b = toDGtal( s[ 1 ]->point() );
                _board.drawLine(a[0],a[1],b[0],b[1]);
              }
            else if ( s.dim() == 2 )
              {
                DPoint a = toDGtal( s[ 0 ]->point() );
                DPoint b = toDGtal( s[ 1 ]->point() );
                DPoint c = toDGtal( s[ 2 ]->point() );
                _board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
              }
          }
      }
  }

  /**
     The border is the set of simplices that borders the affine complex of region \a l.
  */
  void viewBorder( Label l, Color pc, Color c, double w )
  {
    typedef typename Border<DigitalAffineComplex>::ConstIterator BorderConstIterator;
    Border<DigitalAffineComplex> border( _dac, l );
    _board.setPenColor( pc );
    _board.setFillColor( c );
    _board.setLineWidth( w );
    for ( BorderConstIterator it = border.begin(), ite = border.end();
          it != ite; ++it )
      {
        FaceHandle f = (*it).first;
        DPoint a = toDGtal( f->vertex( 0 )->point() );
        DPoint b = toDGtal( f->vertex( 1 )->point() );
        DPoint c = toDGtal( f->vertex( 2 )->point() );
        _board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
      }
  }  

  void viewRelativeHull( int l )
  {
    _board.setLineStyle( Board::Shape::SolidStyle );
    DGtal::Color fillc = _label_colors[ l ];
    if ( fillc.red() < 128 ) fillc.red( 200 );
    if ( fillc.green() < 128 ) fillc.green( 200 );
    if ( fillc.blue() < 128 ) fillc.blue( 200 );
    _board.setFillColor( fillc );
    _board.setPenColor( fillc );
    _board.setLineWidth( 0.0 );
    for ( FacesIterator it = _dac.T().finite_faces_begin(), itend = _dac.T().finite_faces_end();
          it != itend; ++it )
      {
        FaceHandle f = it;
        VertexHandle v0 = f->vertex( 0 );
        VertexHandle v1 = f->vertex( 1 );
        VertexHandle v2 = f->vertex( 2 );
        // if ( ( ( _dac.label( v0 ) == l ) || ( _dac.mark( v0 ) == 1 ) )
        //      && ( ( _dac.label( v1 ) == l ) || ( _dac.mark( v1 ) == 1 ) )
        //      && ( ( _dac.label( v2 ) == l ) || ( _dac.mark( v2 ) == 1 ) ) )
        if ( _dac.isFaceExtended( f ) ||
	     ( ( _dac.label( v0 ) == l ) && ( _dac.label( v1 ) == l ) 
	       && ( _dac.label( v2 ) == l ) ) )
          {
            DPoint a = toDGtal( v0->point() );
            DPoint b = toDGtal( v1->point() );
            DPoint c = toDGtal( v2->point() );
            _board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
          }
      }
  }

};
  
  
template <typename Image, typename Predicate,
          typename DAC>
void
insertBands( DAC & dac, const Image & image, const Predicate & predicate,
             int k )
{
  typedef typename Image::Domain Domain;
  typedef typename Domain::Space Space;
  typedef typename Domain::ConstSubRange ConstSubRange;
  typedef typename Space::Point DPoint;
  typedef typename DAC::Point Point;
  // Find bands
  Domain imageDomain = image.domain();
  DPoint d1 = DPoint::diagonal( 1 );
  Domain extDomain( imageDomain.lowerBound() - d1, imageDomain.upperBound() + d1 );
  {
    // Horizontal bands
    ConstSubRange subRangeY = extDomain.subRange( 1, extDomain.lowerBound() );
    DPoint band_p, band_p1;
    for ( typename ConstSubRange::ConstIterator ity = subRangeY.begin(), 
            ityend = subRangeY.end(); ity != ityend; ++ity )
      {
        bool prev_p_inside  = false; // value at (x-1, y)
        bool prev_p1_inside = false; // value at (x-1, y+1)
        bool inside = false;
        ConstSubRange subRangeX = extDomain.subRange( 0, (*ity) + DPoint( 1, 0 ) );
        int sBand = 0;
        for ( typename ConstSubRange::ConstIterator itx = subRangeX.begin(), 
                itxend = subRangeX.end(); itx != itxend; ++itx )
          {
            DPoint p = *itx;
            DPoint p1 = p + DPoint( 0, 1 );
            bool p_inside =  imageDomain.isInside( p ) && predicate( p );
            bool p1_inside = imageDomain.isInside( p1 ) && predicate( p1 );
            if ( prev_p_inside == prev_p1_inside )
              { // was not in a band.
                if ( p_inside != p1_inside ) 
                  { // start a band.
                    band_p = p; 
                    band_p1 = p1; 
                    sBand = 1;
                  }
                else // still not in band, but may be the start of a new band.
                  {
                    inside = p_inside;
                  }
              }
            else
              { // In a band.
                if ( ( p_inside != p1_inside ) && ( p_inside == prev_p_inside ) ) 
                  { // Still in the band.
                    sBand += 1;  // increase length
                  }
                else // end of band
                  {
                    bool ins1 = inside;
                    bool ins2 = ( p_inside == p1_inside )
                      ? p_inside
                      : ( k == 1 );
                    if ( sBand >= 3 )
                      {
                        DPoint prev_p =  p - DPoint( 1, 0 );
                        DPoint prev_p1 = p1 - DPoint( 1, 0 );
                        DPoint mid = ( band_p + prev_p ) / 2;
                        DPoint mid1 = ( band_p1 + prev_p1 ) / 2;
                        dac.addBand( toCGAL<Point>( band_p ), toCGAL<Point>( prev_p ), toCGAL<Point>( mid ),
                                     toCGAL<Point>( band_p1 ), toCGAL<Point>( prev_p1 ), toCGAL<Point>( mid1 ),
                                     ins1 != prev_p_inside, 
                                     ins2 != prev_p_inside );
                      }
                    inside = ins2;
                    if ( p_inside != p1_inside )
                      { // start a band.
                        band_p = p; 
                        band_p1 = p1; 
                        sBand = 1;
                      }
                    else
                      sBand = 0;
                  }
              }
            prev_p_inside = p_inside;
            prev_p1_inside = p1_inside;
          } // end of loop on x-axis.
      } // end of loop on y-axis.
  }
  {
    // Vertical bands
    ConstSubRange subRangeX = extDomain.subRange( 0, extDomain.lowerBound() );
    DPoint band_p, band_p1;
    for ( typename ConstSubRange::ConstIterator itx = subRangeX.begin(), 
            itxend = subRangeX.end(); itx != itxend; ++itx )
      {
        bool prev_p_inside  = false; // value at (x, y-1)
        bool prev_p1_inside = false; // value at (x+1, y-1)
        bool inside = false;
        ConstSubRange subRangeY = extDomain.subRange( 1, (*itx) + DPoint( 0, 1 ) );
        int sBand = 0;
        for ( typename ConstSubRange::ConstIterator ity = subRangeY.begin(), 
                ityend = subRangeY.end(); ity != ityend; ++ity )
          {
            DPoint p = *ity;
            DPoint p1 = p + DPoint( 1, 0 );
            bool p_inside =  imageDomain.isInside( p ) && predicate( p );
            bool p1_inside = imageDomain.isInside( p1 ) && predicate( p1 );
            if ( prev_p_inside == prev_p1_inside )
              { // was not in a band.
                if ( p_inside != p1_inside ) 
                  { // start a band.
                    band_p = p; 
                    band_p1 = p1; 
                    sBand = 1;
                  }
                else // still not in band, but may be the start of a new band.
                  {
                    inside = p_inside;
                  }
              }
            else
              { // In a band.
                if ( ( p_inside != p1_inside ) && ( p_inside == prev_p_inside ) ) 
                  { // Still in the band.
                    sBand += 1;  // increase length
                  }
                else // end of band
                  {
                    bool ins1 = inside;
                    bool ins2 = ( p_inside == p1_inside )
                      ? p_inside
                      : ( k == 1 );
                    if ( sBand >= 3 )
                      {
                        DPoint prev_p =  p - DPoint( 0, 1 );
                        DPoint prev_p1 = p1 - DPoint( 0, 1 );
                        DPoint mid = ( band_p + prev_p ) / 2;
                        DPoint mid1 = ( band_p1 + prev_p1 ) / 2;
                        dac.addBand( toCGAL<Point>( band_p ), toCGAL<Point>( prev_p ), toCGAL<Point>( mid ),
                                     toCGAL<Point>( band_p1 ), toCGAL<Point>( prev_p1 ), toCGAL<Point>( mid1 ),
                                     ins1 != prev_p_inside, 
                                     ins2 != prev_p_inside );
                      }
                    inside = ins2;
                    if ( p_inside != p1_inside )
                      { // start a band.
                        band_p = p; 
                        band_p1 = p1; 
                        sBand = 1;
                      }
                    else
                      sBand = 0;
                  }
              }
            prev_p_inside = p_inside;
            prev_p1_inside = p1_inside;
          } // end of loop on y-axis.
      } // end of loop on x-axis.
  }

}

double randomUniform()
{
  return (double) random() / (double) RAND_MAX;
}
///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef DAC<K> DigitalAffineComplex;
  typedef DigitalAffineComplex::Point Point;
  typedef DGtal::Z2i::Space Space;
  typedef DGtal::Z2i::Domain Domain;
  typedef DGtal::Z2i::Point DPoint;
  typedef DGtal::Z2i::KSpace KSpace;
  typedef DGtal::Z2i::Integer Integer;
  typedef DGtal::FreemanChain<Integer> FreemanChain;
  typedef KSpace::SCell SCell;

  using namespace DGtal;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("point-list,l", po::value<std::string>(), "Specifies the input shape as a list of 2d integer points, 'x y' per line.")
    ("freeman-chain,f", po::value<std::string>(), "Specifies the input shape as a closed contour, coded by a Freeman chaincode: x y 00121001...")
    ("image,i", po::value<std::string>(), "Specifies the input shape as a 2D image filename.")
    ("min,m", po::value<int>()->default_value(1), "Specifies the min threshold of the input 2D image.")
    ("max,M", po::value<int>()->default_value(255), "Specifies the max threshold of the input 2D image.")
    ("dac,d", po::value< std::vector<int> >(), "View the digital affine complex of the points with label [arg], or all of them when arg=-1; may be specified several times on the command line.")
    ("delaunay,D", po::value< std::vector<int> >(), "View the Delaunay triangulation of the points with label [arg], or all of them when arg=-1; may be specified several times on the command line.")
    ("vertices,v", po::value< std::vector<int> >(), "View the vertices (grid mode) of the triangulation with label [arg], or all of them when arg=-1; may be specified several times on the command line.")
    ("paving-vertices,V", po::value< std::vector<int> >(), "View the vertices (paving mode) of the triangulation with label [arg], or all of them when arg=-1; may be specified several times on the command line.")
    ("canonic,c", po::value< std::vector<int> >(), "View the constrained edges of the triangulation with label [arg], or all of them when arg=-1; may be specified several times on the command line.")
    ("border,b", po::value< std::vector<int> >(), "View the border complex of the triangulation with label [arg]; may be specified several times on the command line.")
    ("frontier,F", po::value< std::vector<int> >(), "View the frontier complex of the triangulation with label [arg]; may be specified several times on the command line.")
    ("surface,s", po::value< std::vector<int> >(), "View the boundary surface of the triangulation with label [arg]; may be specified several times on the command line.")
    ("mlp,L", po::value< std::vector<int> >(), "View the minimum length polygon (MLP) of the triangulation with label [arg]; may be specified several times on the command line.")
    ("blue-mlp,B", po::value< std::vector<int> >(), "View in blue the minimum length polygon (MLP) of the triangulation with label [arg]; may be specified several times on the command line.")
    ("random-inside,r",  po::value<double>()->default_value(1.0), "Specifies the % of inside points that are kept. Simulates a degradation noise." )
    ("random-outside,R",  po::value<double>()->default_value(1.0), "Specifies the % of outside points that are kept. Simulates a degradation noise." )
    ("remove-concavities,C", "Remove all concavities (useful for computing a DAC that is not the relative hull." )
    ("optimize,O",  "Flips the border so as to optimize its alignment with the mlp." )
    ("relative-hull,H", po::value< std::vector<int> >(), "Computes and visualize the relative hull of label [arg]." )
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
  if( ! parseOK || vm.count("help") || argc <=1 )
    {
      trace.info()<< "Generate a 2D triangulation from an arbitrary set of points. The triangulation provides a kind of piecewise linear approximation of the digitized set. The digitized set has label 0, the exterior points have label 1." <<std::endl << "Basic usage: " << std::endl
		  << "\2d-triangulation [options] -l <point-list> -v -1"<<std::endl
		  << general_opt << "\n";
      return 0;
    }


  std::set<DPoint> p,q;
  std::set<DPoint> p2,q2;
  DigitalAffineComplex dac;
  KSpace ks;
  double r_in = vm[ "random-inside" ].as<double>();
  double r_out = vm[ "random-outside" ].as<double>();

  trace.beginBlock("Construction of the shape");
  if ( vm.count( "point-list" ) )
    {
      DPoint lo, hi;
      std::vector<DPoint> pts;
      if ( readPointList<DPoint>( pts, lo, hi, vm["point-list"].as<std::string>() ) != 0 )
        {
          trace.error() << "Error reading file <" << vm["point-list"].as<std::string>() << ">." << std::endl;
          return 1;
        }
      // ks.init( lo, hi, true );
      for ( std::vector<DPoint>::const_iterator it = pts.begin(), ite = pts.end();
            it != ite; ++it )
        p.insert( *it );

      trace.beginBlock( "Computing border." );
      computeBorder<Space>( q, p );
      removeInside<Space>( p2, p );
      removeInside<Space>( q2, q );
      trace.endBlock();
    }
  if ( vm.count( "image" ) )
    {
      typedef ImageSelector < Z2i::Domain, unsigned char>::Type Image;
      typedef IntervalThresholder<Image::Value> Binarizer; 
      typedef Image::Domain Domain;
      std::string imageFileName = vm[ "image" ].as<std::string>();
      Image image = GenericReader<Image>::import( imageFileName ); 
      Binarizer b( vm[ "min" ].as<int>(), vm[ "max" ].as<int>() ); 
      PointFunctorPredicate<Image,Binarizer> predicate(image, b); 
      insertBands( dac, image, predicate, 1 );
      for ( Domain::ConstIterator it = image.domain().begin(), ite = image.domain().end();
            it != ite; ++it )
        if ( predicate( *it ) )
          p.insert( *it );

      trace.beginBlock( "Computing border." );
      computeBorder<Space>( q, p );
      removeInside<Space>( p2, p );
      removeInside<Space>( q2, q );
      // computeBorderInDomain<Space,Domain>( q, p, image.domain() );
      // removeInsideInDomain<Space,Domain>( p2, p, image.domain() );
      // removeInsideInDomain<Space,Domain>( q2, q, image.domain() );
      trace.endBlock();
    }
  if ( vm.count( "freeman-chain" ) )
    {
      trace.beginBlock( "Computing points from freeman chains." );
      // read freeman chain(s).
      std::string fileName = vm["freeman-chain"].as<std::string>();
      std::vector< FreemanChain > vectFcs = PointListReader< DPoint >
        ::getFreemanChainsFromFile<Integer>( fileName ); 
      // compute bounding box
      DPoint lo, hi;
      vectFcs[ 0 ].computeBoundingBox( lo[ 0 ], lo[ 1 ], hi[ 0 ], hi[ 1 ] );
      for ( unsigned int i = 1; i < vectFcs.size(); ++i )
        {
          // bool isClosed = vectFcs.at(i).isClosed(); 
          DPoint nlo, nhi;
          vectFcs[ i ].computeBoundingBox( nlo[ 0 ], nlo[ 1 ], nhi[ 0 ], nhi[ 1 ] );
          lo = lo.inf( nlo );
          hi = hi.sup( nhi );
        }
      DGtal::trace.info() << "- " << vectFcs.size() << " contours." << std::endl;
      DGtal::trace.info() << "- bbox=" << lo << ", " << hi << std::endl;
      // init Khalimsky space.
      ks.init( lo - DPoint::diagonal(1), hi + DPoint::diagonal(1), true );
      // extracts inner and outer points.
      for ( unsigned int i = 0; i < vectFcs.size(); ++i )
        {
          std::vector<DPoint> vectPts; 
          FreemanChain::getContourPoints( vectFcs.at(i), vectPts ); 
          GridCurve<KSpace> gc;
          bool ok = gc.initFromPointsVector( vectPts );
          DGtal::trace.info() << "- conversion to GridCurve is " <<  (ok ? "ok" : "KO" ) 
                              << std::endl;
          typedef GridCurve<KSpace>::SCellsRange SCellsRange;
          SCellsRange range = gc.getSCellsRange();
          for ( SCellsRange::ConstIterator it = range.begin(), ite = range.end();
                it != ite; ++it )
            {
              SCell in = ks.sDirectIncident( *it, ks.sOrthDir( *it ) );
              SCell out = ks.sIndirectIncident( *it, ks.sOrthDir( *it ) );
              if ( randomUniform() <= r_in ) 
                p2.insert( toCGAL<DPoint>( ks.sCoords( in ) ) );
              if ( randomUniform() <= r_out ) 
                q2.insert( toCGAL<DPoint>( ks.sCoords( out ) ) );
            }
        }
      trace.endBlock();
 
    }
    

  // typedef Ellipse2D<Z2i::Space> Ellipse; 
  // int a = 5, b = 1;
  // Ellipse2D<Z2i::Space> ellipse(Z2i::Point(0,0), a, b, 0.3 );
  // double h = 0.06125; // 06125; 
  // int N = 16;
  // GaussDigitizer<Z2i::Space,Ellipse> dig;  
  // dig.attach( ellipse );
  // dig.init( ellipse.getLowerBound()+Z2i::Vector(-1,-1),
  //           ellipse.getUpperBound()+Z2i::Vector(1,1), h ); 
  // typedef Flower2D<Z2i::Space> Flower; 
  // int a = 19, b = 9;
  // Flower2D<Z2i::Space> flower(Z2i::Point(0,0), 15, 2, 5, 0);
  // double h = 0.5; 
  // int N = 5;
  // GaussDigitizer<Z2i::Space,Flower> dig;  
  // dig.attach( flower );
  // dig.init( flower.getLowerBound()+Z2i::Vector(-1,-1),
  //           flower.getUpperBound()+Z2i::Vector(1,1), h ); 
  // Z2i::KSpace ks;
  // ks.init( dig.getLowerBound(), dig.getUpperBound(), true );
  // SurfelAdjacency<2> sAdj( true );
  // Z2i::SCell bel = Surfaces<Z2i::KSpace>::findABel( ks, dig, 1000 );
  // std::vector<Z2i::Point> boundaryPoints;
  // Surfaces<Z2i::KSpace>
  //   ::track2DBoundaryPoints( boundaryPoints, ks, sAdj, dig, bel );
  // Z2i::Curve c;
  // c.initFromVector( boundaryPoints );  
  // typedef Z2i::Curve::PointsRange Range; 
  // Range r = c.getPointsRange(); 
  // for ( Range::ConstIterator it = r.begin(), itE = r.end(); it != itE; ++it )
  //   {
  //     Z2i::Point P( *it );
  //     p.insert( P );
  //     p.insert( P + Z2i::Point( 2*N+2, -2*N ) );
  //     p.insert( Z2i::Point( P[ 1 ], P[ 0 ] ) );
  //   }
  // for ( int x = 0; x < 10*N; ++x )
  //   {
  //     p.insert( Z2i::Point( x, (12*x)/13 - 4*N - 2 ) );
  //   }
  trace.endBlock();




  trace.beginBlock( "Creating constrained Delaunay triangulation." );
  std::vector<Point> pts;
  for ( std::set<DPoint>::const_iterator it = p2.begin(), ite = p2.end();
        it != ite; ++it )
    pts.push_back( toCGAL<Point>( *it ) );
  dac.add( pts.begin(), pts.end(), 1, 0 );
  pts.clear();
  for ( std::set<DPoint>::const_iterator it = q2.begin(), ite = q2.end();
        it != ite; ++it )
    pts.push_back( toCGAL<Point>( *it ) );
  dac.add( pts.begin(), pts.end(), 0, 1 );
  dac.setInfiniteLabel( 1 );
  trace.endBlock();

  Board2D board;
  ViewerDAC<DigitalAffineComplex> viewer( board, dac );

  trace.beginBlock( "Visualizing Delaunay complex." );
  if ( vm.count( "paving-vertices" ) ) {
    std::vector<int> labels = vm[ "paving-vertices" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewPavingVertices( *it );
  }
  if ( vm.count( "border" ) ) {
    std::vector<int> labels = vm[ "border" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewBorder( *it, Color::Black, Color( 200, 200, 200, 0 ), 1.0 );
  }
  if ( vm.count( "delaunay" ) ) {
    std::vector<int> labels = vm[ "delaunay" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewTriangulationEdges( *it );
  }
  if ( vm.count( "canonic" ) ) {
    std::vector<int> labels = vm[ "canonic" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewConstrainedEdges( *it );
  }
  if ( vm.count( "surface" ) ) {
    std::vector<int> labels = vm[ "surface" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewBoundarySurface( *it, 4.0 );
  }
  if ( vm.count( "vertices" ) ) {
    std::vector<int> labels = vm[ "vertices" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewVertices( *it );
  }
  board.saveSVG("delaunay.svg");
  board.saveEPS("delaunay.eps");
  board.clear();
  trace.endBlock();

  trace.beginBlock( "Computing Digital Affine Complex." );
  if ( vm.count( "remove-concavities" ) )
    {
      dac.removeAllConcavities( 0 );
      dac.removeAllConcavities( 1 );
    }
  if ( vm.count( "optimize" ) )
    {
      while ( dac.flipBorder( 0 ) ) ;
    }
  if ( vm.count( "relative-hull" ) )
    {
      std::vector<int> labels = vm[ "relative-hull" ].as< std::vector<int> >();
      for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); 
            it != ite; ++it )
	{
	  unsigned int i = 0;
	  do 
	    {
	      viewer.viewRelativeHull( 0 );
	      viewer.viewTriangulationEdges( -1 );
	      viewer.viewVertices( -1 );
	      std::ostringstream sname, snams;
	      sname << "tmp/dac-rhull-" << *it << "-" << i << ".eps";
	      snams << "tmp/dac-rhull-" << *it << "-" << i << ".svg";
	      board.saveEPS( sname.str().c_str() );
	      board.saveSVG( snams.str().c_str() );
	      board.clear();
	      ++i;
	    }
	  while ( dac.relativeHull2( *it ) );
	}
    }
  trace.endBlock();    

  trace.beginBlock( "Visualizing Digital Affine Complex." );
  if ( vm.count( "paving-vertices" ) ) {
    std::vector<int> labels = vm[ "paving-vertices" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewPavingVertices( *it );
  }
  if ( vm.count( "border" ) ) {
    std::vector<int> labels = vm[ "border" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewBorder( *it, Color::Black, Color( 200, 200, 200, 0 ), 1.0 );
  }
  if ( vm.count( "relative-hull" ) )
    {
      std::vector<int> labels = vm[ "relative-hull" ].as< std::vector<int> >();
      for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); 
            it != ite; ++it )
        viewer.viewRelativeHull( *it );
    }
  if ( vm.count( "frontier" ) ) {
    std::vector<int> labels = vm[ "frontier" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewFrontier( *it, Color( 130, 130, 250 ), Color( 200, 200, 250, 0 ), 3.0 );
  }
  if ( vm.count( "dac" ) ) {
    std::vector<int> labels = vm[ "dac" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewTriangulationEdges( *it );
  }
  if ( vm.count( "canonic" ) ) {
    std::vector<int> labels = vm[ "canonic" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewConstrainedEdges( *it );
  }
  if ( vm.count( "surface" ) ) {
    std::vector<int> labels = vm[ "surface" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewBoundarySurface( *it, 4.0 );
  }
  if ( vm.count( "mlp" ) ) {
    std::vector<int> labels = vm[ "mlp" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewMLPContour( *it, Color( 0, 0, 255 ), 2.0, true );
  }
  if ( vm.count( "blue-mlp" ) ) {
    std::vector<int> labels = vm[ "blue-mlp" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewMLPContour( *it, Color::Blue, 3.0, false );
  }
  if ( vm.count( "vertices" ) ) {
    std::vector<int> labels = vm[ "vertices" ].as< std::vector<int> >();
    for ( std::vector<int>::const_iterator it = labels.begin(), ite = labels.end(); it != ite; ++it )
      viewer.viewVertices( *it );
  }

  // DPoint a( 154, 105 );
  // std::string specificStyle = a.className() + "/Paving";
  // board << DGtal::SetMode( a.className(), "Paving" );
  // board << DGtal::CustomStyle( specificStyle, 
  //                              new DGtal::CustomColors( Color::Black, Color::Magenta ) )
  // 	<< a
  // 	<< DPoint( 153, 107 )
  // 	<< DPoint( 156, 100 )
  // 	<< DPoint( 157, 98 );
  board.saveSVG("dac.svg");
  board.saveEPS("dac.eps");
  board.clear();
  trace.endBlock();
  return 0;
}
