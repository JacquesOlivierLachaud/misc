#include <climits>
#include <iostream>
#include <vector>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/Shapes.h>
#include <DGtal/shapes/ShapeFactory.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/geometry/curves/GridCurve.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/geometry/curves/ArithmeticalDSS.h>
#include <DGtal/geometry/curves/SaturatedSegmentation.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef Delaunay::Vertex_circulator       Vertex_circulator;
typedef Delaunay::Vertex_handle           Vertex_handle;
typedef Delaunay::Edge_iterator           Edge_iterator;
typedef Delaunay::Edge                    Edge;
typedef Delaunay::Finite_faces_iterator   Faces_iterator;
typedef Delaunay::Point                   Point;
typedef K::Vector_2                       Vector;
typedef Delaunay::Face_handle             Face_handle;


DGtal::Z2i::Point toDGtal(const Point &p)
{
  return DGtal::Z2i::Point (p.x(),p.y());
}
DGtal::Z2i::Vector toDGtal(const Vector &p)
{
  return DGtal::Z2i::Vector (p.x(),p.y());
}

DGtal::Z2i::Point::Coordinate
operator^( const DGtal::Z2i::Point & p1, 
	   const DGtal::Z2i::Point & p2 )
{
  return p1[ 0 ] * p2[ 1 ] - p1[ 1 ] * p2[ 0 ];
}

using namespace DGtal;

bool
emptyLatticeTriangle( const Delaunay & t,
                      const Vertex_handle & v1,
                      const Vertex_handle & v2,
                      const Vertex_handle & v3 )
{
  if ( t.is_infinite( v1 ) 
       || t.is_infinite( v2 ) 
       || t.is_infinite( v3 ) ) return false;
  Z2i::Point a( toDGtal( v1->point())),
    b(toDGtal( v2->point())),
    c(toDGtal( v3->point()));
  
  Z2i::Vector ab( b - a ), ac( c - a );
  int d = ab[ 0 ] * ac[ 1 ] - ab[ 1 ] * ac[ 0 ];
  return ( d == 1 ) || (d == -1 );
}
bool
emptyLatticeTriangle( const Delaunay & t, const Face_handle & f )
{
  if ( t.is_infinite( f ) ) return false;
  Z2i::Point a( toDGtal(f->vertex(0)->point())),
    b(toDGtal(f->vertex(1)->point())),
    c(toDGtal(f->vertex(2)->point()));
  
  Z2i::Vector ab( b - a ), ac( c - a );
  int d = ab[ 0 ] * ac[ 1 ] - ab[ 1 ] * ac[ 0 ];
  return ( d == 1 ) || (d == -1 );
}

int
twiceNbLatticePointsInTriangle( const Delaunay & t,
				const Vertex_handle & v1,
				const Vertex_handle & v2,
				const Vertex_handle & v3 )
{
  if ( t.is_infinite( v1 ) 
       || t.is_infinite( v2 ) 
       || t.is_infinite( v3 ) ) return 10000000;
  Z2i::Point a( toDGtal( v1->point())),
    b(toDGtal( v2->point())),
    c(toDGtal( v3->point()));
  
  Z2i::Vector ab( b - a ), ac( c - a );
  int d = ab[ 0 ] * ac[ 1 ] - ab[ 1 ] * ac[ 0 ];
  d = (d >= 0) ? (d - 1) : (-d - 1);
  return d;
}

int
twiceNbLatticePointsInTriangle( const Delaunay & t, const Face_handle & f )
{
  return twiceNbLatticePointsInTriangle( t, f->vertex(0), f->vertex(1), f->vertex(2) );
}

bool isEdgeElementary( const Delaunay & t,
		       const Vertex_handle & v1, 
		       const Vertex_handle & v2 )
{
  if ( t.is_infinite( v1 ) || t.is_infinite( v2 ) ) return false;
  Z2i::Point a( toDGtal( v1->point())),
    b(toDGtal( v2->point()));
  return (b-a).norm( Z2i::Point::L_1 ) == 1;
}

bool isQuadrilateral( const Delaunay & t,
		      const Vertex_handle & v1, 
		      const Vertex_handle & v2,
		      const Vertex_handle & v3,
		      const Vertex_handle & v4 )
{
  Z2i::Point a( toDGtal( v1->point())),
    b(toDGtal( v2->point())),
    c(toDGtal( v3->point())),
    d(toDGtal( v4->point()));
  return ( (( b-a )^( c-b )) > 0)
    && ( (( c-b )^( d-c )) > 0)
    && ( (( d-c )^( a-d )) > 0)
    && ( (( a-d )^( b-a )) > 0);
}

namespace DGtal {
  class DigitalCore {
  public:
    typedef Delaunay                      Triangulation;
    typedef Delaunay::Face_handle         FaceHandle;
    typedef Delaunay::All_faces_iterator  AllFacesIterator;
    typedef Delaunay::Finite_faces_iterator FiniteFacesIterator;
    typedef Delaunay::Edge                Edge;
    typedef Delaunay::Edge_iterator       EdgeIterator;
    typedef Delaunay::Edge_circulator     EdgeCirculator;
    typedef Delaunay::Finite_vertices_iterator FiniteVerticesIterator;
    typedef Delaunay::Vertex_circulator   VertexCirculator;
    typedef Delaunay::Vertex_handle       VertexHandle;
    typedef Delaunay::Point               CGALPoint;
    typedef K::Vector_2                   CGALVector;
    typedef std::set<Edge>                EdgeSet;
    typedef std::set<FaceHandle>          FaceSet;
    enum BoundaryAngle { FLAT, CONVEX, CONCAVE, EXTREMITY, NONE, MULTIPLE, ERROR };
    enum Flip { NO_NEED, FLIP_ERROR, FLIP_DONE };
    struct Concavity {
      VertexHandle v;
      CGALVector dir;
      int priority[ 2 ];
      
      inline Concavity( const DigitalCore & core, const VertexHandle & vertex, const CGALVector & direction )
        : v( vertex ), dir( direction )
      {
        const Triangulation & T = core.triangulation();
        Edge e1, e2;
        if ( ! getEdges( e1, e2, core ) )
	  {
	    trace.error() << "[DigitalCore::Concavity::Concavity] Concavity should be correct." << std::endl;
	    dir = CGALVector();
	  }
	else
	  {
	    Z2i::Vector v1 = toDGtal( T.segment( e1 ).to_vector() );
	    Z2i::Vector v2 = toDGtal( T.segment( e2 ).to_vector() );
	    priority[ 0 ] = v1.norm( Z2i::Vector::L_1 );
	    priority[ 1 ] = v2.norm( Z2i::Vector::L_1 );
	    if ( priority[ 1 ] < priority[ 0 ] )  std::swap( priority[ 0 ], priority[ 1 ] );
	  }
      }

      inline Concavity( const Concavity & other )
        : v( other.v ), dir( other.dir)
      {
        priority[ 0 ] = other.priority[ 0 ];
        priority[ 1 ] = other.priority[ 1 ];
      }

      inline Concavity & operator=( const Concavity & other )
      {
        if ( this != &other )
          {
            v = other.v;
            dir = other.dir;
            priority[ 0 ] = other.priority[ 0 ];
            priority[ 1 ] = other.priority[ 1 ];
          }
        return *this;
      }

      bool operator<( const Concavity & other ) const
      {
        return ( priority[ 0 ] > other.priority[ 0 ] )
          || ( ( priority[ 0 ] == other.priority[ 0 ] )
               && ( ( priority[ 1 ] > other.priority[ 1 ] ) 
                    || ( ( priority[ 1 ] == other.priority[ 1 ] ) 
                         && ( ( v < other.v )
                              || ( ( v == other.v )
                                   && ( ( dir[ 0 ] < other.dir[ 0 ] )
                                        || ( ( dir[ 0 ] == other.dir[ 0 ] ) 
                                             && ( dir[ 1 ] < other.dir[ 1 ] ) ) ) ) ) ) ) );
      }

      bool isFillable( const DigitalCore & core ) const
      {
        Edge e1, e2;
        if ( getEdges( e1, e2, core ) )
          {
            VertexHandle v1 = core.source( e1 );
            VertexHandle v2 = core.target( e2 );
            return core.twiceNbLatticePointsInTriangle( v1, v, v2 ) == 0;
          }
        return false;
      }

      bool getEdges( Edge & e1, Edge & e2, const DigitalCore & core ) const
      {
        const Triangulation & T = core.triangulation();
        std::vector<Edge> tEdges;
        std::vector<Edge> sEdges;
        core.getBoundaryPairs( tEdges, sEdges, v );
        int j = -1;
        // int best = 0;
        for ( int i = 0; i < tEdges.size(); ++i )
          {
            CGALVector v0 = T.segment( tEdges[ i ] ).to_vector();
            CGALVector v1 = T.segment( sEdges[ i ] ).to_vector();
            if ( ( ( toDGtal( v0 ) ^ toDGtal( v1 ) ) > 0 ) // only concave vertices.
		 && ( ( v1 - v0 ) == dir ) ) 
              {
		j = i;
                // int val = ( v1 - v0 ) * dir;
                // if ( val > best ) { j = i; best = val; }
              }
          }
        if ( j == -1 )
          {
            // trace.info() << "[DigitalCore::Concavity::getEdges] No edge pair is correct for direction." 
	    // << std::endl;
            return false;
          }
        e1 = tEdges[ j ];
        e2 = sEdges[ j ];
        return true;
      }

      Flip flipExternalEdge( DigitalCore & core )
      {
        const Triangulation & T = core.triangulation();
        Edge e1, e2;
        bool ok = getEdges( e1, e2, core );
        if ( ! ok ) return FLIP_ERROR;
        if ( core.isEdgeInterior( e1 ) || core.isEdgeInterior( e2 ) )
          {
            trace.error() << "[DigitalCore::Concavity::flipExternalEdge] At least one boundary edge is interior." << std::endl;
            return FLIP_ERROR;
          }
        if ( e1.first == e2.first ) // concavity is reduced to a triangle, no more flip needed.
          {
            // trace.info() << "[DigitalCore::Concavity::flipExternalEdge] No more flip are needed." << std::endl;
            return NO_NEED;
          }
        Edge fe = Edge( e1.first, T.ccw( e1.second ) );
	if ( core.isEdgeQuadrilateral( fe ) ) 
          {
            // trace.info() << "[DigitalCore::Concavity::flipExternalEdge] flip edge." << std::endl;
            core.secureFlip( fe );
            return FLIP_DONE;
          }
        trace.info() << "[DigitalCore::Concavity::flipExternalEdge] Edge is not quadrilateral." 
		     << core.source( e1 )->point() << " -> " << core.target( e1 )->point()
		     << std::endl;
        return FLIP_ERROR;
      }

      
    };
  private:
    /// The triangulation. It is copied since it is modified.
    Triangulation T;
    /// Current boundary as a set of edges.
    EdgeSet myBoundary;
    /// Current interior faces as a set of face handles.
    FaceSet myInterior;

  public:
    
    
    DigitalCore() {}

    DigitalCore( const Triangulation & t )
      : T( t )
    {
      computeBasicEdges();
    }
    
    void init( const Triangulation & t )
    {
      myBoundary.clear();
      myInterior.clear();
      T = t;
      computeBasicEdges();
    }

    
    const Triangulation & triangulation() const
    { return T; }
    
    Triangulation & triangulation()
    { return T; }
    
    
    const EdgeSet & boundary() const
    { return myBoundary; }

    
    const FaceSet & interior() const
    { return myInterior; }

    /// The incident face is considered ccw.
    
    VertexHandle source( const Edge & e ) const
    {
      return e.first->vertex( T.ccw( e.second ) );
    }

    /// The incident face is considered ccw.
    
    VertexHandle target( const Edge & e ) const
    {
      return e.first->vertex( T.cw( e.second ) );
    }

    
    Edge nextCCWAroundFace( const Edge & e ) const
    {
      return Edge( e.first, T.ccw( e.second ) );
    }

    
    Edge nextCWAroundFace( const Edge & e ) const
    {
      return Edge( e.first, T.cw( e.second ) );
    }

    /// @return the next edge around the source vertex, such that this
    /// edge has the same source and is CCW around it.
    
    Edge nextCCWAroundSourceVertex( const Edge & e ) const
    {
      Edge e1 = nextCWAroundFace( e );
      Edge n = T.mirror_edge( e1 );
      if ( source( n ) != source( e ) )
        trace.error() << "[DigitalCore::nextCCWAroundSourceVertex] sources are not consistent." << std::endl;
      if ( target( n ) == target( e ) )
        trace.error() << "[DigitalCore::nextCCWAroundSourceVertex] targets are equal." << std::endl;
      return n;
    }

    /// @return the next edge around the source vertex, such that this
    /// edge has the same source and is CW around it.
    
    Edge nextCWAroundSourceVertex( const Edge & e ) const
    {
      return nextCCWAroundFace( T.mirror_edge( e ) );
    }

  public:
    
    bool isEdgeBasic( const Edge & e ) const
    {
      VertexHandle v1 = e.first->vertex( T.ccw( e.second ) );
      VertexHandle v2 = e.first->vertex( T.cw( e.second ) );
      if ( T.is_infinite( v1 ) || T.is_infinite( v2 ) ) 
        return false;
      Z2i::Point a( toDGtal( v1->point())),
        b(toDGtal( v2->point()));
      return (b-a).norm( Z2i::Point::L_1 ) == 1;
    }

    inline
    bool isEdgeBoundary( const Edge & e ) const
    {
      return myBoundary.find( e ) != myBoundary.end();
    }

    inline
    bool isEdgeInterior( const Edge & e ) const
    {
      return myInterior.find( e.first ) != myInterior.end();
    }

    inline
    bool isFaceInterior( const FaceHandle & f ) const
    {
      return myInterior.find( f ) != myInterior.end();
    }

    inline
    int twiceNbLatticePointsInTriangle( const Face_handle & f ) const
    {
      return this->twiceNbLatticePointsInTriangle( f->vertex(0), f->vertex(1), f->vertex(2) );
    }

    int
    twiceNbLatticePointsInTriangle( const Vertex_handle & v1,
                                    const Vertex_handle & v2,
                                    const Vertex_handle & v3 ) const
    {
      if ( T.is_infinite( v1 ) 
           || T.is_infinite( v2 ) 
           || T.is_infinite( v3 ) ) return INT_MAX;
      Z2i::Point a( toDGtal( v1->point())),
        b(toDGtal( v2->point())),
        c(toDGtal( v3->point()));
      Z2i::Vector ab( b - a ), ac( c - a );
      int d = ab[ 0 ] * ac[ 1 ] - ab[ 1 ] * ac[ 0 ];
      d = (d >= 0) ? (d - 1) : (-d - 1);
      return d;
    }

    bool isEdgeQuadrilateral( const Edge & e1 ) const
    {
      Edge e2 = T.mirror_edge( e1 );
      VertexHandle v1 = e1.first->vertex( e1.second );
      VertexHandle v2 = e1.first->vertex( T.ccw( e1.second ) );
      VertexHandle v3 = e2.first->vertex( e2.second );
      VertexHandle v4 = e1.first->vertex( T.cw( e1.second ) );
      Z2i::Point a( toDGtal( v1->point())),
        b(toDGtal( v2->point())),
        c(toDGtal( v3->point())),
        d(toDGtal( v4->point()));
      return ( (( b-a )^( c-b )) > 0)
        && ( (( c-b )^( d-c )) > 0)
        && ( (( d-c )^( a-d )) > 0)
        && ( (( a-d )^( b-a )) > 0);
    }

    void secureFlip( const Edge & e1 )
    {
      std::vector<Edge> save_boundary;
      Edge e2 = T.mirror_edge( e1 );
      Edge e1_a( e1.first, T.cw( e1.second ) );
      Edge e1_b( e1.first, T.ccw( e1.second ) );
      Edge e2_a( e2.first, T.cw( e2.second ) );
      Edge e2_b( e2.first, T.ccw( e2.second ) );
      if ( isEdgeBoundary( e1_a ) ) save_boundary.push_back( e1_a );
      if ( isEdgeBoundary( e1_b ) ) save_boundary.push_back( e1_b );
      if ( isEdgeBoundary( e2_a ) ) save_boundary.push_back( e2_a );
      if ( isEdgeBoundary( e2_b ) ) save_boundary.push_back( e2_b );
      // get mirror edges (since edge flip will change incident faces).
      for ( std::vector<Edge>::iterator it = save_boundary.begin(), itend = save_boundary.end();
            it != itend; ++it )
        {
          Edge e = *it;
          myBoundary.erase( e );
          *it = T.mirror_edge( e );
        }
      // flip edge.
      T.flip( e1.first, e1.second );
      // push back in boundary mirror of mirror edges.
      for ( std::vector<Edge>::iterator it = save_boundary.begin(), itend = save_boundary.end();
            it != itend; ++it )
        myBoundary.insert( T.mirror_edge( *it ) );
    }

    void getBoundaryPairs( std::vector<Edge> & tEdges,
                           std::vector<Edge> & sEdges,
                           const VertexHandle & v ) const
    {
      if ( T.is_infinite( v ) ) return;
      EdgeCirculator ci = T.incident_edges( v );
      if ( ci != 0 ) {
	Edge e = *ci;
	if ( source( e ) != v ) e = T.mirror_edge( e );
	if ( source( e ) != v )
	  trace.error() << "[DigitalCore::getBoundaryPairs] Invalid incident edge." << std::endl;
	Edge start = e;
	bool found = false;
	do {
	  if ( ( ! T.is_infinite( e ) ) && isEdgeBoundary( e ) ) 
	    {
	      found = true;
	      break;
	    }
	  else
	    {
	      e = nextCCWAroundSourceVertex( e );
	    }
	} while ( e != start );
	if ( found ) 
	  {
	    start = e;
	    do {
	      if ( ( ! T.is_infinite( e ) ) && isEdgeBoundary( e ) ) 
		sEdges.push_back( e );
	      Edge ep = nextCWAroundFace( e );
	      if ( ( ! T.is_infinite( ep ) ) && isEdgeBoundary( ep ) )
		tEdges.push_back( ep );
	      e = nextCCWAroundSourceVertex( e );
	    } while ( e != start );
	  }
      }
      unsigned int s = sEdges.size() + tEdges.size();
      if ( sEdges.size() != tEdges.size() )
        {
          trace.error() << "[DigitalCore::getBoundaryPairs] Odd number of boundary edges." << std::endl;
          trace.error() << "#sEdges=" << sEdges.size();
	  for ( int i = 0; i < sEdges.size(); ++i )
	    trace.error() << " sEdges[" << i << "]=" << source( sEdges[ i ] )->point()  
			  << " -> " << target( sEdges[ i ] )->point();
	  trace.error() << std::endl;
          trace.error() << "#tEdges=" << tEdges.size();
	  for ( int i = 0; i < tEdges.size(); ++i )
	    trace.error() << " tEdges[" << i << "]=" << source( tEdges[ i ] )->point()  
			  << " -> " << target( tEdges[ i ] )->point();
	  trace.error() << std::endl;
	  sEdges.clear();
	  tEdges.clear();
          return;
        }
    }

    void makeConcavities( std::vector<Concavity> & concavities, 
                          const VertexHandle & v )
    {
      std::vector<Edge> sEdges;
      std::vector<Edge> tEdges;
      getBoundaryPairs( tEdges, sEdges, v ); 
      unsigned int s = tEdges.size() + sEdges.size();
      if ( s == 0 ) return;
      s >>= 1;
      for ( unsigned int i = 0; i < s; ++i )
        {
          CGALVector v0 = T.segment( tEdges[ i ] ).to_vector();
          CGALVector v1 = T.segment( sEdges[ i ] ).to_vector();
          if ( ( ( toDGtal( v0 ) ^ toDGtal( v1 ) ) > 0 )
	       && ( ! isEdgeInterior( tEdges[ i ] ) )
	       && ( ! isEdgeInterior( sEdges[ i ] ) ) )
            {
              concavities.push_back( Concavity( *this, v, v1 - v0 ) );
            }
        }
    }

    bool fillConcavities()
    {
      // Prepare queue
      std::priority_queue<Concavity> Q;
      std::vector<Concavity> concavities;
      for ( FiniteVerticesIterator it = T.finite_vertices_begin(), itend = T.finite_vertices_end();
            it != itend; ++it )
        {
          makeConcavities( concavities, it );
          for ( std::vector<Concavity>::const_iterator it2 = concavities.begin(), it2end = concavities.end();
                it2 != it2end; ++it2 )
            Q.push( *it2 );
          concavities.clear();
        }
      bool fill = false;
      while ( ! Q.empty() )
        {
          Concavity cc = Q.top();
          Q.pop();
          if ( cc.isFillable( *this ) )
            {
              // trace.info() << "Concavity " << cc.v->point() << " is fillable." << std::endl;
              Flip res;
              while ( ( res = cc.flipExternalEdge( *this ) ) == FLIP_DONE )
                ;
              // Fill face
              if ( res == NO_NEED )
                {
                  Edge e1, e2, ne;
                  cc.getEdges( e1, e2, *this );
                  ne = Edge( e2.first, T.ccw( e2.second ) );
                  // trace.info() << "Filling with " << source( ne )->point() << " -> " << target( ne )->point() << std::endl;
                  FaceHandle f = e1.first;
                  myInterior.insert( f );
                  Edge nedge = T.mirror_edge( ne );
                  myBoundary.erase( e1 );
                  myBoundary.erase( e2 );
		  if ( isEdgeBoundary( nedge ) )
		    trace.error() << "[DigitalCore::fillConcavities] Already in boundary." << std::endl;
		  Edge mirror_nedge = T.mirror_edge( nedge );
		  if ( isEdgeBoundary( mirror_nedge ) )
		    { // Boundaries merge.
		      myBoundary.erase( mirror_nedge );
		    }
                  else
		    { // Boundary moves.
		      myBoundary.insert( nedge );
		      makeConcavities( concavities, source( nedge ) );
		      makeConcavities( concavities, target( nedge ) );
		      // trace.info() << "Adding " << concavities.size() << " concavities." << std::endl;
		      for ( std::vector<Concavity>::const_iterator it2 = concavities.begin(), it2end = concavities.end();
			    it2 != it2end; ++it2 )
			Q.push( *it2 );
		      concavities.clear();
		    }
                  fill = true;
                }
            }
        }
      return fill;
    }

    int countExteriorEdges( VertexHandle v ) const
    {
      EdgeCirculator ci = T.incident_edges( v );
      EdgeCirculator ce = ci;
      int n = 0;
      if ( ci != 0 )
        do {
          if ( ( ! isEdgeBoundary( *ci ) ) 
               && ( ! isEdgeBoundary( T.mirror_edge( *ci ) ) ) 
               && ( ! isEdgeInterior( *ci ) )
               && ( ! isEdgeInterior( T.mirror_edge( *ci ) ) ) )
            ++n;
        } while ( ++ci != ce );
      return n;
    }

    BoundaryAngle boundaryAngle( VertexHandle v ) const
    {
      std::vector<Edge> bEdges;
      EdgeCirculator ci = T.incident_edges( v );
      EdgeCirculator ce = ci;
      int n = 0;
      if ( ci != 0 )
        do {
          bool e_b = isEdgeBoundary( *ci );
          bool em_b = isEdgeBoundary( T.mirror_edge( *ci ) );
          if ( e_b )    bEdges.push_back( *ci );
          if ( em_b )   bEdges.push_back( T.mirror_edge( *ci ) );
        } while ( ++ci != ce );
      if ( bEdges.size() == 0 ) return NONE;
      else if ( bEdges.size() == 1 ) return EXTREMITY;
      else if ( bEdges.size() >= 3 ) return MULTIPLE;
      if ( source( bEdges[ 0 ] ) == v ) 
        std::swap( bEdges[ 0 ], bEdges[ 1 ] );
      if ( target( bEdges[ 0 ] ) != v || source( bEdges[ 1 ] ) != v )
        {
          trace.info() << "V  = " << v->point() << std::endl;
          trace.info() << "E0 = " << source( bEdges[ 0 ] )->point()
                       << " -> " << target( bEdges[ 0 ] )->point() << std::endl;
          trace.info() << "E1 = " << source( bEdges[ 1 ] )->point()
                       << " -> " << target( bEdges[ 1 ] )->point() << std::endl;
          trace.error() << "[DigitalCore::boundaryAngle] Invalid boundary. Probably infinite." << std::endl;
          return ERROR;
        }
      if ( isFaceInterior( bEdges[ 0 ].first ) || isFaceInterior( bEdges[ 1 ].first ) )
        {
          trace.info() << "[DigitalCore::boundaryAngle] Folding of boundary." << std::endl;
          return ERROR;
        }
      Z2i::Point a( toDGtal( source( bEdges[ 0 ] )->point() ) );
      Z2i::Point b( toDGtal( v->point() ) );
      Z2i::Point c( toDGtal( target( bEdges[ 1 ] )->point() ) );
      int d = (b-a)^(c-b);
      if ( d > 0 ) return CONCAVE;
      else if ( d < 0 ) return CONVEX;
      else return FLAT;
    }

    int extend()
    {
      trace.beginBlock("[DigitalCore] Extend boundary.");

      // Queue for computing the digital core.
      FaceSet priorityQ;
      
      // Prepare queue.
      for ( EdgeSet::const_iterator it = myBoundary.begin(), itend = myBoundary.end();
            it != itend; ++it )
        priorityQ.insert( it->first );

      // Start extension
      bool isMarked[ 3 ];
      int nb = 0;
      while ( ! priorityQ.empty() )
        {
          // trace.progressBar( 1.0, 1.0+log( (double) priorityQ.size() ) );
          FaceSet::iterator it = priorityQ.begin();
          FaceHandle face_h = *it;
          priorityQ.erase( it );
          // infinite cells are not processed.
          if ( T.is_infinite( face_h ) ) continue;
          // cells containing integer points in the interior are not processed.
          if ( twiceNbLatticePointsInTriangle( face_h ) != 0 ) continue;
          // checking cell faces for extension.
          if ( checkFaceExtension( isMarked, face_h ) )
            {
              extendFace( priorityQ, isMarked, face_h );
              ++nb;
            }
        }
      trace.info() << "- cells flipped = " << nb << std::endl;
      trace.info() << "- boundary has " << myBoundary.size() << " edges." << std::endl;
      trace.endBlock();
      return nb;
    }
    
    int flipUpdate()
    {
      int nb_flip = 0;
      EdgeIterator itnext;
      for( EdgeIterator it = T.edges_begin(), itend=T.edges_end();
           it != itend; it = itnext )
      {
	// vertex(cw(i)) and vertex(ccw(i)) of f.
	itnext = it; ++itnext;
        //  A ---- D
        //  | \    |
        //  |  \e1 |
        //  | e2\  |
        //  |    \ |
        //  C ---- B
	Edge e1 = *it;
	if ( isEdgeBoundary( e1 ) ) continue;
	if ( isEdgeInterior( e1 ) ) continue;
	Edge e2 = T.mirror_edge( e1 );
	if ( ! isEdgeQuadrilateral( e1 ) ) continue;
        int nb_f1 = twiceNbLatticePointsInTriangle( e1.first ); // ABD
	int nb_f2 = twiceNbLatticePointsInTriangle( e2.first ); // ACB
        VertexHandle A = e1.first->vertex( T.ccw( e1.second ) );
        VertexHandle B = e1.first->vertex( T.cw( e1.second ) );
        VertexHandle C = e2.first->vertex( e2.second );
        VertexHandle D = e1.first->vertex( e1.second );
	int nb_flip_f1 = twiceNbLatticePointsInTriangle( D, A, C );
	int nb_flip_f2 = twiceNbLatticePointsInTriangle( D, C, B );
	int nb_min = nb_f1 <= nb_f2 ? nb_f1 : nb_f2;
	int nb_flip_min = nb_flip_f1 <= nb_flip_f2 ? nb_flip_f1 : nb_flip_f2;
        int nbOutEdge_A = countExteriorEdges( A );
        int nbOutEdge_B = countExteriorEdges( B );
        int nbOutEdge_C = countExteriorEdges( C );
        int nbOutEdge_D = countExteriorEdges( D );
        bool concave_A = boundaryAngle( A ) == CONCAVE;
        bool concave_B = boundaryAngle( B ) == CONCAVE;
        bool concave_C = boundaryAngle( C ) == CONCAVE;
        bool concave_D = boundaryAngle( D ) == CONCAVE;
        // trace.info() << "A(" << concave_A << ":" << nbOutEdge_A << ") "
        //              << "B(" << concave_B << ":" << nbOutEdge_B << ") "
        //              << "C(" << concave_C << ":" << nbOutEdge_C << ") "
        //              << "D(" << concave_D << ":" << nbOutEdge_D << ") "
        //              << std::endl;
        int nbOutEdgeInConcavity = ( concave_A ? nbOutEdge_A : INT_MAX );
        nbOutEdgeInConcavity = std::min( nbOutEdgeInConcavity, concave_B ? nbOutEdge_B : INT_MAX );
        nbOutEdgeInConcavity = std::min( nbOutEdgeInConcavity, concave_C ? nbOutEdge_C : INT_MAX );
        nbOutEdgeInConcavity = std::min( nbOutEdgeInConcavity, concave_D ? nbOutEdge_D : INT_MAX );
        int nbOutEdgeInConcavity_flip = ( concave_A ? nbOutEdge_A-1 : INT_MAX );
        nbOutEdgeInConcavity_flip = std::min( nbOutEdgeInConcavity_flip, concave_B ? nbOutEdge_B-1 : INT_MAX );
        nbOutEdgeInConcavity_flip = std::min( nbOutEdgeInConcavity_flip, concave_C ? nbOutEdge_C+1 : INT_MAX );
        nbOutEdgeInConcavity_flip = std::min( nbOutEdgeInConcavity_flip, concave_D ? nbOutEdge_D+1 : INT_MAX );
        // int nbOutEdgeInConcavity = ( concave_A ? nbOutEdge_A : 0 )
        //   + ( concave_B ? nbOutEdge_B : 0 )
        //   + ( concave_C ? nbOutEdge_C : 0 )
        //   + ( concave_D ? nbOutEdge_D : 0 );
        // int nbOutEdgeInConcavity_flip = ( concave_A ? nbOutEdge_A-1 : 0 )
        //   + ( concave_B ? nbOutEdge_B-1 : 0 )
        //   + ( concave_C ? nbOutEdge_C+1 : 0 )
        //   + ( concave_D ? nbOutEdge_D+1 : 0 );

	// if ( ( nb_flip_min < nb_min )
        //      || ( ( nb_flip_min == nb_min )
        //           && ( nbOutEdgeInConcavity_flip < nbOutEdgeInConcavity ) ) )
        if ( nbOutEdgeInConcavity_flip < nbOutEdgeInConcavity ) // sufficient ?
          {
	    // std::cout << "flipped " << e1.first->vertex( e1.second )->point()
	    //           << "->" << e2.first->vertex( e2.second )->point()
	    //           << std::endl;
            std::vector<Edge> save_boundary;
            Edge e1_a( e1.first, T.cw( e1.second ) );
            Edge e1_b( e1.first, T.ccw( e1.second ) );
            Edge e2_a( e2.first, T.cw( e2.second ) );
            Edge e2_b( e2.first, T.ccw( e2.second ) );
            if ( isEdgeBoundary( e1_a ) ) save_boundary.push_back( e1_a );
            if ( isEdgeBoundary( e1_b ) ) save_boundary.push_back( e1_b );
            if ( isEdgeBoundary( e2_a ) ) save_boundary.push_back( e2_a );
            if ( isEdgeBoundary( e2_b ) ) save_boundary.push_back( e2_b );
            // get mirror edges (since edge flip will change incident faces).
            for ( std::vector<Edge>::iterator it = save_boundary.begin(), itend = save_boundary.end();
                  it != itend; ++it )
              {
                Edge e = *it;
                myBoundary.erase( e );
                *it = T.mirror_edge( e );
              }
            // flip edge.
	    T.flip( e1.first, e1.second );
            // push back in boundary mirror of mirror edges.
            for ( std::vector<Edge>::iterator it = save_boundary.begin(), itend = save_boundary.end();
                  it != itend; ++it )
              myBoundary.insert( T.mirror_edge( *it ) );
            ++nb_flip;
	  }
      }
      return nb_flip;
    }

  protected:
    void computeBasicEdges()
    {
      trace.beginBlock("[DigitalCore] Compute basic edges.");
      DGtal::uint64_t nb = 0;
      for( EdgeIterator it = T.edges_begin(), itend = T.edges_end();
           it != itend; ++it )
        {
          if ( isEdgeBasic( *it ) ) // || checkPlanarVFacet( f ) )
            {
	      myBoundary.insert( *it );
	      myBoundary.insert( T.mirror_edge( *it ) );
              ++nb;
            }
        }
      trace.info() << "- nb basic edges = " << nb << std::endl;
      trace.endBlock();
    }

    bool checkFaceExtension( bool isBoundary[ 3 ], const FaceHandle & face ) const
    {
      unsigned int n = 0;
      for ( int i = 0; i < 3; ++i )
        {
          isBoundary[ i ] = isEdgeBoundary( Edge( face, i ) );
          if ( isBoundary[ i ] ) ++n;
        }
      // At least 2 are marked, by convexity, we can close the gap
      // to the further faces of the cell.
      bool propagate = n >= 2;
      return propagate;
    }

    void extendFace( FaceSet & priorityQ, bool isBoundary[ 3 ], const FaceHandle & face )
    { // Switch face to interior
      myInterior.insert( face );
      for ( unsigned int i = 0; i < 3; ++i ) {
        if ( ! isBoundary[ i ] )
          {
            Edge nedge = T.mirror_edge( Edge( face, i ) );
            priorityQ.insert( nedge.first );
            myBoundary.insert( nedge );
          }
        else
          {
            myBoundary.erase( Edge( face, i ) );
          }
      }
    }

  };
}


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main ( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("point-list,l", po::value<std::string>(), "Specifies the input shape as a list of 2d integer points, 'x y' per line.")
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
		  << "\2d-triangulation [options] -l <point-list>"<<std::endl
		  << general_opt << "\n";
      return 0;
    }



  std::vector<Z2i::Point> error_list;
  // error_list.push_back( Z2i::Point( 30, -15 ) );
  // error_list.push_back( Z2i::Point( 32, -16 ) );
  // error_list.push_back( Z2i::Point( 8, 26 ) );
  // error_list.push_back( Z2i::Point( 16, 25 ) );

  Delaunay t;
  
  trace.beginBlock("Construction the shape");
  typedef Ellipse2D<Z2i::Space> Ellipse; 
  int a = 5, b = 3;
  Ellipse2D<Z2i::Space> ellipse(Z2i::Point(0,0), a, b, 0.3 );
  double h = 0.25; //06125; 
  GaussDigitizer<Z2i::Space,Ellipse> dig;  
  dig.attach( ellipse );
  dig.init( ellipse.getLowerBound()+Z2i::Vector(-1,-1),
            ellipse.getUpperBound()+Z2i::Vector(1,1), h ); 
  // typedef Flower2D<Z2i::Space> Flower; 
  // int a = 19, b = 9;
  // Flower2D<Z2i::Space> flower(Z2i::Point(0,0), 15, 2, 5, 0);
  // double h = 0.5; 
  // GaussDigitizer<Z2i::Space,Flower> dig;  
  // dig.attach( flower );
  // dig.init( flower.getLowerBound()+Z2i::Vector(-1,-1),
  //           flower.getUpperBound()+Z2i::Vector(1,1), h ); 
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
  trace.endBlock();


  trace.beginBlock("Delaunay");
  for(Range::ConstIterator it=r.begin(), itend=r.end(); it != itend;
      ++it)
    { 
      t.insert( Point( (*it)[0], (*it)[1]));
      t.insert( Point( (*it)[0] + 3 + (int) ceil( ((double)b)/h ),
		       (*it)[1] - 3 - (int) ceil( ((double)a)/h ) ));
      t.insert( Point( (*it)[0] + 3 + (int) ceil( ((double)a)/h ),
		       (*it)[1] + 3 + (int) ceil( ((double)b)/h ) ));
    }
  trace.endBlock();

  std::cout << "number of vertices :  " ;
  std::cout << t.number_of_vertices() << std::endl;
  std::cout << "number of faces :  " ;
  std::cout << t.number_of_faces() << std::endl;

  trace.beginBlock( "Computing digital core." );
  DigitalCore core( t );

  bool fill;
  do {
    trace.info() << "------------ fill concavities ---------------------------------" << std::endl;
    fill = core.fillConcavities();
  } while ( fill );

  // int nb_flip;
  // do {
  //   trace.info() << "------------------------------------------------------------" << std::endl;
  //   core.extend();
  //   nb_flip = core.flipUpdate();
  //   trace.info() << "- nb flip = " << nb_flip << std::endl;
  // } while ( nb_flip != 0 );

  trace.endBlock();

  {
    DGtal::Board2D board;
    
    Z2i::Point dP;
    board << CustomStyle( dP.className(), 
                          new CustomPen( Color(0,0,0), Color(230,230,230), 1, 
                                         Board2D::Shape::SolidStyle,
                                         Board2D::Shape::RoundCap,
                                         Board2D::Shape::RoundJoin ));
    for( Range::ConstIterator it=r.begin(), itend=r.end(); it != itend;
         ++it)
      board << *it;
    board << CustomStyle( dP.className(), 
                          new CustomPen( Color(0,0,0), Color(255,0,255), 1, 
                                         Board2D::Shape::SolidStyle,
                                         Board2D::Shape::RoundCap,
                                         Board2D::Shape::RoundJoin ));
    for ( std::vector<Z2i::Point>::const_iterator it = error_list.begin(),
	    itend = error_list.end(); it != itend; ++it )
      board << *it;
    for ( DigitalCore::EdgeSet::const_iterator it = core.boundary().begin(),
            itend = core.boundary().end(); it != itend; ++it )
      {
        Z2i::Point a( toDGtal( core.source( *it )->point() ) );
        Z2i::Point b( toDGtal( core.target( *it )->point() ) );
        board.setPenColor(DGtal::Color::Blue);
        board.setFillColor( DGtal::Color::None );
        board.setLineWidth( 3.0 );
        board.drawLine(a[0],a[1],b[0],b[1]);
      }

    for( Faces_iterator it = core.triangulation().finite_faces_begin(), 
           itend=core.triangulation().finite_faces_end();
         it != itend; ++it)
    {
      Z2i::Point a( toDGtal(it->vertex(0)->point())),
	b(toDGtal(it->vertex(1)->point())),
	c(toDGtal(it->vertex(2)->point()));

      if ( emptyLatticeTriangle( core.triangulation(), it ) ) //( ( d == 1 ) || (d == -1 ) )
        {
          board.setPenColor(DGtal::Color::Green);
          board.setFillColor( DGtal::Color::None );
          board.setLineWidth( 1.0 );
          board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
        }
      else
        {
          board.setPenColor(DGtal::Color::Red);
          board.setFillColor( DGtal::Color::None );
          board.setLineWidth( 1.0 );
          board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
        }
    }
    board.saveSVG("core.svg");
    board.saveEPS("core.eps");

  }


  // GridCurve
  Z2i::Curve gc;
  gc.initFromPointsRange( r.begin(), r.end() );
  typedef Z2i::Curve::PointsRange::ConstIterator ConstIterator;
  typedef ArithmeticalDSS<ConstIterator,int,4> DSS4;
  typedef SaturatedSegmentation<DSS4> Segmentation;
  //Segmentation
  Z2i::Curve::PointsRange range = gc.getPointsRange();
  DSS4 dss4RecognitionAlgorithm;
  Segmentation theSegmentation( range.begin(), range.end(), dss4RecognitionAlgorithm );
         

  DGtal::Board2D board;

  Z2i::Point dP;
  board << CustomStyle( dP.className(), 
                        new CustomPen( Color(0,0,0), Color(230,230,230), 1, 
                                       Board2D::Shape::SolidStyle,
                                       Board2D::Shape::RoundCap,
                                       Board2D::Shape::RoundJoin ));
  for(Range::ConstIterator it=r.begin(), itend=r.end(); it != itend;
      ++it)
    board << *it;

  // // Implements flips on Delaunay on the sole criteria to minimize edge lengths !
  // // Interesting.
  // bool changed = true;
  // while ( changed ) 
  //   {
  //     typedef Edge_iterator EdgeIterator;
  //     typedef Vertex_handle VertexHandle;
  //     changed = false;
  //     EdgeIterator itnext;
  //     for( EdgeIterator it = t.edges_begin(), itend = t.edges_end();
  //          it != itend; it = itnext )
  // 	{
  // 	  itnext = it; ++itnext;
  // 	  Edge e1 = *it;
  // 	  Edge e2 = t.mirror_edge( e1 );
  // 	  if ( t.is_infinite( e1 ) ) continue;
  // 	  if ( t.is_infinite( e2 ) ) continue;
  // 	  //  A ---- D
  // 	  //  | \    |
  // 	  //  |  \e1 |
  // 	  //  | e2\  |
  // 	  //  |    \ |
  // 	  //  C ---- B
  // 	  VertexHandle A = e1.first->vertex( t.ccw( e1.second ) );
  // 	  VertexHandle B = e1.first->vertex( t.cw( e1.second ) );
  // 	  VertexHandle C = e2.first->vertex( e2.second );
  // 	  VertexHandle D = e1.first->vertex( e1.second );
  // 	  double l1 = ( toDGtal( A->point() ) - toDGtal( B->point() ) ).norm( Z2i::Point::L_2 );
  // 	  double l2 = ( toDGtal( C->point() ) - toDGtal( D->point() ) ).norm( Z2i::Point::L_2 );
  // 	  // trace.info() << "No flip for new length " << l2 << " > " << l1 << std::endl;
  // 	  if ( ( l2 < l1 ) && isQuadrilateral( t, A, C, B, D ) )
  // 	    {
  // 	      trace.info() << "Flip for new length " << l2 << " < " << l1 << std::endl;
  // 	      t.flip( e1.first, e1.second );
  // 	      changed = true;
  // 	    }
  // 	}
  //   }
 
  
  for( Edge_iterator it = t.edges_begin(), itend = t.edges_end();
       it != itend; ++it )
    {
      Edge e1 = *it;
      Edge e2 = t.mirror_edge( e1 );
      bool empty_e1 = emptyLatticeTriangle( t, e1.first );
      bool empty_e2 = emptyLatticeTriangle( t, e2.first );
      if ( ( empty_e1 && ! empty_e2 )
	   || ( ! empty_e1 && empty_e2 ) )
	{
	  board.setPenColor(DGtal::Color::Blue);
          board.setFillColor( DGtal::Color::None );
          board.setLineWidth( 3.0 );
	  Z2i::Point a( toDGtal( e1.first->vertex( t.ccw( e1.second ) )->point() ) );
	  Z2i::Point b( toDGtal( e1.first->vertex( t.cw( e1.second ) )->point() ) );
          board.drawLine(a[0],a[1],b[0],b[1]);
	}
    }
  
  for(Faces_iterator it = t.finite_faces_begin(), itend=t.finite_faces_end();
      it != itend; ++it)
    {
      Z2i::Point a( toDGtal(it->vertex(0)->point())),
	b(toDGtal(it->vertex(1)->point())),
	c(toDGtal(it->vertex(2)->point()));

      // Z2i::Vector ab( b - a ), ac( c - a );
      // int d = ab[ 0 ] * ac[ 1 ] - ab[ 1 ] * ac[ 0 ];
      if ( emptyLatticeTriangle( t, it ) ) //( ( d == 1 ) || (d == -1 ) )
        {
          board.setPenColor(DGtal::Color::Green);
          board.setFillColor( DGtal::Color::None );
          board.setLineWidth( 1.0 );
          board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
        }
      else
        {
          board.setPenColor(DGtal::Color::Red);
          board.setFillColor( DGtal::Color::None );
          //          board.setFillColorRGBi(200,200,200,128);
          board.setLineWidth( 1.0 );
          board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
        }
    }

  // Segmentation::SegmentComputerIterator i = theSegmentation.begin();
  // Segmentation::SegmentComputerIterator end = theSegmentation.end();
  // board.setPenColor(DGtal::Color::Green);
  // board.setFillColor( DGtal::Color::None );
  // board << SetMode( "ArithmeticalDSS", "BoundingBox" );
  // std::string aStyleName = "ArithmeticalDSS/BoundingBox";
  // for ( ; i != end; ++i) {
  //   DSS4 current(*i);
  //   board << CustomStyle( aStyleName, 
  //                         new CustomPenColor( DGtal::Color::Magenta ) )
  //         << current;
  // } 

  // Display Voronoi.
  // for(Edge_iterator it = t.edges_begin(), itend=t.edges_end();
  //     it != itend; ++it)
  //   {
  //     // vertex(cw(i)) and vertex(ccw(i)) of f.
  //     Face_handle itf = it->first;
  //     int i = it->second;
  //     Z2i::Point a( toDGtal(itf->vertex( t.cw( i ) )->point()));
  //     Z2i::Point b( toDGtal(itf->vertex( t.ccw( i ) )->point()));

  //     CGAL::Object o = t.dual( it );
  //     if (CGAL::object_cast<K::Segment_2>(&o)) 
  //       {
  //         const K::Segment_2* ptrSegment = CGAL::object_cast<K::Segment_2>(&o);
  //         board.setPenColor(DGtal::Color::Black);
  //         board.setFillColor( DGtal::Color::None );
  //         board.setLineWidth( 2.0 );
  //         board.drawLine( ptrSegment->source().x(),
  //                         ptrSegment->source().y(),
  //                         ptrSegment->target().x(),
  //                         ptrSegment->target().y() );
  //       }
  //     else if (CGAL::object_cast<K::Ray_2>(&o)) 
  //       {
  //         const K::Ray_2* ptrRay = CGAL::object_cast<K::Ray_2>(&o);
  //         board.setPenColor(DGtal::Color::Black);
  //         board.setFillColor( DGtal::Color::None );
  //         board.setLineWidth( 2.0 );
  //         double dx = ptrRay->to_vector().x();
  //         double dy = ptrRay->to_vector().y();
  //         double norm = sqrt( dx*dx+dy*dy );
  //         dx = 5.0 * dx / norm;
  //         dy = 5.0 * dy / norm;
  //         board.drawArrow( ptrRay->source().x(),
  //                          ptrRay->source().y(),
  //                          ptrRay->source().x() + dx, //1*ptrRay->to_vector().x(),
  //                          ptrRay->source().y() + dy ); //1*ptrRay->to_vector().y() );
  //       }
  //   }

  
  board.saveSVG("delaunay.svg");
  board.saveEPS("delaunay.eps");
  
  for ( DigitalCore::FiniteVerticesIterator it = core.triangulation().finite_vertices_begin(),
          ite = core.triangulation().finite_vertices_end(); it != ite; ++it )
    std::cout << it->point() << std::endl;

  return 0;
}
