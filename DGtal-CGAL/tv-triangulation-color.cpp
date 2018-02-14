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
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/shapes/TriangulatedSurface.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/GradientColorMap.h>
#include <DGtal/io/colormaps/GrayscaleColorMap.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/PPMWriter.h"
#include "CairoViewer.h"
#include <DGtal/geometry/helpers/ContourHelper.h>
#include "BasicVectoImageExporter.h"
#include "ImageConnecter.h"


// #include <CGAL/Delaunay_triangulation_2.h>
// #include <CGAL/Constrained_Delaunay_triangulation_2.h>
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Triangulation_2.h>

#include "cairo.h"

// #include "Triangulation2DHelper.h"
// #include "UmbrellaPart2D.h"
// #include "Auxiliary.h"
///////////////////////////////////////////////////////////////////////////////
double randomUniform()
{
  return (double) random() / (double) RAND_MAX;
}

namespace DGtal {
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
  
  struct TVTriangulation
  {

    typedef Z2i::Integer               Integer;
    typedef Z2i::RealPoint             Point;
    typedef Z2i::RealVector            Vector;
    typedef Z2i::RealPoint             RealPoint;
    typedef Z2i::RealVector            RealVector;
    typedef Z2i::Domain                Domain;
    typedef TriangulatedSurface<Point> Triangulation;
    typedef Triangulation::VertexIndex VertexIndex;
    typedef Triangulation::Vertex      Vertex;
    typedef Triangulation::Arc         Arc;
    typedef Triangulation::Face        Face;
    typedef Triangulation::VertexRange VertexRange;
    typedef Triangulation::FaceRange   FaceRange;
    typedef double                     Scalar;
    typedef PointVector< 3, Scalar >   Value;
    typedef std::pair<DGtal::Color, std::vector<std::vector<TVTriangulation::Point> > > ColorContours;

    // typedef std::array< Value, 2 >     VectorValue;
    struct VectorValue {
      Value x;
      Value y;
    };
    typedef std::vector<Scalar>        ScalarForm;
    typedef std::vector<Value>         ValueForm;
    typedef std::vector<VectorValue>   VectorValueForm;

    /// The domain triangulation
    Triangulation        T;
    /// The image values at each vertex
    ValueForm            _I;
    /// Initial vertices
    Integer              _nbV;
    /// Tells if it is a color image (just for optimization).
    bool                 _color;
    /// Power for gradient computation
    Scalar               _power;
    /// List of arcs whose energy may have been changed by a surrounding flip.
    std::vector<Arc>     _Queue;
    /// List of arcs that may be flipped with the same energy.
    std::vector<Arc>     _Q_equal;

    /// The TV-regularized values
    ValueForm            _u;
    /// The TV-regularized vectors
    VectorValueForm      _p;
    /// The vector storing the tv energy of each triangle.
    ScalarForm           _tv_per_triangle;
    /// The total variation energy of T.
    Scalar               _tv_energy;
    /// true if arc is flippable.
    std::vector<bool>    _flippable;
    
    // Data needed for vectorization. Each contour is a succession of
    // arc, where head points outside.
    
    /// Contains the relative positions of the contour at each arc.
    /// A value between 0 and 1: 0 is tail and 1 is head of arc
    std::vector<Scalar>     _t;
    /// Contains the coordinates of the barycenter of each face.
    std::vector<RealPoint>  _b;
    /// Contains the area for each vertex.
    std::vector<Scalar>     _A; 
    
    /// @return the regularized value at vertex v.
    const Value& u( const VertexIndex v ) const
    { return _u[ v ]; }

    /// @return the regularized barycenter at face f.
    const RealPoint& barycenter( const Face f ) const
    { return _b[ f ]; }
    /// @return the contour point at arc \a a.
    RealPoint contourPoint( const Arc a ) const
    {
      RealPoint B = T.position( T.head( a ) );
      RealPoint A = T.position( T.tail( a ) );
      return ( 1.0 - _t[ a ] ) * A + _t[ a ] * B;
    }
    
    /// invalidate a vertex by using specific value (to process image border).
    void invalidate( const VertexIndex v )
    { _u[ v ][0]=std::numeric_limits<double>::min();
      _u[ v ][1]=std::numeric_limits<double>::min();
      _u[ v ][2]=std::numeric_limits<double>::min();
    }
    bool isinvalid(const VertexIndex v)
    {
      return   _u[ v ][0]==std::numeric_limits<double>::min() &&
               _u[ v ][1]==std::numeric_limits<double>::min() &&
               _u[ v ][2]==std::numeric_limits<double>::min();
    }
      
    /// The norm used for the scalars induced by vector-value space (RGB)
    std::function< Scalar( const Value& v ) >       _normX;
    /// The norm used for the 2d vectors induced by vector-value space (RGB)
    std::function< Scalar( const VectorValue& v ) > _normY;
    
    static Scalar square( Scalar x ) { return x*x; }

    // ------------- Discrete operators ---------------------

    /// @return the vector of the norms of each vector in p (one per triangle). 
    ScalarForm norm( const VectorValueForm& p ) const
    {
      ScalarForm S( T.nbFaces() );
      for ( Face f = 0; f < T.nbFaces(); f++ )
	S[ f ] = _normY( p[ f ] );
      return S;
    }
    // Definition of a global gradient operator that assigns vectors to triangles
    VectorValueForm grad( const ValueForm& u ) const
    { // it suffices to traverse all (valid) triangles.
      VectorValueForm G( T.nbFaces() );
      for ( Face f = 0; f < T.nbFaces(); ++f ) {
	G[ f ] = grad( f, u );
      }
      return G;
    }

    // Definition of a (local) gradient operator that assigns vectors to triangles
    VectorValue grad( Face f, const ValueForm& u ) const
    {
      VertexRange V = T.verticesAroundFace( f );
      return grad( V[ 0 ], V[ 1 ], V[ 2 ], u );
    }

    // Definition of a (local) gradient operator that assigns vectors to triangles
    VectorValue grad( VertexIndex i, VertexIndex j, VertexIndex k,
		      const ValueForm& u ) const
    {
      // [ yj-yk yk-yi yi-yk ] * [ ui ]
      // [ xk-xj xi-xk xj-xi ]   [ uj ]
      //                         [ uk ]
      const Point& pi = T.position( i );
      const Point& pj = T.position( j );
      const Point& pk = T.position( k );
      const Value& ui = u[ i ];
      const Value& uj = u[ j ];
      const Value& uk = u[ k ];
      VectorValue G;
      const int d = _color ? 3 : 1;
      for ( int m = 0; m < d; ++m ) {
	G.x[ m ] = ui[ m ] * ( pj[ 1 ] - pk[ 1 ] )
	  +        uj[ m ] * ( pk[ 1 ] - pi[ 1 ] )
	  +        uk[ m ] * ( pi[ 1 ] - pj[ 1 ] );
	G.y[ m ] = ui[ m ] * ( pk[ 0 ] - pj[ 0 ] )
	  +        uj[ m ] * ( pi[ 0 ] - pk[ 0 ] )
	  +        uk[ m ] * ( pj[ 0 ] - pi[ 0 ] );
      }
      G.x *= 0.5;
      G.y *= 0.5;
      return G;
    }

    // Definition of a (global) divergence operator that assigns
    // scalars to vertices from a vector field.
    ValueForm div( const VectorValueForm& G ) const
    {
      const int d = _color ? 3 : 1;
      ValueForm S( T.nbVertices() );
      for ( VertexIndex v = 0; v < T.nbVertices(); ++v )
	S[ v ] = Value{ 0, 0, 0 };
      for ( Face f = 0; f < T.nbFaces(); ++f ) {
	auto V          = T.verticesAroundFace( f );
	const Point& pi = T.position( V[ 0 ] );
	const Point& pj = T.position( V[ 1 ] );
	const Point& pk = T.position( V[ 2 ] );
	const VectorValue& Gf = G[ f ];
	for ( int m = 0; m < d; ++m ) {
	  S[ V[ 0 ] ][ m ] -= ( pj[ 1 ] - pk[ 1 ] ) * Gf.x[ m ]
	    +                 ( pk[ 0 ] - pj[ 0 ] ) * Gf.y[ m ];
	  S[ V[ 1 ] ][ m ] -= ( pk[ 1 ] - pi[ 1 ] ) * Gf.x[ m ]
	    +                 ( pi[ 0 ] - pk[ 0 ] ) * Gf.y[ m ];
	  S[ V[ 2 ] ][ m ] -= ( pi[ 1 ] - pj[ 1 ] ) * Gf.x[ m ]
	    +                 ( pj[ 0 ] - pi[ 0 ] ) * Gf.y[ m ];
	}
      }
      for ( VertexIndex v = 0; v < T.nbVertices(); ++v )
	S[ v ] *= 0.5;
      return S;
    }


    /// @return the scalar form lambda.u
    ValueForm multiplication( Scalar lambda, const ValueForm& u ) const
    {
      ValueForm S( T.nbVertices() );
      for ( VertexIndex v = 0; v < T.nbVertices(); ++v )
	S[ v ] = lambda * u[ v ];
      return S;
    }
    
    /// @return the scalar form u - v
    ValueForm subtraction( const ValueForm& u, const ValueForm& v ) const
    {
      ValueForm S( T.nbVertices() );
      for ( VertexIndex i = 0; i < T.nbVertices(); ++i )
	  S[ i ] = u[ i ] - v[ i ];
      return S;
    }
    /// u -= v
    void subtract( ValueForm& u, const ValueForm& v ) const
    {
      for ( VertexIndex i = 0; i < T.nbVertices(); ++i )
	u[ i ] -= v[ i ];
    }

    /// @return the scalar form a.u + b.v
    ValueForm combination( const Scalar a, const ValueForm& u,
			   const Scalar b, const ValueForm& v ) const
    {
      ValueForm S( T.nbVertices() );
      trace.info() << "[combination] variation is ["
		   << ( b * *( std::min_element( v.begin(), v.end() ) ) )
		   << " " << ( b * *( std::max_element( v.begin(), v.end() ) ) )
		   << std::endl;
      for ( VertexIndex i = 0; i < T.nbVertices(); ++i )
	S[ i ] = a * u[ i ] + b * v[ i ];
      return S;
    }
      
    /// The function that evaluates the energy at each triangle.
    /// It is now just the norm of the gradient.
    Scalar computeEnergyTV( VertexIndex v1, VertexIndex v2, VertexIndex v3 ) const
    {
      return _normY( grad( v1, v2, v3, _u ) );
    }

    /// @return the tv energy stored at this face.
    Scalar computeEnergyTV( const Face f )
    {
      VertexRange V = T.verticesAroundFace( f );
      return ( _tv_per_triangle[ f ] = computeEnergyTV( V[ 0 ], V[ 1 ], V[ 2 ] ));
    }
    
    /// @return the tv energy stored at this face.
    Scalar energyTV( const Face f ) const
    {
      return _tv_per_triangle[ f ];
    }

    /// @return the tv energy stored at this face.

    Scalar& energyTV( const Face f )
    {
      return _tv_per_triangle[ f ];
    }

    /// Compute (and store in _tv_per_triangle) the TV-energy per triangle.
    Scalar computeEnergyTV()
    {
      Scalar E = 0;
      for ( Face f = 0; f < T.nbFaces(); ++f )	{
	E += computeEnergyTV( f );
      }
      _tv_energy = E;
      // trace.info() << "TV(u) = " << E << std::endl;
      return E;
    }

    /// Gets the current TV energy of the triangulation.
    Scalar getEnergyTV()
    {
      return _tv_energy;
    }

    /// @return the aspect ratio of a face (the greater, the most elongated it is.
    Scalar aspectRatio( const Face f ) const
    {
      VertexRange P  = T.verticesAroundFace( f );
      const Point& a = T.position( P[ 0 ] );
      const Point& b = T.position( P[ 1 ] );
      const Point& c = T.position( P[ 2 ] );
      Vector     ab = b - a;
      Vector     bc = c - b;
      Vector     ca = a - c;
      Scalar    dab = ab.norm();
      Scalar    dbc = bc.norm();
      Scalar    dca = ca.norm();
      Vector    uab = ab / dab;
      Vector    ubc = bc / dbc;
      Vector    uca = ca / dca;
      Scalar     ha = ( ab - ab.dot( ubc ) * ubc  ).norm();
      Scalar     hb = ( bc - bc.dot( uca ) * uca  ).norm();
      Scalar     hc = ( ca - ca.dot( uab ) * uab  ).norm();
      return std::max( dab / hc, std::max( dbc / ha, dca / hb ) );
    }

    /// @return the diameter of a face (the greater, the most elongated it is.
    Scalar diameter( const Face f ) const
    {
      VertexRange P  = T.verticesAroundFace( f );
      const Point& a = T.position( P[ 0 ] );
      const Point& b = T.position( P[ 1 ] );
      const Point& c = T.position( P[ 2 ] );
      Vector     ab = b - a;
      Vector     bc = c - b;
      Vector     ca = a - c;
      Scalar    dab = ab.norm();
      Scalar    dbc = bc.norm();
      Scalar    dca = ca.norm();
      return std::max( dab, std::max( dbc, dca ) );
    }

    // -------------- Construction services -------------------------
    
    // Constructor from color image.
    template <typename Image>
    TVTriangulation( const Image&  I, bool color,
		     Scalar p = 0.5, Scalar sim = 0.0 )
    {
      _color = color;
      _power = p;
      // Defining norms.
      if ( color ) {
	_normX = [p] ( const Value& v ) -> Scalar
	  {
	    return pow( square( v[ 0 ] ) + square( v[ 1 ] ) + square( v[ 2 ] ), p);
	  };
	_normY = [p] ( const VectorValue& v ) -> Scalar
	  {
	    return pow( square( v.x[ 0 ] ) + square( v.y[ 0 ] ), p )
	    + pow( square( v.x[ 1 ] ) + square( v.y[ 1 ] ), p )
	    + pow( square( v.x[ 2 ] ) + square( v.y[ 2 ] ), p );
	  };
	// // Standard ColorTV is
	// _normY = [p] ( const VectorValue& v ) -> Scalar
	// {
	//   return pow( square( v.x[ 0 ] ) + square( v.y[ 0 ] )
	// 		+ square( v.x[ 1 ] ) + square( v.y[ 1 ] )
	// 		+ square( v.x[ 2 ] ) + square( v.y[ 2 ] ), p );
	// };
      } else {
	_normX = [p] ( const Value& v ) -> Scalar
	  {
	    return pow( square( v[ 0 ] ), p );
	  };
	_normY = [p] ( const VectorValue& v ) -> Scalar
	  {
	    return pow( square( v.x[ 0 ] ) + square( v.y[ 0 ] ), p );
	  };
      }
      // Creates image form _I
      typedef std::function< int( int ) > ColorConverter;
      ColorConverter converters[ 4 ];
      converters[ 0 ] = ColorToRedFunctor();
      converters[ 1 ] = ColorToGreenFunctor();
      converters[ 2 ] = ColorToBlueFunctor();
      converters[ 3 ] = GrayToGrayFunctor();
      int   red = color ? 0 : 3;
      int green = color ? 1 : 3;
      int  blue = color ? 2 : 3;
      VertexIndex v = 0;
      _I.resize( I.size() );
      for ( unsigned int val : I ) {
	_I[ v++ ] = Value( (Scalar) converters[ red ]  ( val ),
			   (Scalar) converters[ green ]( val ),
			   (Scalar) converters[ blue ] ( val ) );
      }
      
      // Building connections.
      typedef ImageContainerBySTLVector<Domain,Value> ValueImage;
      typedef ImageConnecter<ValueImage> Connecter;
      ValueImage tmpI( I.domain() );
      auto it = tmpI.begin();
      for ( auto v : _I ) *it++ = v;

      Connecter connecter;
      typename Connecter::Comparator comp
	= [this] ( Value v1, Value v2 ) { return _normX( v1 - v2 ); };
      trace.info() << "Compute connections ... ";
      connecter.init2( tmpI, comp, sim );
      trace.info() << "ended." << std::endl;

      // Building triangulation
      const Point taille = I.extent();
      // Creates vertices
      for ( auto p : I.domain() ) T.addVertex( p );
      // Creates triangles
      for ( Integer y = 0; y < taille[ 1 ] - 1; ++y ) {
	for ( Integer x = 0; x < taille[ 0 ] - 1; ++x ) {
	  const VertexIndex v00 = y * taille[ 0 ] + x;
	  const VertexIndex v10 = v00 + 1;
	  const VertexIndex v01 = v00 + taille[ 0 ];
	  const VertexIndex v11 = v01 + 1;
	  auto  how = connecter.howConnected( Point( x, y ) );
	  bool diag00_11 = ( how.diagonal == Connecter::Diagonal00_11 );
	  if ( diag00_11 ) {
	    T.addTriangle( v00, v01, v11 );
	    T.addTriangle( v00, v11, v10 );
	  } else {
	    T.addTriangle( v00, v01, v10 );
	    T.addTriangle( v10, v01, v11 );
	  }
	}
      }
      bool ok = T.build();
      trace.info() << "Build triangulation: "
		   << ( ok ? "OK" : "ERROR" ) << std::endl;
      _nbV = T.nbVertices();
      // Building forms.
      _u = _I;                  // u = image at initialization
      _p.resize( T.nbFaces() ); // p = 0     at initialization
      // TV-energy is computed and stored per face to speed-up computations.
      _tv_per_triangle.resize( T.nbFaces() );
      computeEnergyTV();

      // Fix some arcs;
      _flippable.resize( T.nbArcs() );
      int nbFlippable = 0;
      for ( Arc a = 0; a < T.nbArcs(); ++a ) {
	Point p = T.position( T.head( a ) );
	Point q = T.position( T.tail( a ) );
	Point l = p.inf( q );
	Point u = p.sup( q );
	auto  how = connecter.howConnected( l );
	if ( ( l - u ).dot( l - u ) == 2 )
	  _flippable[ a ] = ( how.diagonal == Connecter::Default );
	else if ( u[ 0 ] != l[ 0 ] )
	  _flippable[ a ] = ! how.horizontal;
	else 
	  _flippable[ a ] = ! how.vertical;
	nbFlippable += _flippable[ a ] ? 1 : 0;
      }
      trace.info() << "Nb arcs flippable = " << nbFlippable
		   << "/" << T.nbArcs() << std::endl;
    }

    template <typename Image>
    bool outputU( Image& J ) const
    {
      VertexIndex v = 0;
      for ( unsigned int & val : J ) {
	val = _color
	  ? ( ( (int) _u[ v ][ 0 ] ) << 16 )
	  + ( ( (int) _u[ v ][ 1 ] ) << 8 )
	  + ( (int) _u[ v ][ 2 ] )
	  : ( ( (int) _u[ v ][ 0 ] ) << 16 )
	  + ( ( (int) _u[ v ][ 0 ] ) << 8 )
	  + ( (int) _u[ v ][ 0 ] );
	v += 1;
      }
      return v == T.nbVertices();
    }

    // @return the determinant of \a pq and \a qr.
    static Scalar det( const Point& pq, const Point& qr )
    {
      return pq[ 0 ] * qr[ 1 ] - pq[ 1 ] * qr[ 0 ];
    }

    static Scalar doesTurnLeft( const Point& p, const Point& q, const Point& r )
    {
      const Point pq = q - p;
      const Point qr = r - q;
      return det( pq, qr );
    }
    static Scalar doesTurnLeft( const Point& pq, const Point& qr )
    {
      return det( pq, qr );
    }
    
    // Check strict convexity of quadrilateron.
    bool isConvex( const VertexRange& V ) const
    {
      Point P[] = { T.position( V[ 1 ] ) - T.position( V[ 0 ] ),
		    T.position( V[ 2 ] ) - T.position( V[ 1 ] ),
		    T.position( V[ 3 ] ) - T.position( V[ 2 ] ),
		    T.position( V[ 0 ] ) - T.position( V[ 3 ] ) };
      bool cvx = ( doesTurnLeft( P[ 0 ], P[ 1 ] ) < 0 )
	&&       ( doesTurnLeft( P[ 1 ], P[ 2 ] ) < 0 )
	&&       ( doesTurnLeft( P[ 2 ], P[ 3 ] ) < 0 )
	&&       ( doesTurnLeft( P[ 3 ], P[ 0 ] ) < 0 );
      return cvx;
    }
    /**
       NB: process only arcs (s,t) where ( t > s ).
       
       @return 1 is energy is lowered (a has been flipped), 0 if arc
       is flippable but it does not change the energy, negative
       otherwise (-1: boundary, -2 s > t, -3 non convex, -4 increases
       the energy.
    */
    int updateArc( const Arc a ) {
      // Checks that edge can be flipped.
      if ( ! _flippable[ a ] ) return -4;
      VertexRange P = T.verticesAroundArc( a );
      if ( P.size() != 4 )   return -1;
      if ( P[ 0 ] < P[ 2 ] ) return -2;
      if ( ! isConvex( P ) ) return -3;
      // Computes energies
      const Face    f012 = T.faceAroundArc( a );
      const Face    f023 = T.faceAroundArc( T.opposite( a ) );
      const Scalar  E012 = energyTV( f012 ); //P[ 0 ], P[ 1 ], P[ 2 ] );
      const Scalar  E023 = energyTV( f023 ); //P[ 0 ], P[ 2 ], P[ 3 ] );
      const Scalar  E013 = computeEnergyTV( P[ 0 ], P[ 1 ], P[ 3 ] );
      const Scalar  E123 = computeEnergyTV( P[ 1 ], P[ 2 ], P[ 3 ] );
      const Scalar Ecurr = E012 + E023;
      const Scalar Eflip = E013 + E123;
      if ( Eflip < Ecurr ) {
	// Save arcs that may be affected.
	queueSurroundingArcs( a );
	T.flip( a );
	_tv_per_triangle[ f012 ] = E123; // f012 -> f123
	_tv_per_triangle[ f023 ] = E013; // f023 -> f013
	_tv_energy += Eflip - Ecurr;
	return 1;
      } else if ( Eflip == Ecurr ) {
	return ( Eflip > 0.0 ) ? 0 : -6;
      } else return -7;
    }

    void queueSurroundingArcs( const Arc a )
    {
      Arc around[ 4 ];
      around[ 0 ] = T.next( a );
      around[ 1 ] = T.next( around[ 0 ] );
      around[ 2 ] = T.next( T.opposite( a ) );
      around[ 3 ] = T.next( around[ 2 ] );
      for ( int i = 0; i < 4; ++i ) {
	_Queue.push_back( around[ i ] );
	_Queue.push_back( T.opposite( around[ i ] ) );
      }
    }

    /// Quantify the regularized image _u.
    void quantify( const int level )
    {
      const Scalar factor = 255.0 / ( level - 1);  
      for ( VertexIndex i = 0; i < T.nbVertices(); ++i )
	for ( int m = 0; m < 3; ++m ) {
	  _u[ i ][ m ] = round( ( _u[ i ][ m ] ) / factor ) * factor;
	  _u[ i ][ m ] = std::min( 255.0, std::max( 0.0, _u[ i ][ m ] ) );
	}
      computeEnergyTV();
    }
    
    /// Does one pass of TV regularization (u, p and I must have the
    /// meaning of the previous iteration).
    Scalar tvPass( Scalar lambda, Scalar dt, Scalar tol, int N = 10 )
    {
      trace.info() << "TV( u ) = " << getEnergyTV() << std::endl;
      //trace.info() << "lambda.f" << std::endl;
      VectorValueForm p( _p.size() ); // this is p^{n+1}
      ValueForm      lf = multiplication( lambda, _I );
      Scalar     diff_p = 0.0;
      int             n = 0; // iteration number
      do {
	// trace.info() << "div( p ) - lambda.f" << std::endl;
	ValueForm        dp_lf = subtraction( div( _p ), lf );
	// trace.info() << "G := grad( div( p ) - lambda.f)" << std::endl;
	VectorValueForm gdp_lf = grad( dp_lf );
	// trace.info() << "N := | G |" << std::endl;
	ScalarForm     ngdp_lf = norm( gdp_lf );
	// trace.info() << "p^n+1 := ( p + dt * G ) / ( 1 + dt | G | )" << std::endl;
	diff_p = 0.0;
	for ( Face f = 0; f < T.nbFaces(); f++ ) {
	  Scalar alpha = 1.0 / ( 1.0 + dt * ngdp_lf[ f ] );
	  if ( alpha <= 0.0 )
	    trace.warning() << "Face " << f << " alpha=" << alpha << std::endl;
	  p[ f ].x  = alpha * ( _p[ f ].x + dt * gdp_lf[ f ].x );
	  p[ f ].y  = alpha * ( _p[ f ].y + dt * gdp_lf[ f ].y );
	  VectorValue delta = { p[ f ].x - _p[ f ].x,
				p[ f ].y - _p[ f ].y };
	  diff_p = std::max( diff_p, _normY( delta ) );
	}
	trace.info() << "Iter n=" << (n++) << " diff_p=" << diff_p
		     << " tol=" << tol << std::endl;
	std::swap( p, _p );
      } while ( ( diff_p > tol ) && ( n < N ) );
      _u = combination( 1.0, _I, -1.0/lambda, div( _p ) );
      if ( ! _color ) {
	for ( VertexIndex i = 0; i < _u.size(); ++i )
	  _u[ i ][ 2 ] = _u[ i ][ 1 ] = _u[ i ][ 0 ];
      }
      trace.info() << "TV( u ) = " << computeEnergyTV() << std::endl;
      return diff_p;
    }
    
    // equal_strategy:
    // 0: do nothing
    // 1: subdivide all
    // 2: flip all
    // 3: flip all only if flipped = 0
    // @return either (nbflip,0), (nbflip, nbsub_eq) or (nbflip, nbflip_eq).
    std::pair<Integer,Integer>
    onePass( Scalar& total_energy, int equal_strategy = 0 )
    {
      Integer nbflipped = 0;
      Integer   nbequal = 0;
      total_energy      = 0;
      std::vector<Arc> Q_process;
      std::swap( _Queue, Q_process );
      // Taking care of first pass
      if ( Q_process.size() == 0 )
	for ( Arc a = 0; a < T.nbArcs(); ++a )
	  if ( _flippable[ a ] ) Q_process.push_back( a );
      // Processing arcs
      for ( Arc a : Q_process ) {
	int update = updateArc( a );
	if ( update > 0 ) nbflipped++;
	else if ( update == 0 ) _Q_equal.push_back( a );
      }
      total_energy = getEnergyTV();
      trace.info() << "TV( u ) = " << total_energy
		   << " nbflipped=" << nbflipped
		   << "/" << Q_process.size();
      if ( equal_strategy == 1 ) {
	nbequal = subdivide( _Q_equal );
	trace.info() << " nbsubequal=" << nbequal
		     << "/" << _Q_equal.size();
	_Q_equal.clear();
      } else if ( equal_strategy == 2 ) {
	nbequal = flipEqual( _Q_equal );
	trace.info() << " nbflipequal=" << nbequal
		     << "/" << _Q_equal.size();
	_Q_equal.clear();
      } else if ( ( equal_strategy == 3 ) && ( nbflipped == 0 ) ) {
	nbequal = flipEqual( _Q_equal );
	trace.info() << " nbflipequal=" << nbequal
		     << "/" << _Q_equal.size();
	_Q_equal.clear();
      } else if ( ( ( equal_strategy == 4 ) || ( equal_strategy == 5 ) )
		  && ( nbflipped == 0 ) ) {
      nbequal = flipEqualWithProb( _Q_equal, 0.5 );
	trace.info() << " nbflipequalP=" << nbequal
		     << "/" << _Q_equal.size();
	_Q_equal.clear();
      }
      trace.info() << std::endl;
      return std::make_pair( nbflipped, nbequal );
    }
    
    template <typename Range>
    Integer flipEqual( const Range& range )
    {
      Integer nbflip = 0;
      for ( Arc a : range ) {
	int update = updateArc( a );
	if ( update == 0 ) {
	  // Save arcs that may be affected.
	  queueSurroundingArcs( a );
	  T.flip( a );
	  nbflip++;
	}
      }
      return nbflip;
    }
    
    template <typename Range>
    Integer flipEqualWithProb( const Range& range, double p )
    {
      Integer nbflip = 0;
      for ( Arc a : range ) {
	int update = updateArc( a );
	if ( update == 0 ) {
	  // Put arc back to potentially process it again.
	  _Queue.push_back( a );
	  _Queue.push_back( T.opposite( a ) );
	  if ( randomUniform() < p ) {
	    // Save arcs that may be affected.
	    queueSurroundingArcs( a );
	    T.flip( a );
	    nbflip++;
	  } 
	}
      }
      return nbflip;
    }
    
    template <typename Range>
    Integer subdivide( const Range& range )
    {
      Scalar  energy = 0.0;
      Integer nbsubdivided = 0;
      for ( Arc a : range ) {
	int update = updateArc( a );
	if ( update == 0 ) {
	  VertexRange P = T.verticesAroundArc( a );
	  // Allow one level of subdivision.
	  if ( std::max( std::max( P[ 0 ], P[ 1 ] ),
			 std::max( P[ 2 ], P[ 3 ] ) ) >= _nbV ) continue;
	  // Save arcs that may be affected.
	  queueSurroundingArcs( a );
	  // Remember faces.
	  const Face    f012 = T.faceAroundArc( a );
	  const Face    f023 = T.faceAroundArc( T.opposite( a ) );
	  Scalar     Ebefore = energyTV( f012 ) + energyTV( f023 );
	  Scalar      Eafter = 0.0; 
	  Point B = ( T.position( P[ 0 ] ) + T.position( P[ 1 ] )
		      + T.position( P[ 2 ] ) + T.position( P[ 3 ] ) ) * 0.25;
	  VertexIndex v = T.split( a, B );
	  FaceRange   F = T.facesAroundVertex( v );
	  Face    new_f = _tv_per_triangle.size();
	  _tv_per_triangle.resize( new_f + 2 );
	  _p.resize( new_f + 2 );
	  for ( Face f : F )
	    Eafter    += computeEnergyTV( f );
	  _tv_energy   += Eafter - Ebefore;
	  
	  Value V = ( _u[ P[ 0 ] ] + _u[ P[ 1 ] ]
		      + _u[ P[ 2 ] ] + _u[ P[ 3 ] ] ) * 0.25;
	  Value VI = ( _I[ P[ 0 ] ] + _I[ P[ 1 ] ]
		       + _I[ P[ 2 ] ] + _I[ P[ 3 ] ] ) * 0.25;
	  _u.push_back( V );
	  _I.push_back( VI );
	  ++nbsubdivided;
	}
      }
      return nbsubdivided;
    }

    // -------------------------- Regularization services --------------------
  public:

    /// Initializes the regularization process for vectorization.
    void regularizeContours( Scalar max_dt = 0.001, int max_iter = 10 )
    {
      initContours();
      for ( int n = 0; n < max_iter; ++n ) {
	std::vector<Scalar> former_t = _t;
	updateBarycenters();
	updateContours();
	Scalar dt = enforceContoursArea( former_t );
	if ( dt < max_dt ) break;
	std::cout << n << " dt = " << dt << " max_dt = " << max_dt << std::endl;
      }
    }

    /// Initializes the regularization process for vectorization.
    void initContours()
    {
      // Computes the barycenter for each valid face.
      _b.resize( T.nbFaces() );
      for ( Face f = 0; f < T.nbFaces(); ++f )	{
	VertexRange V = T.verticesAroundFace( f );
	_b[ f ] = ( T.position( V[ 0 ] )
		    + T.position( V[ 1 ] )
		    + T.position( V[ 2 ] ) ) / 3.0;
      }
      // Computes the contour intersections at each arc.
      _t.resize( T.nbArcs() );
      for ( Arc a = 0; a < T.nbArcs(); ++a ) {
	Arc    opp_a = T.opposite( a );
	if ( T.isArcBoundary( a ) || T.isArcBoundary( opp_a ) ) {
	  _t[ a ] = 0.5;
	  continue;
	}
	Face     f_a = T.faceAroundArc( a );
	Face f_opp_a = T.faceAroundArc( opp_a );
	RealPoint  B = T.position( T.head( a ) );
	RealPoint  A = T.position( T.head( opp_a ) );
	auto       I = intersect( A, B, _b[ f_a ], _b[ f_opp_a ] );
	_t[ a ]      = std::min( 0.999, std::max( 0.001, I.first ) );
	// std::cout << "t[" << a << "]=" << _t[ a ] << std::endl;
      }
      // Computes the areas associated with each vertex
      _A.resize( T.nbVertices() );
      for ( Vertex v = 0; v < T.nbVertices(); ++v ) {
	// std::cout << v << " area=" << areaAtVertex( v ) << std::endl;
	// _A[ v ] = areaAtVertex( v );
	_A[ v ] = (Scalar) T.degree( v ) / 6.0;//areaAtVertex( v );
      }
    }

    bool isArcContour( const Arc a ) const {
      return _u[ T.head( a ) ] != _u[ T.tail( a ) ];
    }
    Scalar arcSimilarity( const Arc a ) const {
      Value diff = _u[ T.head( a ) ] - _u[ T.tail( a ) ];
      return diff.dot( diff );
      // return ( _u[ T.head( a ) ] - _u[ T.tail( a ) ] ).norm();
    }

    
    /// Computes the barycenters for each triangle. The idea is to
    /// count only arcs with different values.
    void updateBarycenters()
    {
      for ( Face f = 0; f < T.nbFaces(); ++f )	{
	VertexRange V = T.verticesAroundFace( f );
	Scalar      w = 0.0;
	Point       B = Point::zero;
	// Point       P[ 3 ] = { T.position( V[ 0 ] ), T.position( V[ 1 ] ),
	// 		       T.position( V[ 2 ] ) };
	Arc arcs[ 3 ];
	Scalar s[ 3 ];
	// Sort similarities...
	arcs[ 0 ] = T.arc( V[ 0 ], V[ 1 ] );
	arcs[ 1 ] = T.next( arcs[ 0 ] );
	arcs[ 2 ] = T.next( arcs[ 1 ] );
	for ( int i = 0; i < 3; ++i ) {
	  s[ i ] = arcSimilarity( arcs[ i ] );
	}
	int  m  = std::min_element( s, s + 3 ) - s;
	bool id = ( ( 2*s[ m ] ) < s[ (m+1)%3 ] )
	  &&      ( ( 2*s[ m ] ) < s[ (m+2)%3 ] );
	if ( id ) s[ m ] = 0;
	for ( int i = 0; i < 3; ++i ) {
	  B += s[ i ] * contourPoint ( arcs[ i ] );
	  w += s[ i ];
	}
	if ( w  > 0.0 ) _b[ f ] += 0.5 * ( B / w - _b[ f ] );
	// Arc    a = T.arc( V[ 0 ], V[ 1 ] );
	// Scalar s = arcSimilarity( a );
	// B       += s * ( ( 1 - _t[ a ] ) * P[ 0 ] + _t[ a ] * P[ 1 ] );
	// w       += s;
	// a        = T.next( a );
	// s        = arcSimilarity( a );
	// B       += s * ( ( 1 - _t[ a ] ) * P[ 1 ] + _t[ a ] * P[ 2 ] );
	// w       += s;
	// a        = T.next( a );
	// s        = arcSimilarity( a );
	// B       += s * ( ( 1 - _t[ a ] ) * P[ 2 ] + _t[ a ] * P[ 0 ] );
	// w       += s;
	// if ( w  > 0.0 ) _b[ f ] += 0.5 * ( B / w - _b[ f ] );
      }
    }

    /// Computes the contour points for each arc from the barycenters.
    /// @return the maximum displacement.
    void updateContours()
    {
      for ( Arc a = 0; a < T.nbArcs(); ++a ) {
	Arc    opp_a = T.opposite( a );
	// if ( T.isArcBoundary( a ) || T.isArcBoundary( opp_a ) ) continue;
	Face     f_a = T.faceAroundArc( a );
	Face f_opp_a = T.faceAroundArc( opp_a );
	RealPoint  B = T.position( T.head( a ) );
	RealPoint  A = T.position( T.head( opp_a ) );
	if ( ( f_a != T.INVALID_FACE ) && ( f_opp_a != T.INVALID_FACE ) ) {
	  auto       I = intersect( A, B, _b[ f_a ], _b[ f_opp_a ] );
	  Scalar     t = std::min( 0.999, std::max( 0.001, I.first ) );
	  _t[ a ]     += 0.5 * ( t - _t[ a ] );
	} else {
	  _t[ a ]      = 0.5;
	}
      }
    }

    /// Moves contour points to keep constant areas
    /// @return the maximum displacement.
    Scalar enforceContoursArea( const std::vector<Scalar>& prev_t )
    {
      Scalar max_t = 0.0;
      // Update t to keep volume constant
      for ( Vertex v = 0; v < T.nbVertices(); ++v ) {
	auto  out_arcs = T.outArcs( v );
	for ( Arc a : out_arcs ) {
	   // Checks the evolution of the arc area.
	  Scalar ratio = areaAtArc( a ) / ( _A[ v ] / T.degree( v ) );
	  if ( ratio <= 0.001 ) {
	    trace.warning() << "Negative ratio " << ratio
			    << " for area[ " << v << " ]=" << _A[ v ]
			    << " at_arc=" << areaAtArc( a )
			    << " d=" << T.degree( v )
			    << " A/d=" << ( _A[ v ] / T.degree( v ) )
			    << std::endl;
	    ratio = 0.001;
	  }
	  // _t[ a ] /= ( 0.2 + 0.8 * ratio );
	  _t[ a ] = 0.5 * _t[ a ] * ( 1.0 + 1.0 / ratio );
	}
	// // Checks the evolution of the vertex area.
	// Scalar ratio = areaAtVertex( v ) / _A[ v ];
	// if ( ratio <= 0.0 )
	//   trace.warning() << "Negative ratio " << ratio
	// 		  << " for area[ " << v << " ]="
	// 		  << _A[ v ] << std::endl;
	// // Rescale all points to recover the correct area.
	// auto  out_arcs = T.outArcs( v );
	// // Todo: not the right formula
	// for ( Arc a : out_arcs ) _t[ a ] /= ( 0.5 + 0.5 * ratio );
      }
      // Averages movements on each edge.
      for ( Arc a = 0; a < T.nbArcs(); ++a ) {
	const Arc opp_a = T.opposite( a );
	// if ( T.isArcBoundary( a ) || T.isArcBoundary( opp_a ) ) continue;
	if ( T.head( a ) < T.head( opp_a ) ) continue;
	Scalar     t = 0.5 * ( _t[ a ] + ( 1.0 - _t[ opp_a ] ) );
	t            = std::min( 0.999, std::max( 0.001, t ) );
	max_t        = std::max( max_t, fabs( t - prev_t[ a ] ) );
	// max_t        = std::max( max_t, fabs( 1.0 - t - _t[ opp_a ] ) );
	_t[ a ]      = t;
	_t[ opp_a ]  = 1.0 - t;
      }
      return max_t;
    }
    
    
    /// The area at an arc \a a is the area associated to the vertex
    /// tail of \a a within the face of \a a. It corresponds to the
    /// zone between the tail, the barycenter and the two intersection
    /// points at the edges containing the tail vertex.
    Scalar areaAtArc( Arc a ) const {
      const Face    f = T.faceAroundArc( a );
      const Arc opp_a = T.opposite( a );
      if ( T.isArcBoundary( a ) || T.isArcBoundary( opp_a ) ) return 1.0/6.0;
      
      const Arc    a2 = T.next( opp_a );
      const Face   f2 = T.faceAroundArc( a2 );
      const Point   B = T.position( T.head( a ) );
      const Point   A = T.position( T.head( opp_a ) );
      const Scalar  t = _t[ a ];
      return -0.5 * ( det( t * ( B - A ), _b[ f ] - A )
		     + det( _b[ f2 ] - A , t * ( B - A ) ) );
      // const Face   f = T.faceAroundArc( a );
      // const Arc   an = T.next( a );
      // const Arc  ann = T.next( an );
      // const Point  B = T.position( T.head( a ) );
      // const Point  C = T.position( T.head( an ) );
      // const Point  A = T.position( T.head( ann ) );
      // const Arc   a2 = T.opposite( ann );
      // const Scalar t = _t[ a ];
      // const Scalar u = _t[ a2 ];
      // return - 0.5 * ( det( t * ( B - A ), _b[ f ] - A )
      // 		       +  det( _b[ f ] - A , u * ( C - A ) ) );
    }

    Scalar areaAtVertex( const Vertex v ) const {
      auto  out_arcs = T.outArcs( v );
      Scalar    area = 0;
      for ( Arc a : out_arcs )
	area += areaAtArc( a );
      return area;
    }

    /// Given two straight lines (AB) and (CD), returns their
    /// intersection point as two scalars (t,u) such that:
    /// I = A + t AB = C + u CD
    static 
    std::pair<Scalar,Scalar> intersect( const Point& A, const Point& B,
					const Point& C, const Point& D )
    {
      Vector AB = B - A;
      Vector DC = C - D;
      Vector AC = C - A;
      Scalar  d = det( AB, DC );
      if ( d == 0 ) return std::make_pair( 0.5, 0.5 );
      Scalar  t = ( DC[ 1 ] * AC[ 0 ] - DC[ 0 ] * AC[ 1 ] ) / d;
      Scalar  u = ( AB[ 0 ] * AC[ 1 ] - AB[ 1 ] * AC[ 0 ] ) / d;
      return std::make_pair( t, u );
    }
  };

  // Useful function for viewing triangulations.

  /**
     This class is intended for visualizing Affine Valued
     triangulation with CAIRO.
  */
  class CairoViewerTV : public CairoViewer<Z2i::Space>
  {
  public:
    typedef CairoViewer<Z2i::Space>  Base;
    typedef TVTriangulation          TVT;
    typedef TVT::Value               Value;
    typedef TVT::Triangulation       Triangulation;
    typedef TVT::VertexIndex         VertexIndex;
    typedef TVT::Vertex              Vertex;
    typedef TVT::Arc                 Arc;
    typedef TVT::Face                Face;
    typedef TVT::VertexRange         VertexRange;
    typedef TVT::Scalar              Scalar;
    typedef TVT::Point               Point;
    typedef Z2i::RealPoint           RealPoint;

  public:

    /**
       Constructor. 
    */
    CairoViewerTV( int x0, int y0, int width, int height, 
		   double xfactor = 1.0, double yfactor = 1.0,
		   int shading = 0,
		   bool color = true,
		   double disc_stiffness = 0.5,
		   double disc_amplitude = 0.75 )
      : Base( x0, y0, width, height, xfactor, yfactor, shading, color,
	      disc_stiffness, disc_amplitude )
    {}
    
    /// Destructor.
    ~CairoViewerTV() {}

    void viewTVTLinearGradientTriangle( TVT & tvT, Face f )
    {
      VertexRange V = tvT.T.verticesAroundFace( f );
      Point a = tvT.T.position( V[ 0 ] );
      Point b = tvT.T.position( V[ 1 ] );
      Point c = tvT.T.position( V[ 2 ] );
      drawLinearGradientTriangle( RealPoint( a[ 0 ], a[ 1 ] ),
				  RealPoint( b[ 0 ], b[ 1 ] ),
				  RealPoint( c[ 0 ], c[ 1 ] ),
				  tvT.u( V[ 0 ] ),
				  tvT.u( V[ 1 ] ),
				  tvT.u( V[ 2 ] ) );
    }
    void viewTVTNonLinearGradientTriangle( TVT & tvT, Face f )
    {
      VertexRange V = tvT.T.verticesAroundFace( f );
      Point a = tvT.T.position( V[ 0 ] );
      Point b = tvT.T.position( V[ 1 ] );
      Point c = tvT.T.position( V[ 2 ] );
      drawNonLinearGradientTriangle( RealPoint( a[ 0 ], a[ 1 ] ),
				     RealPoint( b[ 0 ], b[ 1 ] ),
				     RealPoint( c[ 0 ], c[ 1 ] ),
				     tvT.u( V[ 0 ] ),
				     tvT.u( V[ 1 ] ),
				     tvT.u( V[ 2 ] ) );
    }
    void viewTVTGouraudTriangle( TVT & tvT, Face f )
    {
      VertexRange V = tvT.T.verticesAroundFace( f );
      Point a = tvT.T.position( V[ 0 ] );
      Point b = tvT.T.position( V[ 1 ] );
      Point c = tvT.T.position( V[ 2 ] );
      drawGouraudTriangle( RealPoint( a[ 0 ], a[ 1 ] ),
			   RealPoint( b[ 0 ], b[ 1 ] ),
			   RealPoint( c[ 0 ], c[ 1 ] ),
			   tvT.u( V[ 0 ] ),
			   tvT.u( V[ 1 ] ),
			   tvT.u( V[ 2 ] ) );
    }
    void viewTVTFlatTriangle( TVT & tvT, Face f )
    {
      VertexRange V = tvT.T.verticesAroundFace( f );
      Point       a = tvT.T.position( V[ 0 ] );
      Point       b = tvT.T.position( V[ 1 ] );
      Point       c = tvT.T.position( V[ 2 ] );
      Value     val = tvT.u( V[ 0 ] )
	+ tvT.u( V[ 1 ] ) + tvT.u( V[ 2 ] );
      val          /= 3.0;
      drawFlatTriangle( RealPoint( a[ 0 ], a[ 1 ] ),
			RealPoint( b[ 0 ], b[ 1 ] ),
			RealPoint( c[ 0 ], c[ 1 ] ), val );
    }

    void viewTVTTriangleDiscontinuity( TVT & tvT, Face f )
    {
      VertexRange V = tvT.T.verticesAroundFace( f );
      Point       a = tvT.T.position( V[ 0 ] );
      Point       b = tvT.T.position( V[ 1 ] );
      Point       c = tvT.T.position( V[ 2 ] );
      Value     val = { 255, 0, 0 };
      drawFlatTriangle( RealPoint( a[ 0 ], a[ 1 ] ),
			RealPoint( b[ 0 ], b[ 1 ] ),
			RealPoint( c[ 0 ], c[ 1 ] ), val );
    }


    /**
       Displays the AVT with flat or Gouraud shading.
    */
    void view( TVT & tvT )
    {
      cairo_set_operator( _cr,  CAIRO_OPERATOR_ADD );
      cairo_set_line_width( _cr, 0.0 ); 
      cairo_set_line_cap( _cr, CAIRO_LINE_CAP_BUTT );
      cairo_set_line_join( _cr, CAIRO_LINE_JOIN_BEVEL );
      for ( Face f = 0; f < tvT.T.nbFaces(); ++f )
	{
	  if ( _shading == 1 )      viewTVTGouraudTriangle( tvT, f );
	  else if ( _shading == 2 ) viewTVTLinearGradientTriangle( tvT, f );
	  else                      viewTVTFlatTriangle   ( tvT, f );
	}
    }

    /**
       Displays the AVT with flat or Gouraud shading, and displays a
       set of discontinuities as a percentage of the total energy.
    */
    void view( TVT & tvT, Scalar discontinuities )
    {
      // We need first to sort faces according to their energyTV.
      std::vector<Face> tv_faces( tvT.T.nbFaces() );
      for ( Face f = 0; f < tvT.T.nbFaces(); ++f )
	tv_faces[ f ] = f;
      std::sort( tv_faces.begin(), tv_faces.end(),
		 [ & tvT ] ( Face f1, Face f2 ) -> bool
		 // { return ( tvT.energyTV( f1 ) ) > ( tvT.energyTV( f2 ) ); } 
		 // { return ( tvT.aspectRatio( f1 ) )
                 //     > ( tvT.aspectRatio( f2 ) ); } 
		 { return ( tvT.energyTV( f1 ) * tvT.diameter( f1 ) )
		     > ( tvT.energyTV( f2 ) * tvT.diameter( f2 ) ); } 
		 // { return ( tvT.energyTV( f1 ) * tvT.aspectRatio( f1 ) )
		 //     > ( tvT.energyTV( f2 ) * tvT.aspectRatio( f2 ) ); } 
                 );
      Scalar Etv = tvT.getEnergyTV();
      Scalar Ctv = 0.0;
      Scalar Otv = Etv * discontinuities;
      
      cairo_set_operator( _cr,  CAIRO_OPERATOR_ADD );
      cairo_set_line_width( _cr, 0.0 ); 
      cairo_set_line_cap( _cr, CAIRO_LINE_CAP_BUTT );
      cairo_set_line_join( _cr, CAIRO_LINE_JOIN_BEVEL );
      for ( int i = 0; i < tv_faces.size(); ++i )
	{
	  Face f = tv_faces[ i ];
	  Ctv   += tvT.energyTV( f );
	  if ( Ctv < Otv ) { // display discontinuity
	    // viewTVTTriangleDiscontinuity( tvT, f );
	    viewTVTNonLinearGradientTriangle( tvT, f );
	  }
	  else {
	    if ( _shading == 1 )      viewTVTGouraudTriangle( tvT, f );
	    else if ( _shading == 2 ) viewTVTLinearGradientTriangle( tvT, f );
	    else                      viewTVTFlatTriangle   ( tvT, f );
	  }
	}
      // cairo_set_operator( _cr,  CAIRO_OPERATOR_OVER );
      // for ( int idx = 0; idx < tvT.T.nbVertices(); ++idx ) {
      // 	Point       a = tvT.T.position( idx );
      // 	Value     val = tvT.u( idx );
      // 	//cairo_set_source_rgb( _cr, 1.0, 0.0, 0.0 );
      // 	cairo_set_source_rgb( _cr, val[ 0 ] * _redf, val[ 1 ] * _greenf, val[ 2 ] * _bluef );
      // 	cairo_set_line_width( _cr, 0.0 );
      // 	cairo_rectangle( _cr, i( a[ 0 ] ), j( a[ 1 ] )-1, 1.0, 1.0 );
      // 	cairo_fill( _cr );
      // }
    }

  };
  
  // shading; 0:flat, 1:gouraud, 2:linear gradient.
  void viewTVTriangulation
  ( TVTriangulation& tvT, double b, double x0, double y0, double x1, double y1,
	int shading, bool color, std::string fname, double discontinuities,
    double stiffness, double amplitude )
  {
    CairoViewerTV cviewer
      ( x0, y0, x1, y1,
	b, b, shading, color, stiffness, amplitude );
    // CairoViewerTV cviewer
    //   ( (int) round( x0 ), (int) round( y0 ), 
    // 	(int) round( (x1+1 - x0) * b ), (int) round( (y1+1 - y0) * b ), 
    // 	b, b, shading, color, stiffness, amplitude );
    cviewer.view( tvT, discontinuities );
    cviewer.save( fname.c_str() );
  }
  
  void viewTVTriangulationAll
  ( TVTriangulation& tvT, double b, double x0, double y0, double x1, double y1,
    bool color, std::string fname, int display, double discontinuities,
    double stiffness, double amplitude )
  {
    if ( display & 0x1 )
      viewTVTriangulation( tvT, b, x0, y0, x1, y1, 0, color, fname + ".png",
			   discontinuities, stiffness, amplitude );
    if ( display & 0x2 )
      viewTVTriangulation( tvT, b, x0, y0, x1, y1, 1, color, fname + "-g.png",
			   discontinuities, stiffness, amplitude );
    if ( display & 0x4 )
      viewTVTriangulation( tvT, b, x0, y0, x1, y1, 2, color, fname + "-lg.png",
			   discontinuities, stiffness, amplitude );
  }

  void exportEPSMesh(TVTriangulation& tvT, const std::string &name, unsigned int width,
                     unsigned int height, bool displayMesh=true, double scale=1.0)
  {
    BasicVectoImageExporter exp( name, width, height, displayMesh, 100);
    BasicVectoImageExporter expMean( "mean.eps", width, height, displayMesh, scale);
    
    for(TVTriangulation::Face f = 0; f < tvT.T.nbFaces(); f++)
    {
      TVTriangulation::VertexRange V = tvT.T.verticesAroundFace( f );
      std::vector<TVTriangulation::Point> tr;
      tr.push_back(tvT.T.position(V[0]));
      tr.push_back(tvT.T.position(V[1]));
      tr.push_back(tvT.T.position(V[2]));
      tr.push_back(tvT.T.position(V[0]));
      TVTriangulation::Value   valMedian;
      TVTriangulation::Value   valMean = tvT.u( V[ 0 ] )+ tvT.u( V[ 1 ] ) + tvT.u( V[ 2 ] );
      valMean /= 3.0;
      double val1 = sqrt(tvT.u( V[ 0 ] )[0]*tvT.u( V[ 0 ] )[0]+
                         tvT.u( V[ 0 ] )[1]*tvT.u( V[ 0 ] )[1]+
                         tvT.u( V[ 0 ] )[2]*tvT.u( V[ 0 ] )[2]);
      double val2 = sqrt(tvT.u( V[ 1 ] )[0]*tvT.u( V[ 1 ] )[0]+
                         tvT.u( V[ 1 ] )[1]*tvT.u( V[ 1 ] )[1]+
                         tvT.u( V[ 1 ] )[2]*tvT.u( V[ 1 ] )[2]);
      double val3 = sqrt(tvT.u( V[ 2 ] )[0]*tvT.u( V[ 2 ] )[0]+
                         tvT.u( V[ 2 ] )[1]*tvT.u( V[ 2 ] )[1]+
                         tvT.u( V[ 2 ] )[2]*tvT.u( V[ 2 ] )[2]);
      if ((val1>=val2 && val1 <= val3) ||
         (val1<=val2 && val1 >= val3))
      {
        valMedian = tvT.u( V[ 0 ] );
      }
      else if ((val2<=val1 && val2 >= val3)||
               (val2>=val1 && val2 <= val3))
      {
        valMedian = tvT.u( V[ 1 ] );
      }else{
         valMedian = tvT.u( V[ 2 ] );
      }

      exp.addRegion(tr, DGtal::Color(valMedian[0], valMedian[1], valMedian[2]), 0.001);
      expMean.addRegion(tr, DGtal::Color(valMean[0], valMean[1], valMean[2]), 0.001);  
    }
    
    
  }

  TVTriangulation::Arc pivotNext(TVTriangulation& tvT, TVTriangulation::Arc a,
                 const TVTriangulation::Value &valTrack)
  {
    TVTriangulation::Value currentHead =  tvT.u(tvT.T.head(a));
    while( currentHead[0] == valTrack[0] && currentHead[1] == valTrack[1]  && currentHead[2] == valTrack[2]  )
      {
        a = tvT.T.next(a); 
        currentHead =  tvT.u(tvT.T.head(a));
      }
      return tvT.T.opposite(a);
  }



  std::vector<TVTriangulation::Point> trackBorderFromFace(TVTriangulation& tvT,  TVTriangulation::Face startArc,
                                                          TVTriangulation::Value valInside, std::vector<bool> &markedArcs)
  {
    std::vector<TVTriangulation::Point> res;
    
    // starting ext point: arc tail
    TVTriangulation::Face faceIni = tvT.T.faceAroundArc(startArc);

    TVTriangulation::Arc currentArc = startArc;
    TVTriangulation::Face currentFace = faceIni;
    markedArcs[startArc] = true;
    
    do 
    {
      TVTriangulation::VertexRange V = tvT.T.verticesAroundFace( currentFace );
      TVTriangulation::Point center = tvT.barycenter(currentFace);//(tvT.T.position(V[0])+tvT.T.position(V[1])+tvT.T.position(V[2]))/3.0;
      res.push_back(center);
      currentArc = pivotNext(tvT, currentArc, valInside);          
      currentFace = tvT.T.faceAroundArc(currentArc);
     if(currentFace == TVTriangulation::Triangulation::INVALID_FACE )
        {
            break;
        }
        markedArcs[currentArc] = true;
    } while(currentFace != faceIni);
    return res;
  }
  
  
  std::vector<std::vector<TVTriangulation::Point> > trackBorders(TVTriangulation& tvT, unsigned int num)
  {
    typedef std::map<DGtal::Color, std::vector<unsigned int> >  MapColorContours;
    MapColorContours mapContours;
    std::vector<std::vector<TVTriangulation::Point> > resAll;
    std::vector<bool> markedArcs(tvT.T.nbArcs());
    for(unsigned int i = 0; i< markedArcs.size(); i++){ markedArcs[i]=false; }
    bool found = true;
    while(found){
      found = false;
      for(unsigned int a = 0; a< markedArcs.size(); a++)
      {
          // tracking Head color
        TVTriangulation::Value valH = tvT.u(tvT.T.head(a));
        TVTriangulation::Value valT = tvT.u(tvT.T.tail(a));
        
        found = !markedArcs[a] && (valH[0]!=valT[0] || valH[1]!=valT[1] || valH[2]!=valT[2]) && !tvT.T.isArcBoundary(a);
        if(found)
        {
          resAll.push_back( trackBorderFromFace(tvT, a, valH, markedArcs));
          if (mapContours.count(DGtal::Color(valH[0], valH[1], valH[2]))==0)
          {
            std::vector<unsigned int> indexC;
            indexC.push_back(resAll.size()-1);
            mapContours[DGtal::Color(valH[0], valH[1], valH[2])]=indexC;
          }
          else
          {
            mapContours[DGtal::Color(valH[0], valH[1], valH[2])].push_back(resAll.size()-1);
          }
        }
      }
    }
    auto itMap = mapContours.begin();
      for(unsigned int i=0;i < num; i++)  itMap++;
    std::vector<std::vector<TVTriangulation::Point> > res;
    for(unsigned int i=0; i< (itMap->second).size(); i++){
      res.push_back(resAll[(itMap->second)[i]]);
    }
    return res;
  }



  /**
   * Invalidate contour border and compute the median image color.
   * 
   **/
  DGtal::Color invalidateImageBorder(TVTriangulation& tvT )
  {
    std::vector<TVTriangulation::Value> vectColors;
    for(TVTriangulation::Arc a=0; a < tvT.T.nbArcs(); a++)
    {
      if(tvT.T.isArcBoundary(a))
      {
        vectColors.push_back(tvT.u(tvT.T.head(a)));
      }
    }
    for(TVTriangulation::Arc a=0; a < tvT.T.nbArcs(); a++)
    {
      if(tvT.T.isArcBoundary(a))
      {
        tvT.invalidate(tvT.T.head(a));
        tvT.invalidate(tvT.T.tail(a));
      }
    }
    std::sort(vectColors.begin(), vectColors.end(), [](const TVTriangulation::Value &a,
                                                       const TVTriangulation::Value &b){
        return a[0]*a[0]+a[1]*a[1]+a[2]*a[2] > b[0]*b[0]+b[1]*b[1]+b[2]*b[2];});
    TVTriangulation::Value valMed = vectColors[vectColors.size()/2]; 
    return DGtal::Color(valMed[0], valMed[1], valMed[2]);
  }

  

  

    
  std::vector<TVTriangulation::ColorContours> trackAllBorders(TVTriangulation& tvT, unsigned int width, unsigned int height)
  {
    typedef std::map<DGtal::Color, std::vector<unsigned int> >  MapColorContours;
    std::vector<TVTriangulation::ColorContours> res;
    MapColorContours mapContours;
    std::vector<std::vector<TVTriangulation::Point> > resAll;
    std::vector<bool> markedArcs(tvT.T.nbArcs());
    for(unsigned int i = 0; i< markedArcs.size(); i++){ markedArcs[i]=false; }
    DGtal::Color med = invalidateImageBorder(tvT);

    // Adding background:
    TVTriangulation::ColorContours c;
    c.first = med;
    std::vector<std::vector<TVTriangulation::Point>> bg = {{TVTriangulation::Point(0, height),
                                                            TVTriangulation::Point(width, height),
                                                            TVTriangulation::Point(width, 0 ),
                                                            TVTriangulation::Point(0,0)}};
    c.second = bg;
    res.push_back(c);
   
    
    bool found = true;
    while(found){
      found = false;
      for(unsigned int a = 0; a< markedArcs.size(); a++)
      {
        // tracking Head color
        TVTriangulation::Value valH = tvT.u(tvT.T.head(a));
        TVTriangulation::Value valT = tvT.u(tvT.T.tail(a));
        
        found = !markedArcs[a] && (valH[0]!=valT[0] || valH[1]!=valT[1] || valH[2]!=valT[2]) && !tvT.T.isArcBoundary(a)
          && !tvT.isinvalid(tvT.T.head(a));
        if(found)
        {
          resAll.push_back( trackBorderFromFace(tvT, a, valH, markedArcs));
          if (mapContours.count(DGtal::Color(valH[0], valH[1], valH[2]))==0)
          {
            std::vector<unsigned int> indexC;
            indexC.push_back(resAll.size()-1);
            mapContours[DGtal::Color(valH[0], valH[1], valH[2])]=indexC;
          }
          else
          {
            mapContours[DGtal::Color(valH[0], valH[1], valH[2])].push_back(resAll.size()-1);
          }
        }
      }
    }

    for(auto cc : mapContours)
    {
      TVTriangulation::ColorContours c;
      c.first = cc.first;
      for(unsigned int i=0; i< (cc.second).size(); i++){
        c.second.push_back(resAll[(cc.second)[i]]);
      }
      res.push_back(c);
    }
    return res;
  }



  
  void exportEPSMeshDual(TVTriangulation& tvT, const std::string &name, unsigned int width,
                         unsigned int height, bool displayMesh, unsigned int numColor, double scale)
  {
    BasicVectoImageExporter exp( name, width, height, displayMesh, scale);    
    for(TVTriangulation::VertexIndex v = 0; v < tvT.T.nbVertices(); v++)
    {
      std::vector<TVTriangulation::RealPoint> tr;
      auto outArcs = tvT.T.outArcs( v );
      for ( auto rit = outArcs.rbegin(), ritEnd = outArcs.rend();
	    rit != ritEnd; ++rit ) {
	auto a = *rit;
	tr.push_back( tvT.contourPoint( a ) );
	if ( tvT.T.isArcBoundary( a ) ) {
	  trace.warning() << "Boundary arc" << std::endl;
	} else {
	  auto f = tvT.T.faceAroundArc( a );
	  tr.push_back( tvT.barycenter( f ) );
	}
      }
      // TVTriangulation::FaceRange F = tvT.T.facesAroundVertex( v );
      // for(auto f: F)
      // {
      //   TVTriangulation::Point center = tvT.barycenter( f );
	
      //   // TVTriangulation::VertexRange V = tvT.T.verticesAroundFace( f );
      //   // TVTriangulation::Point center = tvT.T.position(V[0])+tvT.T.position(V[1])+tvT.T.position(V[2]);
      //   // center /= 3.0;
      //   tr.push_back(center);
      // }
      TVTriangulation::Value val = tvT.u(v);
      exp.addRegion(tr, DGtal::Color(val[0], val[1], val[2]), 0.001);        
      
    }
    if(displayMesh)
      {
        for(TVTriangulation::Face f = 0; f < tvT.T.nbFaces(); f++)
	  {
	    TVTriangulation::VertexRange V = tvT.T.verticesAroundFace( f );
	    std::vector<TVTriangulation::Point> tr;
	    tr.push_back(tvT.T.position(V[0]));
	    tr.push_back(tvT.T.position(V[1]));
	    tr.push_back(tvT.T.position(V[2]));
	    tr.push_back(tvT.T.position(V[0]));
	    
	    exp.addContour(tr, DGtal::Color(0, 200, 200), 0.01);        
	  }
        std::vector<std::vector<TVTriangulation::Point> > contour = trackBorders(tvT, numColor);
        for (auto c: contour){
            DGtal::Color col;
            if(ContourHelper::isCounterClockWise(c) )
            {
                col = DGtal::Color(200, 20, 200);
            }
            else
            {
                col = DGtal::Color(200, 200, 20);
            }
            exp.addContour(c, col, 0.1);}
      } 
  }



  
  void exportEPSContoursDual(TVTriangulation& tvT, const std::string &name, unsigned int width,
                             unsigned int height, double scale)
  {
    BasicVectoImageExporter exp( name, width, height, false, scale);    
    std::vector<TVTriangulation::ColorContours> contourCol = trackAllBorders(tvT, width, height);
    for (auto c: contourCol){
      DGtal::Color col = c.first;
      exp.addRegions(c.second, col);
    }
    
  }



  


} // namespace DGtal




///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  using namespace DGtal;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "Specifies the input shape as a 2D image PPM filename.")
    ("limit,L", po::value<int>()->default_value(100), "Gives the maximum number of passes (default is 100).")  
    ("bitmap,b", po::value<double>()->default_value( 2.0 ), "Rasterization magnification factor [arg] for PNG export." )
    ("strategy,s", po::value<int>()->default_value(4), "Strategy for quadrilatera with equal energy: 0: do nothing, 1: subdivide, 2: flip all, 3: flip all when #flipped normal = 0, 4: flip approximately half when #flipped normal = 0, 5: flip approximately half when #flipped normal = 0, subdivide if everything fails." )
    ("tv-power,p", po::value<double>()->default_value( 0.5 ), "The power coefficient used to compute the gradient ie |Grad I|^{2p}. " )
    ("lambda,l", po::value<double>()->default_value( 0.0 ), "The data fidelity term in TV denoising (if lambda <= 0, then the data fidelity is exact" ) 
    ("dt", po::value<double>()->default_value( 0.248 ), "The time step in TV denoising (should be lower than 0.25)" ) 
    ("tolerance,t", po::value<double>()->default_value( 0.01 ), "The tolerance to stop the TV denoising." ) 
    ("quantify,q", po::value<int>()->default_value( 256 ), "The quantification for colors (number of levels, q=2 means binary." ) 
    ("tv-max-iter,N", po::value<int>()->default_value( 10 ), "The maximum number of iteration in TV's algorithm." )
    ("nb-alt-iter,A", po::value<int>()->default_value( 1 ), "The number of iteration for alternating TV and TV-flip." )
    ("display-tv,d", po::value<int>()->default_value( 0 ), "Tells the display mode after TV of output files per bit: 0x1 : output Flat colored triangles, 0x2 : output Gouraud colored triangles, 0x4: output Linear Gradient triangles." )
    ("display-flip,D", po::value<int>()->default_value( 4 ), "Tells the display mode after flips of output files per bit: 0x1 : output Flat colored triangles, 0x2 : output Gouraud colored triangles, 0x4: output Linear Gradient triangles." )
    ("discontinuities", po::value<double>()->default_value( 0.0 ), "Tells to display a % of the TV discontinuities (the triangles with greatest energy)." ) 
    ("stiffness", po::value<double>()->default_value( 0.9 ), "Tells how to stiff the gradient around discontinuities (amplitude value is changed at stiffness * middle)." ) 
    ("amplitude", po::value<double>()->default_value( 0.75 ), "Tells the amplitude of the stiffness for the gradient around discontinuities." )
    ("similarity", po::value<double>()->default_value( 0.0 ), "Tells when two colors are considered identical for connectedness." )
    ("displayMesh", "display mesh of the eps display." )
    ("exportEPSMesh,e", po::value<std::string>(), "Export the triangle mesh." )
    ("exportEPSMeshDual,E", po::value<std::string>(), "Export the triangle mesh." )
    ("exportEPSContoursDual,C", po::value<std::string>(), "Export the image regions filled." )
    ("epsScale", po::value<double>()->default_value( 1.0 ), "Change the default eps scale to increase display size on small images (using 10 will display easely small images while 1.0 is more adapted to bigger images) . " )
    ("numColorExportEPSDual", po::value<unsigned int>()->default_value(0), "num of the color of the map." )
    ("regularizeContour,R", po::value<int>()->default_value( 0 ), "regularizes the dual contours for <nb> iterations." )
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
  if( ! parseOK || vm.count("help") || argc <= 1 || (! vm.count( "input") ) )
    {
      trace.info()<< "Generate the best piecewise linear TV from a color image. It works by simply flipping edges of some initial triangulation to optimize the given energy." <<std::endl << "Basic usage: " << std::endl
		  << "\ttv-triangulation-color [options] -i <image.ppm> -b 4"<<std::endl
		  << general_opt << "\n";
      return 0;
    }

  // Useful types
  typedef DGtal::Z2i::Domain Domain;

  using namespace DGtal;

  trace.beginBlock("Construction of the triangulation");
  typedef ImageSelector < Z2i::Domain, unsigned int>::Type Image;
  typedef ImageSelector < Z2i::Domain, Color>::Type ColorImage;
  
  std::string img_fname = vm[ "input" ].as<std::string>();
  Image image           = GenericReader<Image>::import( img_fname ); 
  std::string extension = img_fname.substr(img_fname.find_last_of(".") + 1);
  bool            color = false;
  if ( extension == "ppm" ) color = true;
  trace.info() << "Image <" << img_fname
	       << "> size=" << image.extent()[ 0 ]
	       << "x" << image.extent()[ 1 ]
	       << " color=" << ( color ? "True" : "False" ) << std::endl;
  double    p = vm[ "tv-power" ].as<double>();
  double  sim = vm[ "similarity" ].as<double>();
  TVTriangulation TVT( image, color, p, sim );
  trace.info() << TVT.T << std::endl;
  trace.endBlock();

  trace.info() << std::fixed;
  trace.beginBlock("TV regularization");
  double lambda = vm[ "lambda" ].as<double>();
  double     dt = vm[ "dt" ].as<double>();
  double    tol = vm[ "tolerance" ].as<double>();
  int     quant = vm[ "quantify" ].as<int>();
  int         N = vm[ "tv-max-iter" ].as<int>();
  if ( lambda > 0.0 ) {
    TVT.tvPass( lambda, dt, tol, N );
  }
  if ( quant > 0 ) TVT.quantify( quant );
  trace.endBlock();
  
  trace.beginBlock("Output TV image (possibly quantified)");
  Image J( image.domain() );
  bool ok = TVT.outputU( J );
  struct UnsignedInt2Color {
    Color operator()( unsigned int val ) const { return Color( val ); }
  };
  PPMWriter<Image, UnsignedInt2Color>::exportPPM( "output-tv.ppm", J );
  trace.endBlock();

  trace.beginBlock("Displaying triangulation");
  {
    int  display = vm[ "display-tv" ].as<int>();
    double     b = vm[ "bitmap" ].as<double>();
    double  disc = vm[ "discontinuities" ].as<double>();
    double    st = vm[ "stiffness" ].as<double>();
    double    am = vm[ "amplitude" ].as<double>();
    double    x0 = 0.0;
    double    y0 = 0.0;
    double    x1 = (double) image.domain().upperBound()[ 0 ];
    double    y1 = (double) image.domain().upperBound()[ 1 ];
    viewTVTriangulationAll( TVT, b, x0, y0, x1, y1, color, "after-tv",
			    display, disc, st, am );
  }
  trace.endBlock();
  
  trace.beginBlock("Optimizing the triangulation");
  int miter = vm[ "limit" ].as<int>();
  int strat = vm[ "strategy" ].as<int>();
  int nbAlt = vm[ "nb-alt-iter" ].as<int>();
  if ( quant != 256 && nbAlt != 1 ) {
    nbAlt = 1;
    trace.warning() << "Quantification is not compatible with alternating TV + flips" << std::endl;
  }
  std::pair<int,int> nbs;
  for ( int n = 0; n < nbAlt; ++n ) {
    if ( n > 0 && lambda > 0.0 ) {
      TVT.tvPass( lambda, dt, tol, N );
    }
    int       iter = 0;
    int       last = 1;
    bool subdivide = false;
    trace.info() << "TV( u ) = " << TVT.getEnergyTV() << std::endl;
    while ( true ) {
      if ( iter++ > miter ) break;
      double energy = 0.0;
      nbs = TVT.onePass( energy, strat );
      if ( ( last == 0 ) && ( nbs.first == 0 ) ) {
	if ( subdivide || strat != 5 ) break;
	subdivide = true;
	nbs = TVT.onePass( energy, 1 );
      }
      last = nbs.first;
    }
  }
  trace.endBlock();

  trace.beginBlock("Displaying triangulation");
  {
    int  display = vm[ "display-flip" ].as<int>();
    double     b = vm[ "bitmap" ].as<double>();
    double  disc = vm[ "discontinuities" ].as<double>();
    double    st = vm[ "stiffness" ].as<double>();
    double    am = vm[ "amplitude" ].as<double>();
    double    x0 = 0.0;
    double    y0 = 0.0;
    double    x1 = (double) image.domain().upperBound()[ 0 ];
    double    y1 = (double) image.domain().upperBound()[ 1 ];
    viewTVTriangulationAll( TVT, b, x0, y0, x1, y1, color, "after-tv-opt",
			    display, disc, st, am );
  }
  trace.endBlock();

  trace.beginBlock("regularizing contours");
  {
    int         N = vm[ "regularizeContour" ].as<int>();
    TVT.regularizeContours( 0.001, N );
  }
  trace.endBlock();
  trace.beginBlock("Export base triangulation");

  double epsScale = vm["epsScale"].as<double>();
  if(vm.count("exportEPSMesh"))
    {
      unsigned int w = image.extent()[ 0 ];
      unsigned int h = image.extent()[ 1 ];
      std::string name = vm["exportEPSMesh"].as<std::string>();
      exportEPSMesh(TVT, name, w, h ,vm.count("displayMesh"), epsScale);
      
    }
  if(vm.count("exportEPSMeshDual"))
    {
      unsigned int w = image.extent()[ 0 ];
      unsigned int h = image.extent()[ 1 ];
      std::string name = vm["exportEPSMeshDual"].as<std::string>();
      unsigned int numColor = vm["numColorExportEPSDual"].as<unsigned int>();
      exportEPSMeshDual(TVT, name, w, h, vm.count("displayMesh"), numColor, epsScale);
      
    }
  if(vm.count("exportEPSContoursDual"))
  {
    unsigned int w = image.extent()[ 0 ];
    unsigned int h = image.extent()[ 1 ];
    std::string name = vm["exportEPSContoursDual"].as<std::string>();
    exportEPSContoursDual(TVT, name, w, h, epsScale);
  }
  trace.endBlock();

  
  // trace.beginBlock("Merging triangles");
  // {
    
  // }
  // trace.endBlock();
  return 0;
}
