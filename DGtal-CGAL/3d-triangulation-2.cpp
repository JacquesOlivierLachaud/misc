#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <QtGui/qapplication.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>


#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/Shapes.h>
#include <DGtal/shapes/ShapeFactory.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>

#include "Auxiliary.h"
#include "Triangulation3DHelper.h"
#include "SimplicialStrip3D.h"
#include "RelativeConvexHull.h"

typedef DGtal::SpaceND<3, DGtal::int64_t> Z3;
typedef Z3::Point PointZ3;
typedef DGtal::HyperRectDomain<Z3> Domain;
typedef Domain::ConstIterator DomainConstIterator;

template <typename TKernel3>
struct toCGALFunctor {
  typedef typename TKernel3::Point_3 Point3;
  inline Point3 operator()( const PointZ3 & p ) const
  {
    return Point3( p[0], p[1], p[2] );
  }
};

template <typename TKernel3>
struct toDGtalFunctor {
  typedef typename TKernel3::Point_3 Point3;
  inline PointZ3 operator()( const Point3 & p ) const
  {
    return PointZ3( p.x(), p.y(), p.z() );
  }
};

inline
PointZ3 operator^( const PointZ3 & a, const PointZ3 & b )
{
  BOOST_STATIC_ASSERT( PointZ3::dimension == 3 );
  return PointZ3( a[ 1 ] * b[ 2 ] - a[ 2 ] * b[ 1 ],
		  a[ 2 ] * b[ 0 ] - a[ 0 ] * b[ 2 ],
		  a[ 0 ] * b[ 1 ] - a[ 1 ] * b[ 0 ] );
}



template <typename Viewer, typename ToDGtal, typename Facet>
void displayFacet( Viewer & viewer, const ToDGtal & toDGtal,
		   const Facet & f, const DGtal::Color & col )
{
  typedef typename ToDGtal::Point3 Point;
  Point a( f.first->vertex( (f.second+1)%4 )->point() );
  Point b( f.first->vertex( (f.second+2)%4 )->point() );
  Point c( f.first->vertex( (f.second+3)%4 )->point() );
  viewer.addTriangle( a.x(), a.y(), a.z(),
		      c.x(), c.y(), c.z(),
		      b.x(), b.y(), b.z(),
		      col );
}

template <typename Viewer, typename ToDGtal, typename Cell>
void displayCell( Viewer & viewer, const ToDGtal & toDGtal,
		  const Cell & c, const DGtal::Color & col )
{
  for ( int i = 0; i < 4; ++i )
    displayFacet( viewer, toDGtal, std::make_pair( c, i ), col );
}



// A modifier creating a triangle with the incremental builder.
template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> {
public:
    Build_triangle() {}
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( 3, 1, 6);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        B.add_vertex( Point( 0, 0, 0));
        B.add_vertex( Point( 1, 0, 0));
        B.add_vertex( Point( 0, 1, 0));
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.end_surface();
    }
};


DGtal::int64_t
countLatticePointsInTetrahedra( const PointZ3 & a, const PointZ3 & b, const PointZ3 & c, const PointZ3 & d )
{
  DGtal::int64_t nb = 0;
  PointZ3 ab = b - a;
  PointZ3 bc = c - b;
  PointZ3 cd = d - c;  
  PointZ3 da = a - d;
  PointZ3 abc = ab^bc;
  PointZ3 bcd = bc^cd;
  PointZ3 cda = cd^da;
  PointZ3 dab = da^ab;
  PointZ3::Component abc_shift = abc.dot( a );
  PointZ3::Component bcd_shift = bcd.dot( b );
  PointZ3::Component cda_shift = cda.dot( c );
  PointZ3::Component dab_shift = dab.dot( d );
  if ( abc.dot( d ) < abc_shift )  { abc = -abc; abc_shift = -abc_shift; }
  if ( bcd.dot( a ) < bcd_shift )  { bcd = -bcd; bcd_shift = -bcd_shift; }
  if ( cda.dot( b ) < cda_shift )  { cda = -cda; cda_shift = -cda_shift; }
  if ( dab.dot( c ) < dab_shift )  { dab = -dab; dab_shift = -dab_shift; }
  PointZ3 inf = a.inf( b ).inf( c ).inf( d );
  PointZ3 sup = a.sup( b ).sup( c ).sup( d );
  Domain domain( inf, sup );
  for ( DomainConstIterator it = domain.begin(), itE = domain.end();
	it != itE; ++it )
    {
      PointZ3 p = *it;
      if ( ( abc.dot( p ) >= abc_shift )
	   && ( bcd.dot( p ) >= bcd_shift )
	   && ( cda.dot( p ) >= cda_shift )
	   && ( dab.dot( p ) >= dab_shift ) )
	++nb;
    }
  return nb;
}

/**
   @return 4 iff there are 4 lattice points in tetrahedra, otherwise returns 5.
*/
DGtal::int64_t
countLatticePointsInTetrahedraIf4( const PointZ3 & a, const PointZ3 & b, const PointZ3 & c, const PointZ3 & d )
{
  DGtal::int64_t nb = 0;
  PointZ3 ab = b - a;
  PointZ3 bc = c - b;
  PointZ3 cd = d - c;  
  PointZ3 da = a - d;
  PointZ3 abc = ab^bc;
  PointZ3 bcd = bc^cd;
  PointZ3 cda = cd^da;
  PointZ3 dab = da^ab;
  PointZ3::Component abc_shift = abc.dot( a );
  PointZ3::Component bcd_shift = bcd.dot( b );
  PointZ3::Component cda_shift = cda.dot( c );
  PointZ3::Component dab_shift = dab.dot( d );
  if ( abc.dot( d ) < abc_shift )  { abc = -abc; abc_shift = -abc_shift; }
  if ( bcd.dot( a ) < bcd_shift )  { bcd = -bcd; bcd_shift = -bcd_shift; }
  if ( cda.dot( b ) < cda_shift )  { cda = -cda; cda_shift = -cda_shift; }
  if ( dab.dot( c ) < dab_shift )  { dab = -dab; dab_shift = -dab_shift; }
  PointZ3 inf = a.inf( b ).inf( c ).inf( d );
  PointZ3 sup = a.sup( b ).sup( c ).sup( d );
  Domain domain( inf, sup );
  for ( DomainConstIterator it = domain.begin(), itE = domain.end();
	it != itE; ++it )
    {
      PointZ3 p = *it;
      if ( ( abc.dot( p ) >= abc_shift )
	   && ( bcd.dot( p ) >= bcd_shift )
	   && ( cda.dot( p ) >= cda_shift )
	   && ( dab.dot( p ) >= dab_shift ) )
	{ 
          ++nb;
          if ( nb > 4 ) return 5;
        }
    }
  return nb;
}


namespace DGtal {
  template <typename Kernel>
  class DigitalCore {
  public:
    typedef typename CGAL::Delaunay_triangulation_3<Kernel> Delaunay;
    typedef typename Delaunay::Cell_iterator  CellIterator;
    typedef typename Delaunay::Cell_handle    CellHandle;
    typedef typename Delaunay::Facet          Facet;
    typedef typename Delaunay::Facet_iterator FacetIterator;
    typedef typename Delaunay::Vertex_handle  VertexHandle;
    typedef DGtal::uint64_t               Size;
    typedef typename std::map< CellHandle, Size >  MapCell2Size;
    typedef typename std::set< Facet >             FacetSet;
    typedef typename std::set< CellHandle >        CellSet;
    
  private:
    /// The Delaunay triangulation that carries the digital core.
    const Delaunay* myDelaunay;
    /// Stores the number of lattice of points within each triangulation Cell.
    MapCell2Size myNbLatticePoints;
    /// Current boundary as a set of facets.
    FacetSet myBoundary;
    /// Current interior cells.
    CellSet myInterior;

    toDGtalFunctor<Kernel> toDGtal;
   
    //--------------------------------------------------------------------------
  public:

    inline 
    DigitalCore() : myDelaunay( 0 ) 
    {}
    
    inline 
    DigitalCore( ConstAlias<Delaunay> t ) 
      : myDelaunay( t )
    {
      countLatticePoints();
      computeBasicFacets();
    }

    //--------------------------------------------------------------------------
  public:

    inline
    const FacetSet & boundary() const
    { return myBoundary; }

    inline
    FacetSet & boundary()
    { return myBoundary; }

    inline
    const CellSet & interior() const
    { return myInterior; }

    inline
    bool isFacetBasic( const Facet & f ) const
    {
      if ( ! myDelaunay->is_infinite( f ) )
        {
          int i = f.second;
          PointZ3 a( toDGtal( f.first->vertex( (i+1)%4 )->point() ) );
          PointZ3 b( toDGtal( f.first->vertex( (i+2)%4 )->point() ) );
          PointZ3 c( toDGtal( f.first->vertex( (i+3)%4 )->point() ) );
          return ( ( (b-a).norm( PointZ3::L_infty ) == 1 )
                   && ( (c-b).norm( PointZ3::L_infty ) == 1 )
                   && ( (a-c).norm( PointZ3::L_infty ) == 1 ) );
        }
      else
        return false;
    }

    inline 
    Size nbLatticePoints( const CellHandle & cell ) const
    {
      typename MapCell2Size::const_iterator it = myNbLatticePoints.find( cell );
      ASSERT( it != myNbLatticePoints.end() );
      return it->second;
    }

    inline 
    bool isCellInterior( const CellHandle & c ) const
    {
      return myInterior.find( c ) != myInterior.end();
    }

    inline 
    bool isFacetBoundary( const Facet & f ) const
    {
      return myBoundary.find( f ) != myBoundary.end();
    }

    inline 
    bool isFacetInterior( const Facet & f ) const
    {
      return myInterior.find( f.first ) != myInterior.end();
    }

    inline 
    bool isFacetExterior( const Facet & f ) const
    {
      return ( myInterior.find( f.first ) == myInterior.end() )
        && ! isFacetBoundary( f );
    }

    inline
    bool isVertexInterior( const VertexHandle & vh ) const
    {
      std::vector<CellHandle> incident_cells;
      myDelaunay->incident_cells( vh, std::back_inserter( incident_cells ) );
      for ( typename std::vector<CellHandle>::const_iterator it = incident_cells.begin(),
              itend = incident_cells.end(); it != itend; ++it )
        {
          if ( ! isCellInterior( *it ) ) return false;
        }
      return true;
    }

    inline
    double area( const Facet & f ) const
    {
      return sqrt( myDelaunay->triangle( f ).squared_area() );
    }

    void extend()
    {
      ASSERT( myDelaunay != 0 );
      trace.beginBlock("[DigitalCore] Extend boundary.");

      // Queue for computing the digital core.
      CellSet priorityQ;
      
      // Prepare queue.
      for ( typename FacetSet::const_iterator it = myBoundary.begin(), 
	      itend = myBoundary.end();
            it != itend; ++it )
        priorityQ.insert( it->first );
      
      // Start extension
      bool isMarked[ 4 ];
      unsigned int nb = 0;
      while ( ! priorityQ.empty() )
        {
          typename CellSet::iterator it = priorityQ.begin();
          CellHandle cell_h = *it;
          priorityQ.erase( it );
          // infinite cells are not processed.
          if ( myDelaunay->is_infinite( cell_h ) ) continue;
          // cells containing integer points in the interior are not processed.
          if ( nbLatticePoints( cell_h ) != 4 )         continue;
          // checking cell faces for extension.
	  extendCell( priorityQ, isMarked, cell_h );
	  ++nb;
        }
      trace.info() << "- cells flipped = " << nb << std::endl;
      trace.info() << "- boundary has " << myBoundary.size() << " facets." << std::endl;
      trace.endBlock();
    }

        
    //--------------------------------------------------------------------------
  protected:

    // bool checkCellExtension( bool isBoundary[ 4 ], const CellHandle & cell ) const
    // {
    //   unsigned int n = 0;
    //   for ( int i = 0; i < 4; ++i )
    //     {
    //       isBoundary[ i ] = isFacetBoundary( Facet( cell, i ) );
    //       if ( isBoundary[ i ] ) ++n;
    //     }
    //   // At least 2 are marked, by convexity, we can close the gap
    //   // to the further faces of the cell.
    //   bool propagate = n >= 3;
    //   if ( n == 2 )
    //     { // 2 are marked. We check their area.
    //       double area_boundary = 0.0;
    //       double area_exterior = 0.0;
    //       for ( int i = 0; i < 4; ++i ) {
    //         if ( isBoundary[ i ] ) 
    //           area_boundary += area( Facet( cell, i ) );
    //         else
    //           area_exterior += area( Facet( cell, i ) );
    //       }
    //       if ( area_exterior < myAreaFactorExtend * area_boundary )
    //         propagate = true;
    //     }
    //   return propagate;
    // }
    
    void extendCell( CellSet & priorityQ, bool isBoundary[ 4 ], const CellHandle & cell )
    { // Switch cell to interior
      myInterior.insert( cell );
      for ( unsigned int i = 0; i < 4; ++i ) {
        if ( ! isBoundary[ i ] )
          {
            Facet nfacet = myDelaunay->mirror_facet( Facet( cell, i ) );
            priorityQ.insert( nfacet.first );
            myBoundary.insert( nfacet );
          }
        else
          {
            myBoundary.erase( Facet( cell, i ) );
          }
      }
    }


    void computeBasicFacets()
    {
      ASSERT( myDelaunay != 0 );
      trace.beginBlock("[DigitalCore] Compute basic facets.");
      DGtal::uint64_t nb = 0;
      for( FacetIterator it = myDelaunay->facets_begin(), itend = myDelaunay->facets_end();
           it != itend; ++it )
        {
          if ( isFacetBasic( *it ) ) // || checkPlanarVFacet( f ) )
            {
	      myBoundary.insert( *it );
	      myBoundary.insert( myDelaunay->mirror_facet( *it ) );
              ++nb;
            }
        }
      trace.info() << "- nb basic facets = " << nb << std::endl;
      trace.endBlock();
    }

    void countLatticePoints()
    {
      ASSERT( myDelaunay != 0 );
      trace.beginBlock("[DigitalCore] Counting lattice points.");
      double maxbar = myDelaunay->number_of_cells();
      double bar = 0.0;
      int i = 0;
      for( CellIterator it = myDelaunay->cells_begin(), 
             itend = myDelaunay->cells_end();
           it != itend; ++it, ++bar, ++i )
        {
          if ( i % 1000 == 0 ) trace.progressBar( bar, maxbar );
          if ( myDelaunay->is_infinite( it ) ) 
            myNbLatticePoints[ it ] = -1;
          else
            {
              PointZ3 a( toDGtal(it->vertex(0)->point())),
                b(toDGtal(it->vertex(1)->point())),
                c(toDGtal(it->vertex(2)->point())),
                d(toDGtal(it->vertex(3)->point()));
              myNbLatticePoints[ it ] = countLatticePointsInTetrahedraIf4( a, b, c, d );
            }
        }
      trace.endBlock();
    }
    
  };

}





///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

/**
   Main function.

   @param argc the number of parameters given on the line command.

   @param argv an array of C-string, such that argv[0] is the name of
   the program, argv[1] the first parameter, etc.

  "29*z+47*y+23*x-5"
  "10*z-x^2+y^2-100"
  "z*x*y+x^4-5*x^2+2*y^2*z-z^2-1000"
  "(15*z-x^2+y^2-100)*(x^2+y^2+z^2-1000)" nice
  "(x^2+y^2+(z+5)^2-100)^2+10*(x^2+y^2+(z-3)^2)-2000" bowl
  "x^2+y^2+2*z^2-x*y*z+z^3-100" dragonfly
  "(x^2+y^2+(z+5)^2)^2-x^2*(z^2-x^2-y^2)-100" joli coeur
  "0.5*(z^2-4*4)^2+(x^2-7*7)^2+(y^2-7*7)^2-7.2*7.2*7.2*7.2" convexites et concavites
  "(1-(x^2+y^2))^2*(1+(x^2+y^2))^2-z" "cup"
*/
int main( int argc, char ** argv ) {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  // typedef CGAL::Simple_cartesian<double>     Kernel;
  typedef CGAL::Polyhedron_3<K>               Polyhedron;
  typedef CGAL::Triangulation_3<K>            Triangulation;
  typedef CGAL::Delaunay_triangulation_3<K>   Delaunay;
  typedef Polyhedron::HalfedgeDS              HalfedgeDS;
  typedef Delaunay::Finite_cells_iterator     FiniteCellsIterator;
  
  using namespace DGtal;

  typedef KhalimskySpaceND<3,DGtal::int64_t> K3;
  typedef Z3::Vector Vector;
  typedef Z3::RealPoint RealPoint;
  typedef K3::SCell SCell;
  typedef K3::SCellSet SCellSet;
  typedef SCellSet::const_iterator SCellSetConstIterator;
  // typedef ImplicitRoundedHyperCube<Z3> Shape;
  typedef RealPoint::Coordinate Ring;
  typedef toCGALFunctor<K>                    ToCGAL;
  typedef toDGtalFunctor<K>                   ToDGtal;
  typedef DigitalCore<K>                      DigCore;
  typedef DigCore::FacetSet                   FacetSet;
  typedef DigCore::Facet                      Facet;

  QApplication application(argc,argv); // remove Qt arguments.

  po::options_description general_opt("Specific allowed options (for Qt options, see Qt official site) are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("vol,v", po::value<std::string>(), "specifies the shape as some subset of a .vol file [arg]" )
    ("min,m", po::value<int>()->default_value( 1 ), "the minimum threshold in the .vol file: voxel x in shape iff m < I(x) <= M" )
    ("max,M", po::value<int>()->default_value( 255 ), "the maximum threshold in the .vol file: voxel x in shape iff m < I(x) <= M" )
    ("poly,p", po::value<std::string>(), "specifies the shape as the zero-level of the multivariate polynomial [arg]" )
    ("gridstep,g", po::value<double>()->default_value( 1.0 ), "specifies the digitization grid step ([arg] is a double) when the shape is given as a polynomial." )
    ("bounds,b", po::value<int>()->default_value( 20 ), "specifies the diagonal integral bounds [-arg,-arg,-arg]x[arg,arg,arg] for the digitization of the polynomial surface." )
    ("point-list,l", po::value<std::string>(), "specifies the shape as a list of digital points given in the filename [arg], x y z per line." )
    ("prune,P", po::value<double>(), "prunes the resulting the complex of the tetrahedra were the length of the dubious edge is [arg] times the length of the opposite edge." )
    ("geometry,G", po::value<int>()->default_value( 0 ), "specifies how digital points are defined from the digital surface: arg=0: inner voxels, arg=1: inner and outer voxels, arg=2: pointels" )
    ("area-factor-extend,A", po::value<double>()->default_value( sqrt(3.0) ), "specifies the convexity/concavity area ratio when extending the core." )
    ("area-factor-retract,a", po::value<double>()->default_value( 1.0/sqrt(3.0) ), "specifies the convexity/concavity area ratio when retracting the core." )
    ("view,V", po::value<int>()->default_value( 1 ), "specifies what is displayed. The bits of [arg] indicates what is displayed:  0: digital core, 1: empty unreachable tetrahedra, 2: retracted tetrahedra, 3: digital points" )
    ; 
  
  // parse command line ----------------------------------------------
  bool parseOK=true;
  po::variables_map vm;
  try {
    po::command_line_parser clp( argc, argv );
    clp.options( general_opt );
    po::store( clp.run(), vm );
  } catch( const std::exception& ex ) {
    parseOK = false;
    trace.info() << "Error checking program options: "<< ex.what() << std::endl;
  }
  po::notify( vm );    
  if( !parseOK || vm.count("help")||argc<=1
      || !( vm.count("poly") || vm.count("vol") || vm.count("point-list") ) )
    {
      std::cout << "Usage: " << argv[0] << " [options] {--vol <vol-file> || --poly <polynomial-string>}\n"
		<< "Computes the linear reconstruction of the given digital surface, specified either as a thresholded .vol file, or a zero-level of a multivariate polynomial, or a list of digital points.\n"
		<< general_opt << "\n\n";
      std::cout << "Example:\n"
		<< argv[0] << " -p \"x^2+y^2+2*z^2-x*y*z+z^3-100\" -g " << (double) 0.5 << std::endl;
      return 0;
    }
  
  
  // Construction of the digital surface ----------------------------------------------
  trace.beginBlock("Construction of the digital surface");
  K3 ks;
  SCellSet boundary;
  //  SurfelAdjacency<3> sAdj( true );
  int ok = 4;
  if ( vm.count( "vol" ) )
    { // .vol file
      ok = makeSpaceAndBoundaryFromVolFile( ks, boundary,
					    vm["vol"].as<std::string>(),
					    vm["min"].as<int>(),
					    vm["max"].as<int>() );
     }
  else if ( vm.count( "poly" ) )
    {
      ok = makeSpaceAndBoundaryFromPolynomialString( ks, boundary,
						     vm[ "poly" ].as<std::string>(),
						     vm[ "gridstep" ].as<double>(),
						     vm[ "bounds" ].as<int>() );
    }
  else if ( vm.count( "point-list" ) )
    {
      ok = makeSpaceAndBoundaryFromPointList( ks, boundary,
                                              vm[ "point-list" ].as<std::string>() );
    }
  trace.endBlock();
  if ( ok != 0 ) return ok; // error.


  trace.beginBlock("Constructing the set of points");
  std::vector<PointZ3> digital_points;
  int geometry = vm["geometry"].as<int>();
  if ( geometry == 0 )
    getInnerVoxelCoordinates( digital_points, ks, boundary.begin(), boundary.end() );
  else if ( geometry == 1 )
    getInnerAndOuterVoxelCoordinates( digital_points, ks, boundary.begin(), boundary.end() );
  else if ( geometry == 2 )
    getPointelCoordinates( digital_points, ks, boundary.begin(), boundary.end() );

  trace.beginBlock("Shuffle points.");
  random_shuffle( digital_points.begin(), digital_points.end() );
  trace.endBlock();
  trace.endBlock();

  trace.beginBlock("Creating the Delaunay complex.");
  ToCGAL toCGAL;
  ToDGtal toDGtal;
  Delaunay T;
  double setsize = (double) digital_points.size()-1;
  trace.info() << "Vertices to process: " << setsize << std::endl;
  double step = 0.0;
  int i = 0;
  for ( std::vector<PointZ3>::const_iterator it = digital_points.begin(), 
	  itE = digital_points.end();
	it != itE; ++it, ++step, ++i )
    {
      if ( i % 1000 == 0 ) trace.progressBar( step, setsize );
      T.insert( toCGAL( *it ) );
    }
  trace.endBlock();

  DigCore core( T );
  for ( unsigned int i = 0; i < 10; ++i )
    core.extend();


  // start viewer
  int view = vm["view"].as<int>();
  Viewer3D viewerCore;
  viewerCore.show();
  Color colBasicFacet2( 0, 255, 255, 255 );
  Color colBasicFacet1( 0, 255, 0, 255 );
  Color colSpuriousTetrahedra( 255, 0, 0, 100 );
  if ( view & 0x1 ) { // View digital core.
    for ( FacetSet::const_iterator it = core.boundary().begin(), itend = core.boundary().end();
          it != itend; ++it )
      {
        // we display it.
        // Triangle triangle = t.triangle( *it );
	PointZ3 a( toDGtal( it->first->vertex( (it->second+1)%4 )->point() ) );
	PointZ3 b( toDGtal( it->first->vertex( (it->second+2)%4 )->point() ) );
	PointZ3 c( toDGtal( it->first->vertex( (it->second+3)%4 )->point() ) );
        Facet f2 = T.mirror_facet( *it );
        if ( core.isFacetBoundary( f2 ) )
          { // the mirror facet is also in the triangulation. 
            // We need to move vertices a little bit when two triangles are at the same position.
            PointZ3 n = (b-a)^(c-a);
            double norm = n.norm(PointZ3::L_2);
            double dx[ 3 ];
            for ( unsigned int j = 0; j < 3; ++j )
              dx[ j ] = 0.001*((double) n[j])/norm;
            viewerCore.addTriangle( (double) a[ 0 ] + dx[ 0 ], (double) a[ 1 ] +  dx[ 1 ], (double) a[ 2 ] + dx[ 2 ],
                                    (double) c[ 0 ] + dx[ 0 ], (double) c[ 1 ] +  dx[ 1 ], (double) c[ 2 ] + dx[ 2 ],
                                    (double) b[ 0 ] + dx[ 0 ], (double) b[ 1 ] +  dx[ 1 ], (double) b[ 2 ] + dx[ 2 ],
                                    colBasicFacet2 );
          }
        else
          viewerCore.addTriangle( a[ 0 ], a[ 1 ], a[ 2 ],
                                  c[ 0 ], c[ 1 ], c[ 2 ],
                                  b[ 0 ], b[ 1 ], b[ 2 ],
                                  colBasicFacet1 );
      }
  } //  if ( view & 0x1 ) {

  // start viewer
  // int view = vm["view"].as<int>();
  // Viewer3D viewerRCH;
  // viewerRCH.show();
  // Color colBasicFacet1( 0, 255, 0, 255 );
  // if ( view & 0x1 ) { // View digital core.
  //   for ( FiniteCellsIterator it = T.finite_cells_begin(), 
  //           itend = T.finite_cells_end(); it != itend; ++it )
  //     {
  //       // we display it.
  // 	PointZ3 a( toDGtal( it->vertex( 0 )->point() ) );
  // 	PointZ3 b( toDGtal( it->vertex( 1 )->point() ) );
  // 	PointZ3 c( toDGtal( it->vertex( 2 )->point() ) );
  // 	PointZ3 d( toDGtal( it->vertex( 3 )->point() ) );
  // 	if ( countLatticePointsInTetrahedraIf4( a, b, c, d ) == 4 )
  // 	  displayCell( viewerRCH, toDGtal, it, colBasicFacet1 );
  //     }
  // }
  viewerCore << Viewer3D::updateDisplay;
  application.exec();

  std::cout << "number of vertices :  " ;
  std::cout << T.number_of_vertices() << std::endl;
  std::cout << "number of edges :  " ;
  std::cout << T.number_of_edges() << std::endl;
  std::cout << "number of facets :  " ;
  std::cout << T.number_of_facets() << std::endl;
  std::cout << "number of cells :  " ;
  std::cout << T.number_of_cells() << std::endl;

  
  
  Polyhedron P;
  Build_triangle<HalfedgeDS> triangle;
  P.delegate( triangle);
  CGAL_assertion( P.is_triangle( P.halfedges_begin()));
  return 0;
}
