#include "DGtal/helpers/StdDefs.h"
#include "DGtal/helpers/Shortcuts.h"

int main( int argc, char* argv[] )
{
  using namespace DGtal;
  typedef Shortcuts<Z3i::KSpace> SH3;

  int idx = argc > 2 ? atoi( argv[ 2 ] ) : 0;
  auto params = SH3::defaultParameters();
  params( "surfaceTraversal", "BreadthFirst" ) // specifies breadth-first traversal
    ( "colormap", "Custom" ); // specifies the colormap
  auto vol       = SH3::makeBinaryImage( argv[ 1 ], params );
  auto K         = SH3::getKSpace( vol );
  auto surface   = SH3::makeLightDigitalSurface( vol, K, params );
  auto surfels   = SH3::getSurfelRange( surface, params );
  surfels        = SH3::getSurfelRange( surface, surfels[ idx ], params );
  auto cmap      = SH3::getColorMap( 0, surfels.size() / 20, params );
  SH3::Colors colors( surfels.size() );
  for ( unsigned int i = 0; i < surfels.size(); ++i )
    colors[ i ] = SH3::Color( 200, 200, 200 );
  for ( unsigned int i = 0; i < surfels.size() / 20; ++i ) colors[ i ] = cmap( i );
  colors[ 0 ] = SH3::Color::Red;
  bool ok        = SH3::saveOBJ( surface, SH3::RealVectors(), colors, "vol-primal-bft.obj" );
  for ( unsigned int i = 1; i < 5; ++i )
    colors[ i ] = SH3::Color::Green;
  for ( unsigned int i = 5; i < surfels.size(); ++i )
    colors[ i ] = SH3::Color( 200, 200, 200 );
  bool ok2       = SH3::saveOBJ( surface, SH3::RealVectors(), colors, "vol-primal-neigh.obj" );
  return (ok && ok2) ? 0 : 1;
}

