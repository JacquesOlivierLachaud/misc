#include <cfloat>
#include <iostream>
#include <vector>
#include <string>
#include <DGtal/base/Common.h>
#include <DGtal/base/ConstAlias.h>
#include <DGtal/base/Alias.h>
#include <DGtal/helpers/StdDefs.h>
#include "BezierTriangle2.h"
#include "CairoViewer.h"

int main( int argc, char** argv )
{
  using namespace DGtal;
  using namespace Z2i;

  typedef BezierTriangle2< Space, 3 > BT;
  typedef BT::Value                   Value;
  typedef CairoViewer< Space >        Viewer;
  RealPoint a( 0, 0 );
  RealPoint b( 10, 0 );
  RealPoint c( 0, 10 );
  RealVector o_ab( 0, 1 );
  RealVector o_bc( 1, 1 );
  RealVector o_ca( 1, 1 );
  Value va( 255, 255, 255 );
  Value vb( 0, 255, 0 );
  Value vc( 0, 0, 255 );
  BT bt( { a, b, c }, { o_ab, o_bc, o_ca }, { va, vb, vc } );

  Viewer viewer( 0.0, 0.0, 15.0, 15.0, 16.0, 16.0 );
  //viewer.drawLinearGradientTriangle( a, b, c, va, vb, vc );
  viewer.drawBezierTriangle( bt );
  viewer.save( "bt.png" );

  return 0;
}
