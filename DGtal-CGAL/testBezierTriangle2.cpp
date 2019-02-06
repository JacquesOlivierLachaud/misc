#include <cfloat>
#include <iostream>
#include <vector>
#include <string>
#include <DGtal/base/Common.h>
#include <DGtal/base/ConstAlias.h>
#include <DGtal/base/Alias.h>
#include <DGtal/helpers/StdDefs.h>
#include "BezierTriangle2.h"
#include "BezierTriangle3.h"
#include "BezierCurve.h"
#include "CairoViewer.h"

int main( int argc, char** argv )
{
  using namespace DGtal;
  using namespace Z2i;

  // typedef BezierTriangle3< Space, 3 > BT;
  // typedef BT::Value                   Value;
  // typedef CairoViewer< Space >        Viewer;
  // RealPoint a( 10, 0 );
  // RealPoint b( 20, 0 );
  // RealPoint c( 10, 10 );
  // RealPoint d( 0, 5 );
  // Value va( 255, 255, 255 );
  // Value vb( 0, 0, 0 );
  // Value vc( 0, 0, 0 );
  // Value vd( 100, 100, 100 );
  // Value vabc = BT::midValue( { va, vb, vc } );
  // Value vacd = BT::midValue( { va, vc, vd } );
  // RealVector o_ab = 0.5 * BT::grad( { a, b, c }, { va[ 0 ], vb[ 0 ], vc[ 0 ] } )
  //   + 0.5 * BT::grad( { b, a, RealPoint( 15, -5 ) }, { vb[ 0 ], va[ 0 ], 85 } );
  // RealVector o_bc = 0.5 * BT::grad( { a, b, c }, { va[ 0 ], vb[ 0 ], vc[ 0 ] } )
  //   + 0.5 * BT::grad( { c, b, RealPoint( 18, 8 ) }, { vc[ 0 ], vb[ 0 ], 0 } );
  // RealVector o_ca = 0.5 * BT::grad( { a, b, c }, { va[ 0 ], vb[ 0 ], vc[ 0 ] } )
  //   + 0.5 * BT::grad( { a, c, (a+c+d/3.0) }, { va[ 0 ], vc[ 0 ], vacd[ 0 ] } );
  // std::cout << "o_ab=" << o_ab << std::endl;
  // std::cout << "o_bc=" << o_bc << std::endl;
  // std::cout << "o_ca=" << o_ca << std::endl;
  // BT bt1( { a, b, c }, { va, vb, vc },
  // 	  { o_ab, o_bc, o_ca },
  // 	  vabc,
  // 	  { RealPoint( 15, -5 ), RealPoint( 18, 8 ), (a+c+d)/3.0 },
  // 	  { Value( 85, 85, 85 ), Value( 0, 0, 0 ), vacd } );
  // // BT bt2( { a, c, d }, { va, vc, vd },
  // // 	  { -o_ca, o_cd, o_da },
  // // 	  vacd,
  // // 	  { (a+b+c)/3.0, RealPoint( 5, 10 ), RealPoint( 4, 0 ) }, 
  // // 	  { vabc, Value( 20, 20, 20 ), Value( 40, 40, 40) } );

  // Viewer viewer( 0.0, 0.0, 20.0, 10.0, 4.0, 4.0 );
  // //viewer.drawLinearGradientTriangle( a, b, c, va, vb, vc );
  // viewer.drawColorBezierTriangle( bt1 );
  // // viewer.drawColorBezierTriangle( bt2 );
  // viewer.save( "bt.png" );

  std::vector<RealPoint> P { {0,0}, {4,5}, {8,6}, {15, 3} };
  std::vector<Point>  Q;
  std::vector<double> T;

  BezierCurve<Space> B( P );
  std::vector<RealPoint> rp;
  B.trace( [] (RealPoint p) { return Point( round( p[ 0 ] ), round( p[ 1 ] ) ); },
	   Q, rp, T );
  for ( int i = 0; i < Q.size(); i++ )
    std::cout << i << " Q[i]=(" << Q[i][0] << "," << Q[i][1] << ")"
	      << " t[i]=" << T[ i ] << " B(t[i])="
	      << B( T[ i ] ) << std::endl;
  return 0;
}
