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

int main ()
{
  
  Delaunay t;
  
  trace.beginBlock("Construction the shape");
  typedef Ellipse2D<Z2i::Space> Ellipse; 
  Ellipse2D<Z2i::Space> ellipse(Z2i::Point(0,0), 5, 3, 0.3 );
  double h = 0.02; 
  GaussDigitizer<Z2i::Space,Ellipse> dig;  
  dig.attach( ellipse );
  dig.init( ellipse.getLowerBound()+Z2i::Vector(-1,-1),
            ellipse.getUpperBound()+Z2i::Vector(1,1), h ); 
  // typedef Flower2D<Z2i::Space> Flower; 
  // Flower2D<Z2i::Space> flower(Z2i::Point(0,0), 15, 2, 5, 0);
  // double h = 0.75; 
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
    t.insert( Point( (*it)[0], (*it)[1]));
  trace.endBlock();

  std::cout << "number of vertices :  " ;
  std::cout << t.number_of_vertices() << std::endl;
  std::cout << "number of faces :  " ;
  std::cout << t.number_of_faces() << std::endl;
  
  trace.beginBlock("Area minimizing triangulation");
  Edge_iterator itnext;
  for( Edge_iterator it = t.edges_begin(), itend=t.edges_end();
       it != itend; it = itnext )
    {
      // vertex(cw(i)) and vertex(ccw(i)) of f.
      itnext = it; ++itnext;
      Edge e1 = *it;
      Edge e2 = t.mirror_edge( e1 );
      bool empty_f1 = emptyLatticeTriangle( t, e1.first );
      bool empty_f2 = emptyLatticeTriangle( t, e2.first );
      if ( ( empty_f1 == false )
           && ( empty_f2 == false ) )
        { // try if flip is better.
          bool empty_flip_f1 
            = emptyLatticeTriangle( t,
                                    e1.first->vertex( e1.second ), 
                                    e1.first->vertex( t.ccw( e1.second ) ),
                                    e2.first->vertex( e2.second ) );
          bool empty_flip_f2 
            = emptyLatticeTriangle( t,
                                    e2.first->vertex( e2.second ), 
                                    e2.first->vertex( t.ccw( e2.second ) ),
                                    e1.first->vertex( e1.second ) );
          if ( empty_flip_f1 || empty_flip_f2 )
            {
              std::cout << "flipped " << e1.first->vertex( e1.second )->point()
                        << "->" << e1.first->vertex( e1.second )->point()
                        << std::endl;
              //t.flip( e1.first, e1.second );
            }
        }
    }  
  trace.endBlock();


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
          board.setPenColor(DGtal::Color::Blue);
          board.setFillColor( DGtal::Color::None );
          board.setLineWidth( 3.0 );
          board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
        }
      else
        {
          board.setPenColor(DGtal::Color::Red);
          board.setFillColor( DGtal::Color::None );
          //          board.setFillColorRGBi(200,200,200,128);
          board.setLineWidth( 2.0 );
          board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
        }
    }

  Segmentation::SegmentComputerIterator i = theSegmentation.begin();
  Segmentation::SegmentComputerIterator end = theSegmentation.end();
  board.setPenColor(DGtal::Color::Green);
  board.setFillColor( DGtal::Color::None );
  board << SetMode( "ArithmeticalDSS", "BoundingBox" );
  std::string aStyleName = "ArithmeticalDSS/BoundingBox";
  for ( ; i != end; ++i) {
    DSS4 current(*i);
    board << CustomStyle( aStyleName, 
                          new CustomPenColor( DGtal::Color::Green ) )
          << current;
  } 

  for(Edge_iterator it = t.edges_begin(), itend=t.edges_end();
      it != itend; ++it)
    {
      // vertex(cw(i)) and vertex(ccw(i)) of f.
      Face_handle itf = it->first;
      int i = it->second;
      Z2i::Point a( toDGtal(itf->vertex( t.cw( i ) )->point()));
      Z2i::Point b( toDGtal(itf->vertex( t.ccw( i ) )->point()));
      std::cout << "Edge: " << a << " -> " << b << std::endl;

      CGAL::Object o = t.dual( it );
      if (CGAL::object_cast<K::Segment_2>(&o)) 
        {
          const K::Segment_2* ptrSegment = CGAL::object_cast<K::Segment_2>(&o);
          board.setPenColor(DGtal::Color::Black);
          board.setFillColor( DGtal::Color::None );
          board.setLineWidth( 2.0 );
          board.drawLine( ptrSegment->source().x(),
                          ptrSegment->source().y(),
                          ptrSegment->target().x(),
                          ptrSegment->target().y() );
        }
      else if (CGAL::object_cast<K::Ray_2>(&o)) 
        {
          const K::Ray_2* ptrRay = CGAL::object_cast<K::Ray_2>(&o);
          board.setPenColor(DGtal::Color::Black);
          board.setFillColor( DGtal::Color::None );
          board.setLineWidth( 2.0 );
          double dx = ptrRay->to_vector().x();
          double dy = ptrRay->to_vector().y();
          double norm = sqrt( dx*dx+dy*dy );
          dx = 5.0 * dx / norm;
          dy = 5.0 * dy / norm;
          board.drawArrow( ptrRay->source().x(),
                           ptrRay->source().y(),
                           ptrRay->source().x() + dx, //1*ptrRay->to_vector().x(),
                           ptrRay->source().y() + dy ); //1*ptrRay->to_vector().y() );
        }
    }
  //board << dig.getDomain();

  
  board.saveSVG("delaunay.svg");
  board.saveEPS("delaunay.eps");
  
  return 0;
}
