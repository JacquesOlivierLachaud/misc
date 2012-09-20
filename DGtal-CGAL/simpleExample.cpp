#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/Shapes.h>
#include <DGtal/shapes/ShapeFactory.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/geometry/curves/GridCurve.h>
#include <DGtal/io/boards/Board2D.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef Delaunay::Vertex_circulator Vertex_circulator;
typedef Delaunay::Edge_iterator  Edge_iterator;
typedef Delaunay::Finite_faces_iterator  Faces_iterator;
typedef Delaunay::Point             Point;

int main ()
{
  
  Delaunay t;
  
  typedef Flower2D<Z2i::Space> Flower; 
  Flower2D<Z2i::Space> flower(Z2i::Point(0,0), 200, 50, 5, 0);
  double h = 1; 
  GaussDigitizer<Z2i::Space,Flower> dig;  
  dig.attach( flower );
  dig.init( flower.getLowerBound()+Z2i::Vector(-1,-1),
            flower.getUpperBound()+Z2i::Vector(1,1), h ); 
  Z2i::KSpace ks;
  ks.init( dig.getLowerBound(), dig.getUpperBound(), true );
  SurfelAdjacency<2> sAdj( true );
  Z2i::SCell bel = Surfaces<Z2i::KSpace>::findABel( ks, dig, 1000 );
  vector<Z2i::Point> boundaryPoints;
  Surfaces<Z2i::KSpace>
    ::track2DBoundaryPoints( boundaryPoints, ks, sAdj, dig, bel );
  Z2i::Curve c;
  c.initFromVector( boundaryPoints );  
  typedef Z2i::Curve::PointsRange Range; 
  Range r = c.getPointsRange(); 
  
  for(Range::ConstIterator it=r.begin(), itend=r.end(); it != itend;
      ++it)
    t.insert( Point( (*it)[0], (*it)[1]));
  

  Vertex_circulator vc = t.incident_vertices(t.infinite_vertex()),
    done(vc);
  if (vc != 0) 
    {
      do 
	{ 
	  std::cout << vc->point() << std::endl;
	}
      while(++vc != done);
    }


  std::cout << "number of vertices :  " ;
  std::cout << t.number_of_vertices() << std::endl;
  std::cout << "number of faces :  " ;
  std::cout << t.number_of_faces() << std::endl;
  
  DGtal::Board2D board;
  
  for(Range::ConstIterator it=r.begin(), itend=r.end(); it != itend;
      ++it)
    board << *it;
  
  Z2i::Point dP;
  board << CustomStyle( dP.className(), 
                        new CustomPen( Color(0,0,0), Color(180,18,18), 1, 
                                       Board2D::Shape::SolidStyle,
                                       Board2D::Shape::RoundCap,
                                       Board2D::Shape::RoundJoin ));
    
  board.setPenColor(DGtal::Color::Red);
  board.setFillColorRGBi(255,255,255,0);
  for(Faces_iterator it = t.finite_faces_begin(), itend=t.finite_faces_end();
      it != itend; ++it)
    {
      Z2i::Point a( it->vertex(0)->point().x(),
		    it->vertex(0)->point().y()),
	b(it->vertex(1)->point().x(),
	  it->vertex(1)->point().y()),
	c(it->vertex(2)->point().x(),
	  it->vertex(2)->point().y());
      
      
      board.drawTriangle(a[0],a[1],b[0],b[1],c[0],c[1]);
    }
  
  board.saveSVG("delaunay.svg");
  
  return 0;
}
