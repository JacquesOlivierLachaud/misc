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
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

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
  typedef typename Triangulation::Edge                    Edge;
  typedef typename Triangulation::Finite_faces_iterator   FacesIterator;
  typedef typename Triangulation::Face_handle             FaceHandle;
  typedef typename Triangulation::Point                   Point;
  typedef typename Kernel::Vector_2                  Vector;
  typedef int                                        Label;
  typedef typename std::map<VertexHandle,Label>      VertexLabeling;

private:
  Triangulation _T;
  VertexLabeling _vLabeling;
  std::set<Point> _points;

public:

  /// Constructor.
  DAC() {}

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

  /// The incident face is considered ccw.
  VertexHandle source( const Edge & e ) const
  {
    return e.first->vertex( _T.ccw( e.second ) );
  }
  
  /// The incident face is considered ccw.
  VertexHandle target( const Edge & e ) const
  {
    return e.first->vertex( _T.cw( e.second ) );
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

  
};

template <typename DigitalAffineComplex>
class ViewerDAC
{
public:
  typedef typename DigitalAffineComplex::Triangulation Triangulation;
  typedef typename Triangulation::Point         CPoint;
  typedef typename Triangulation::Edge          Edge;
  typedef typename Triangulation::Edge_iterator EdgeIterator;
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
        DPoint a = toDGtal( _dac.source( edge )->point() );
        DPoint b = toDGtal( _dac.target( edge )->point() );
        int l1 = _dac.labeling().find( _dac.source( edge ) )->second;
        int l2 = _dac.labeling().find( _dac.target( edge ) )->second;
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
  
};


int main( int argc, char** argv )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef DAC<K> DigitalAffineComplex;
  typedef DigitalAffineComplex::Point Point;
  typedef DGtal::Z2i::Domain Domain;
  typedef DGtal::Z2i::Point DPoint;

  std::vector<Point> p;
  DigitalAffineComplex dac;
  p.push_back( Point( 0,0 ) );
  p.push_back( Point( 1,1 ) );
  p.push_back( Point( 2,1 ) );
  p.push_back( Point( 3,1 ) );
  p.push_back( Point( 4,0 ) );
  p.push_back( Point( 3,-1 ) );
  p.push_back( Point( 2,-1 ) );
  p.push_back( Point( 1,-1 ) );
  dac.add( p.begin(), p.end(), 1, 0 );
  p.clear();
  p.push_back( Point( 1,0 ) );
  p.push_back( Point( 2,0 ) );
  p.push_back( Point( 3,0 ) );
  p.push_back( Point( 4,1 ) );
  p.push_back( Point( 5,1 ) );
  p.push_back( Point( 5,0 ) );
  dac.add( p.begin(), p.end(), 0, 1 );

  DGtal::Board2D board;
  ViewerDAC<DigitalAffineComplex> viewer( board, dac );
  board << Domain( DPoint( -2, -2 ), DPoint( 6, 2 ) );
  viewer.viewAll();
  board.saveSVG("dac.svg");
  board.saveEPS("dac.eps");
  return 0;
}
