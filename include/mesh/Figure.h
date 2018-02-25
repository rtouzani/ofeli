/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

   This file is part of OFELI.

   OFELI is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OFELI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

                    Definition of classes to define figures

  ==============================================================================*/

#ifndef __FIGURE_H
#define __FIGURE_H

#include "mesh/Domain.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Figure.h
 *  \brief Definition file for figure classes.
 */

/** \class Figure
 *  \ingroup Mesh
 *  \brief To store and treat a figure (or shape) information.
 *
 * This class is essentially useful to construct data for mesh generators
 * and for distance calculations.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class Figure
{

 public:

  enum {
    UNION         = 1,
    INTERSECTION  = 2,
    MINUS         = 3,
    SHIFT         = 4,
    ROTATION      = 5
  };

/// \brief Default constructor
    Figure() { _code = 1; _nb_figs = 0; }

/// \brief Copy constructor
    Figure(const Figure& f)
    {
       for (size_t i=0; i<10; i++)
          _v[i] = f._v[i];
       for (size_t j=0; j<20; j++) {
          _fig[j] = f._fig[j];
          _oper[j] = f._oper[j];
       }
       _code = f._code;
       _nb_figs = f._nb_figs;
       _fig[0] = &f;
    }

/// \brief Destructor
    virtual ~Figure() { }

/// \brief Choose a code for the domain defined by the figure
    void setCode(int code) { _code = code; }

/// \brief Return signed distance from a given point to current figure
/// @param [in] p Point instance from which distance is computed
    virtual real_t getSignedDistance(const Point<real_t> &p) const
    {
       real_t dist=_fig[0]->getSignedDistance(p);
       for (size_t i=1; i<_nb_figs; i++) {
          if (_oper[i-1]==UNION)
             dist = std::min(_fig[i]->getSignedDistance(p),dist);
          else if (_oper[i-1]==INTERSECTION)
             dist = std::max(_fig[i]->getSignedDistance(p),dist);
          else if (_oper[i-1]==MINUS)
             dist = std::max(-_fig[i]->getSignedDistance(p),dist);
          else if (_oper[i-1]==SHIFT)
             dist = 0;
       }
       return dist;
    }

/// \brief Operator =
    Figure& operator=(const Figure& f)
    {
       for (size_t i=0; i<10; i++)
          _v[i] = f._v[i];
       for (size_t j=0; j<20; j++) {
          _fig[j] = f._fig[j];
          _oper[j] = f._oper[j];
       }
       _code = f._code;
       _nb_figs = f._nb_figs;
       return *this;
    }

/** \brief Calculate signed distance to current figure with respect to grid points
 *  @param [in] g Grid instance
 *  @param [in] d Vect instance containing calculated distance from each grid index
 *  to Figure
 *  @remark Vector <tt>d</tt> doesn't need to be sized before invoking this function
 */
    void getSignedDistance(const Grid&   g,
                           Vect<real_t>& d) const
    {
       if (g.getDim()==2) {
          d.setSize(g.getNx()+1,g.getNy()+1);
          for (size_t i=1; i<=g.getNx()+1; i++)
             for (size_t j=1; j<=g.getNy()+1; j++)
                d(i,j) = getSignedDistance(g.getCoord(i,j));
       }
       else if (g.getDim()==3) {
          d.setSize(g.getNx()+1,g.getNy()+1,g.getNz()+1);
          for (size_t i=1; i<=g.getNx()+1; i++)
             for (size_t j=1; j<=g.getNy()+1; j++)
                for (size_t k=1; k<=g.getNz()+1; k++)
                   d(i,j,k) = getSignedDistance(g.getCoord(i,j,k));
       }
    }

    friend Figure operator&&(const Figure& f1, const Figure& f2);
    friend Figure operator||(const Figure& f1, const Figure& f2);
    friend Figure operator-(const Figure& f1, const Figure& f2);

/** \brief Compute signed distance from a line
 *  @param [in] p Point for which distance is computed
 *  @param [in] a First vertex of line
 *  @param [in] b Second vertex of line
 *  @return Signed distance
 */
    real_t dLine(const Point<real_t>& p,
                 const Point<real_t>& a,
                 const Point<real_t>& b) const
    {
       real_t rn = (p-a)*(b-a), rd=(b-a)*(b-a);
       real_t r = rn/rd;
       Point<real_t> q(a+r*(b-a));
       real_t s = ((a.y-p.y)*(b.x-a.x)-(a.x-p.x)*(b.y-a.y) ) / rd;
       real_t ds, dl = fabs(s)*sqrt(rd);
       Point<real_t> xx = q;
       if ((r>=0) && (r<=1))
          ds = dl;
       else {
          real_t dist1=SqrDistance(a,p), dist2=SqrDistance(b,p);
          if (dist1 < dist2) {
             xx = a;
             ds = sqrt(dist1);
          }
          else {
             xx = b;
             ds = sqrt(dist2);
          }
       }
       return ds;
    }

 protected:
   Point<real_t> _v[10];
   int _code;
   size_t _nb_figs;
   const Figure *_fig[20];
   int _oper[20];
};


/** \fn Figure operator&&(const Figure &f1, const Figure &f2)
 *  \brief Function to define a Figure instance as the intersection of two Figure instances
 *  \ingroup Mesh
 *  @param [in] f1 First Figure instance
 *  @param [in] f2 Second Figure instance
 *  @return Updated resulting Figure instance
 */
Figure operator&&(const Figure& f1,
                  const Figure& f2)
{
   Figure f;
   f._oper[f._nb_figs] = Figure::INTERSECTION;
   f._fig[f._nb_figs++] = &f1;
   f._fig[f._nb_figs++] = &f2;
   return f;
}


/** \fn Figure operator||(const Figure &f1, const Figure &f2)
 *  \brief Function to define a Figure instance as the union of two Figure instances
 *  \ingroup Mesh
 *  @param [in] f1 First Figure instance
 *  @param [in] f2 Second Figure instance
 *  @return Updated resulting Figure instance
 */
Figure operator||(const Figure& f1,
                  const Figure& f2)
{
   Figure f;
   f._oper[f._nb_figs] = Figure::UNION;
   f._fig[f._nb_figs++] = &f1;
   f._fig[f._nb_figs++] = &f2;
   return f;
}


/** \fn Figure operator-(const Figure& f1, const Figure& f2)
 *  \brief Function to define a Figure instance as the set subtraction of two Figure instances
 *  \ingroup Mesh
 *  @param [in] f1 First Figure instance to subtract from
 *  @param [in] f2 Second Figure instance to subtract
 *  @return Updated resulting Figure instance
 */
Figure operator-(const Figure& f1,
                 const Figure& f2)
{
   Figure f;
   f._oper[f._nb_figs] = Figure::MINUS;
   f._fig[f._nb_figs++] = &f1;
   f._fig[f._nb_figs++] = &f2;
   return f;
}


/** \class Rectangle
 *  \ingroup Mesh
 *  \brief To store and treat a rectangular figure.
 */
class Rectangle : public Figure
{

 public:

   using Figure::getSignedDistance;

/// \brief Default constructor
    Rectangle()
    {
       _v[0].x = 0; _v[0].y = 0;
       _v[1].x = 1; _v[1].y = 1;
       _nb_figs++;
    }

/** \brief Constructor
 *  @param [in] bbm Left Bottom point of rectangle
 *  @param [in] bbM Right Top point of rectangle
 *  @param [in] code Code to assign to rectangle
 */
    Rectangle(const Point<real_t>& bbm,
              const Point<real_t>& bbM,
              int                  code=1)
    {
       _nb_figs++;
       _v[0] = bbm;
       _v[1] = bbM;
       _code = code;
    }

/** \brief Assign bounding box of the rectangle
 *  @param [in] bbm Left Bottom point
 *  @param [in] bbM Right Top point
 */
    void setBoundingBox(const Point<real_t>& bbm,
                        const Point<real_t>& bbM) 
    {
      _v[0] = bbm;
      _v[1] = bbM;
    }

/// \brief Return first point of bounding box
    Point<real_t> getBoundingBox1() const { return _v[0]; }

/// \brief Return second point of bounding box
    Point<real_t> getBoundingBox2() const { return _v[1]; }

/** \brief Return signed distance of a given point from the current rectangle
 *  \details The computed distance is negative if <tt>p</tt> lies in the rectangle,
 *  negative if it is outside, and <tt>0</tt> on its boundary
 *  @param [in] p Point<double> instance
 */
    real_t getSignedDistance(const Point<real_t>& p) const
    {
       return -std::min(std::min(std::min(p.y-_v[0].y,_v[1].y-p.y),p.x-_v[0].x+p.x),_v[1].x-p.x);
    }

/// \brief Operator +=
/// \details Translate rectangle by a vector <tt>a</tt>
    Rectangle & operator+=(Point<real_t> a)
    {
       _v[0] += a;
       _v[1] += a;
       return *this;
    }

/// \brief Operator *=
/// \details Scale rectangle by a factor <tt>a</tt>
    Rectangle & operator+=(real_t a)
    {
       _v[0] *= a;
       _v[1] *= a;
       return *this;
    }

};

/** \class Brick
 *  \ingroup Mesh
 *  \brief To store and treat a brick (parallelepiped) figure.
 */
class Brick : public Figure
{

 public:

   using Figure::getSignedDistance;

/// \brief Default constructor
    Brick()
    {
       _v[0] = Point<real_t>(0.,0.,0.);
       _v[1] = Point<real_t>(1.,1.,1.);
       _nb_figs++;
    }

/** \brief Constructor
 *  @param [in] bbm first point (xmin,ymin,zmin)
 *  @param [in] bbM second point (xmax,ymax,zmax)
 *  @param [in] code Code to assign to rectangle
 */
    Brick(const Point<real_t>& bbm,
          const Point<real_t>& bbM,
          int                  code=1)
    {
       _nb_figs++;
       _v[0] = bbm;
       _v[1] = bbM;
       _code = code;
    }

/** \brief Assign bounding box of the brick
 *  @param [in] bbm first point (xmin,ymin,zmin)
 *  @param [in] bbM second point (xmax,ymax,zmax)
 */
    void setBoundingBox(const Point<real_t>& bbm,
                        const Point<real_t>& bbM)
    {
      _v[0] = bbm;
      _v[1] = bbM;
    }

/// \brief Return first point of bounding box (xmin,ymin,zmin)
    Point<real_t> getBoundingBox1() const { return _v[0]; }

/// \brief Return second point of bounding box (xmax,ymax,zmax)
    Point<real_t> getBoundingBox2() const { return _v[1]; }

/** \brief Return signed distance of a given point from the current brick
 *  \details The computed distance is negative if <tt>p</tt> lies in the brick,
 *  negative if it is outside, and 0 on its boundary
 *  @param [in] p Point<double> instance
 */
    real_t getSignedDistance(const Point<real_t>& p) const
    {
       return -std::min(std::min(std::min(std::min(std::min(p.y-_v[0].y, _v[1].y-p.y)
                                                           ,p.x-_v[0].x),_v[1].x-p.x)
                                                           ,p.z-_v[0].z),_v[1].z-p.z);
    }

/// \brief Operator +=
/// \details Translate brick by a vector <tt>a</tt>
    Brick & operator+=(Point<real_t> a)
    {
       _v[0] += a;
       _v[1] += a;
       _v[2] += a;
       return *this;
    }

/// \brief Operator *=
/// \details Scale brick by a factor <tt>a</tt>
    Brick & operator+=(real_t a)
    {
       _v[0] *= a;
       _v[1] *= a;
       _v[2] *= a;
       return *this;
    }

 private:

};


/** \class Circle
 *  \ingroup Mesh
 *  \brief To store and treat a circular figure.
 */
class Circle : public Figure
{

 public:

    using Figure::getSignedDistance;

/// Default construcor
    Circle() 
    {
       _r = 1;
       _v[0] = 0.0;
       _nb_figs++;
    }

/** \brief Constructor
 *  @param [in] c Coordinates of center of circle 
 *  @param [in] r Radius
 *  @param [in] code Code to assign to the generated domain [Default: 1]
 */
    Circle(const Point<real_t>& c,
           real_t               r,
           int                  code=1)
    {
       _nb_figs++;
       _v[0] = c;
       _r = r;
       _code = code;
    }
 
/// \brief Assign radius of circle
    void setRadius(real_t r) { _r = r; }

/// \brief Return radius of circle
    real_t getRadius() const { return _r; }

/// \brief Assign coordinates of center of circle
    void setCenter(const Point<real_t>& c) { _v[0] = c; }

/// \brief Return coordinates of center of circle
    Point<real_t> getCenter() const { return _v[0]; }

/** \brief Return signed distance of a given point from the current circle
 *  \details The computed distance is negative if <tt>p</tt> lies in the disk,
 *  positive if it is outside, and 0 on the circle
 *  @param [in] p Point<double> instance
 */
    real_t getSignedDistance(const Point<real_t>& p) const
    {
       return sqrt((p.x-_v[0].x)*(p.x-_v[0].x)+(p.y-_v[0].y)*(p.y-_v[0].y))-_r;
    }

/// \brief Operator +=
/// \details Translate circle by a vector <tt>a</tt>
    Circle & operator+=(Point<real_t> a)
    {
       _v[0] += a;
       return *this;
    }

/// \brief Operator *=
/// \details Scale circle by a factor <tt>a</tt>
    Circle & operator+=(real_t a)
    {
       _r *= a;
       return *this;
    }

 private:
    real_t _r;
};


/** \class Sphere
 *  \ingroup Mesh
 *  \brief To store and treat a sphere.
 */
class Sphere : public Figure
{

 public:

    using Figure::getSignedDistance;

/// Default construcor
    Sphere() 
    {
       _r = 1;
       _v[0] = 0.0;
       _nb_figs++;
    }

/** \brief Constructor
 *  @param [in] c Coordinates of center of sphere 
 *  @param [in] r Radius
 *  @param [in] code Code to assign to the generated sphere [Default: 1]
 */
    Sphere(const Point<real_t>& c,
                 real_t         r,
                 int            code=1)
    {
       _nb_figs++;
       _v[0] = c;
       _r = r;
       _code = code;
    }
 
/// \brief Assign radius of sphere
    void setRadius(real_t r) { _r = r; }

/// \brief Return radius of sphere
    real_t getRadius() const { return _r; }

/// \brief Assign coordinates of center of sphere
    void setCenter(const Point<real_t>& c) { _v[0] = c; }

/// \brief Return coordinates of center of sphere
    Point<real_t> getCenter() const { return _v[0]; }

/** \brief Return signed distance of a given point from the current sphere
 *  \details The computed distance is negative if <tt>p</tt> lies in the ball,
 *  positive if it is outside, and 0 on the sphere
 *  @param [in] p Point<double> instance
 */
    real_t getSignedDistance(const Point<real_t>& p) const
    {
       return sqrt( (p.x-_v[0].x)*(p.x-_v[0].x)
                   +(p.y-_v[0].y)*(p.y-_v[0].y)
                   +(p.y-_v[0].z)*(p.y-_v[0].z))-_r;
    }

/// \brief Operator +=
/// \details Translate sphere by a vector <tt>a</tt>
    Sphere & operator+=(Point<real_t> a)
    {
       _v[0] += a;
       return *this;
    }

/// \brief Operator *=
/// \details Scale sphere by a factor <tt>a</tt>
    Sphere & operator+=(real_t a)
    {
       _r *= a;
       return *this;
    }

 private:
    real_t _r;
};


/** \class Ellipse
 *  \ingroup Mesh
 *  \brief To store and treat an ellipsoidal figure.
 */
class Ellipse : public Figure
{

 public:

    using Figure::getSignedDistance;

/// \brief Default constructor
/// \details Constructs an ellipse with semimajor axis = 1, and semiminor axis =  1
    Ellipse()
    {
       _nb_figs++;
       _v[0] = Point<real_t>(0.,0.);
       _a = 2.;
       _b = 1.;
    }

/** \brief Constructor with given ellipse data
 *  @param [in] c Coordinates of center
 *  @param [in] a Semimajor axis
 *  @param [in] b Semiminor axis
 *  @param [in] code Code to assign to the generated figure [Default: 1]
 */
    Ellipse(Point<real_t> c,
            real_t        a,
            real_t        b,
            int           code=1)
    {
       _v[0] = c;
       _a = a;
       _b = b;
       _code = code;
    }

/** \brief Return signed distance of a given point from the current ellipse
 *  \details The computed distance is negative if <tt>p</tt> lies in the ellipse,
 *  positive if it is outside, and 0 on its boundary
 *  @param [in] p Point<double> instance
 */
    real_t getSignedDistance(const Point<real_t>& p) const
    {
       return sqrt((p.x-_v[0].x)*(p.x-_v[0].x)/(_a*_a)+(p.y-_v[0].y)*(p.y-_v[0].y)/(_b*_b))-1.0;
    }

/** \brief Operator <tt>+=</tt>
 *  \details Translate ellipse by a vector <tt>a</tt>
 *  @param [in] a Vector defining the translation
 */
    Ellipse & operator+=(Point<real_t> a)
    {
       _v[0] += a;
       return *this;
    }

/** \brief Operator <tt>*=</tt>
 *  \details Scale ellipse by a factor <tt>a</tt>
 *  @param [in] a Scaling value
 */
    Ellipse& operator+=(real_t a)
    {
       _a *= a;
       _b *= a;
       return *this;
    }

 private:
    real_t _a, _b;
};


/** \class Triangle 
 *  \ingroup Mesh
 *  \brief To store and treat a triangle.
 */
class Triangle : public Figure
{

 public:

    using Figure::getSignedDistance;
    using Figure::dLine;

/// \brief Default constructor
/// \details Constructs a unit %triangle with vertices <tt>(0,0)</tt>, <tt>(1,0)</tt> and <tt>(0,1)</tt>
    Triangle()
    {
       _nb_figs++;
       _v[0] = Point<real_t>(0.,0.);
       _v[1] = Point<real_t>(1.,0.);
       _v[2] = Point<real_t>(0.,1.);
       _inv_det = 2;
    }

/** \brief Constructor with vertices and code
 *  @param [in] v1 Coordinates of first vertex of %triangle
 *  @param [in] v2 Coordinates of second vertex of %triangle
 *  @param [in] v3 Coordinates of third vertex of %triangle
 *  @param [in] code Code to assign to the generated figure [Default: 1]
 *  @remark Vertices must be given in couterclockwise order
 */
    Triangle(const Point<real_t>& v1,
             const Point<real_t>& v2,
             const Point<real_t>& v3,
             int                  code=1)
    {
       _nb_figs++;
       _v[0] = v1;
       _v[1] = v2;
       _v[2] = v3;
       _inv_det = 1.0/((_v[1].x-_v[0].x)*(_v[2].y-_v[0].y) -
                       (_v[1].y-_v[0].y)*(_v[2].x-_v[0].x));
       _code = code;
    }

/// \brief Assign first vertex of triangle
    void setVertex1(const Point<real_t>& v) { _v[0] = v; }

/// \brief Assign second vertex of triangle
    void setVertex2(const Point<real_t>& v) { _v[1] = v; }

/// \brief Assign third vertex of triangle
    void setVertex3(const Point<real_t>& v) { _v[2] = v; }

/** \brief Return signed distance of a given point from the current triangle
 *  \details The computed distance is negative if <tt>p</tt> lies in the triangle,
 *  positive if it is outside, and 0 on its boundary
 *  @param [in] p Point<double> instance
 */
    real_t getSignedDistance(const Point<real_t>& p) const
    {
       real_t d=dLine(p,_v[0],_v[1]);
       d = std::min(d,dLine(p,_v[1],_v[2]));
       d = std::min(d,dLine(p,_v[2],_v[0]));
       real_t s = _inv_det*((_v[2].y-_v[0].y)*(p.x-_v[0].x)+(_v[0].x-_v[2].x)*(p.y-_v[0].y));
       real_t t = _inv_det*((_v[0].y-_v[1].y)*(p.x-_v[0].x)+(_v[1].x-_v[0].x)*(p.y-_v[0].y));
       if (s>=0 && s<=1 && t>=0 && s+t<=1)
          return -d;
       else
          return d;
    }

/// \brief Operator +=
/// \details Translate triangle by a vector <tt>a</tt>
    Triangle & operator+=(Point<real_t> a)
    {
       _v[0] += a;
       _v[1] += a;
       _v[2] += a;
       return *this;
    }

/// \brief Operator *=
/// \details Scale triangle by a factor <tt>a</tt>
    Triangle & operator+=(real_t a)
    {
       _v[0] *= a;
       _v[1] *= a;
       _v[2] *= a;
       return *this;
    }

 private:
    real_t _inv_det;

};


/** \class Polygon
 *  \ingroup Mesh
 *  \brief To store and treat a polygonal figure.
 */
class Polygon : public Figure
{

 public:

    using Figure::getSignedDistance;
    using Figure::dLine;

/// \brief Default constructor
    polygon() { }

/** \brief Constructor
 *  @param [in] v Vect instance containing list of coordinates of polygon vertices
 *  @param [in] code Code to assign to the generated domain (Default value = 1)
 */
    Polygon(const Vect<Point<real_t> >& v,
            int                         code=1)
    {
       _nb_figs++;
       _nb_vertices = v.size();
       for (size_t i=0; i<_nb_vertices; i++)
          _v[i] = v[i];
       _code = code;
    }

/// \brief Assign vertices of polygon
/// @param [in] v Vector containing vertices coordinates in counter clockwise order
    void setVertices(const Vect<Point<real_t> >& v)
    {
       for (size_t i=0; i<_nb_vertices; i++)
          _v[i] = v[i];
    }

/** \brief Return signed distance of a given point from the current polygon
 *  \details The computed distance is negative if <tt>p</tt> lies in the polygon,
 *  negative if it is outside, and 0 on its boundary
 *  @param [in] p Point<double> instance
 */
    real_t getSignedDistance(const Point<real_t>& p) const
    {
       real_t d=dLine(p,_v[0],_v[1]);
       for (size_t i=1; i<_nb_vertices-1; i++)
          d = std::min(d,dLine(p,_v[i],_v[i+1]));
       d = std::min(d,dLine(p,_v[_nb_vertices-1],_v[0]));
       return d;
    }

/// \brief Operator +=
/// \details Translate polygon by a vector <tt>a</tt>
    Polygon & operator+=(Point<real_t> a)
    {
       for (size_t i=0; i<_nb_vertices; i++)
          _v[i] += a;
       return *this;
    }

/// \brief Operator *=
/// \details Scale polygon by a factor <tt>a</tt>
    Polygon & operator+=(real_t a)
    {
       for (size_t i=0; i<_nb_vertices; i++)
          _v[i] *= a;
       return *this;
    }

 private:
    size_t _nb_vertices;

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
