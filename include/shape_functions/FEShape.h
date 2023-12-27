/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

         Definition of parent class 'FEShape' from which inherit all finite 
                               element shape classes

  ==============================================================================*/


#ifndef __FE_SHAPE_H
#define __FE_SHAPE_H

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>

#include "OFELI_Config.h"
#include "OFELIException.h"
#include "linear_algebra/Point.h"
#include <algorithm>
using std::vector;


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file FEShape.h
 *  \brief Definition file for class FEShape.
 */

//template<class T_> struct Point;
class Element;
class Side;
class Edge;

/*! \class FEShape
 *  \ingroup Shape
 *  \brief Parent class from which inherit all finite element shape classes
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
class FEShape
{

 public:

/// \brief Default Constructor
    FEShape() { }

/// \brief Constructor for an element
/// @param [in] el Pointer to element
    FEShape(const Element* el);

/// \brief Constructor for a side
/// @param [in] sd Pointer to side
    FEShape(const Side* sd);

/// \brief Destructor
    virtual ~FEShape() { }

// Initialize local point coordinates in element.
// <tt>s</tt> is a point in the reference element
// This function computes jacobian, shape functions and their partial
// derivatives at <tt>s</tt>. Other member functions only return these values.
//    void setLocal(const Point<double> &s);

/// \brief Return shape function of node <tt>i</tt> at given point.
    real_t Sh(size_t i) const { return _sh[i-1]; }

/** \brief Calculate shape function of node <tt>i</tt> at a given point <tt>s</tt>.
 *  @param [in] i Local node label
 *  @param [in] s Point in the reference triangle where the shape function is evaluated
 */
    real_t Sh(size_t        i,
              Point<real_t> s) const
       { s = 0; return _sh[i-1]; }

/** \brief Return determinant of jacobian.
 *  \details If the transformation (Reference element -> Actual element) is not affine, member function \b setLocal()
 *  must have been called before in order to calcuate relevant quantities.
 */
    real_t getDet() const { return _det; }

/// \brief Return coordinates of center of element
    Point<real_t> getCenter() const { return _c; }

/** \brief Localize a point in the element.
 *  \details Return actual coordinates in the reference element.
 *  If the transformation (Reference element -> Actual element) is not affine, member function \b setLocal()
 *  must have been called before in order to calcuate relevant quantities.
 */
    Point<real_t> getLocalPoint() const
    {
       Point<real_t> s = 0.;
       for (size_t i=0; i<_x.size(); i++)
          s += _sh[i] * _x[i];
       return s;
    }

/// \brief Localize a point in the element.
/// \details Return actual coordinates where <tt>s</tt> are coordinates in the reference element.
    Point<real_t> getLocalPoint(const Point<real_t>& s) const
    {
       Point<real_t> t = 0.;
       for (size_t i=0; i<_x.size(); i++)
          t += Sh(i+1,s) * _x[i];
       return t;
    }

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   bool                   _localized;
   real_t                 _det, _measure;
   const Element*         _el;
   const Side*            _sd;
   const Edge*            _ed;
   vector<real_t>         _sh;
   vector<size_t>         _node;
   vector<Point<real_t> > _x, _dshl, _dsh;
   Point<real_t>          _c;
   size_t                 _label;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};


/*! \class triangle
 *  \ingroup Shape
 *  \brief Defines a triangle.
 *  The reference element is the rectangle triangle with two unit edges.
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */
class triangle : public FEShape
{

 public:

/// \brief Default Constructor
    triangle()
    {
    }

/// \brief Constructor for an element.
/// \details The constructed triangle is an element in a 2-D mesh.
    triangle(const Element* el);

/// \brief Constructor for a side.
/// \details The constructed triangle is a side in a 3-D mesh.
    triangle(const Side* sd);

/// \brief Destructor
    virtual ~triangle() { }

/// \brief Return element area.
    real_t getArea() { return _area=0.5*_det; }

/// \brief Return coordinates of center of element.
    Point<real_t> getCenter() const { return _c; }

/// \brief Return coordinates of circumcenter of element.
    Point<real_t> getCircumcenter() const
    {
       real_t z = 0.5/((_x[0].x-_x[2].x)*(_x[2].y-_x[1].y) - (_x[2].x-_x[1].x)*(_x[0].y-_x[2].y));
       z = -0.5/_det;
       Point<real_t> a;
       a.x = -z*(_x[0].NNorm()*(_x[1].y-_x[2].y) + _x[1].NNorm()*(_x[2].y-_x[0].y) + _x[2].NNorm()*(_x[0].y-_x[1].y));
       a.y =  z*(_x[0].NNorm()*(_x[1].x-_x[2].x) + _x[1].NNorm()*(_x[2].x-_x[0].x) + _x[2].NNorm()*(_x[0].x-_x[1].x));
       return a;
   }

/// \brief Return radius of circumscribed circle of triangle.
    real_t getCircumRadius() const
    {
       return 0.25*_h1*_h2*_h3/_area;
    }

/// \brief Return radius of inscribed circle of triangle.
    real_t getInRadius() const
    {
       return 2*_area/(_h1+_h2+_h3);
    }

/// \brief Return reference coordinates of a point <tt>x</tt> in element.
    Point<real_t> getRefCoord(const Point<real_t>& x) const
    {
      real_t s = (_x[2].y-_x[0].y)*(x.x-_x[0].x)+(_x[0].x-_x[2].x)*(x.y-_x[0].y);
      real_t t = (_x[0].y-_x[1].y)*(x.x-_x[0].x)+(_x[1].x-_x[0].x)*(x.y-_x[0].y);
      return Point<real_t> (s/_det,t/_det);
    }

/// \brief Return maximal edge length of triangle
    real_t getMaxEdgeLength() const
    {
       return std::max(_h1,std::max(_h2,_h3));
    }

/// \brief Return minimal edge length of triangle
    real_t getMinEdgeLength() const
    {
       return std::min(_h1,std::min(_h2,_h3));
    }

/// \brief Check whether point <tt>x</tt> is in current triangle or not
    bool isIn(const Point<real_t> &x) const
    {
       Point<real_t> s = getRefCoord(x);
       if (s.x>=0.0 && s.x<=1.0 && s.y>=0.0 && (fabs(s.x+s.y-1.0)<=OFELI_EPSMCH || s.x+s.y<=1.))
          return true;
       return false;
    }

/// \brief Check whether point <tt>x</tt> is strictly in current triangle (not on the boundary) or not
    bool isStrictlyIn(const Point<real_t>& x) const
    {
       Point<real_t> s = getRefCoord(x);
       if (s.x>0 && s.x<1 && s.y>0 && s.x+s.y<1)
          return true;
       return false;
    }

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   real_t        _area, _h1, _h2, _h3;
   Point<real_t> _c;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
