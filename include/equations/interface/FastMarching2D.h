/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

           Definition of Template Class FastMarching2D for Fast Marching
                 algorithm on 2-D structured uniform grids

  ==============================================================================*/


#ifndef __FAST_MARCHING_2D_H
#define __FAST_MARCHING_2D_H

#include <stdlib.h>
#include <math.h>

#include <valarray>
using std::valarray;

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Point2D.h"
#include "mesh/Grid.h"

/*! \file FastMarching2D.h
 *  \brief Definition file for class FastMarching2D.
 */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct Pt
{
   size_t ix, iy;
   real_t val;
   Pt() { }
   Pt(int a) { ix = iy = static_cast<size_t>(a); }
   Pt & operator=(int a) { ix = iy = 0; val = a; return *this; }
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \class FastMarching2D
 *  \ingroup Interface
 *  \brief To run a Fast Marching Method on 2-D structured uniform grids.
 *
 *  This class enables running a Fast Marching procedure to calculate the signed distance function
 *  and extend a given front speed.
 *
 * \author M. Sylla, B. Meden
 * \copyright GNU Lesser Public License
 */

class FastMarching2D
{
 public:

/// \brief Default constructor
    FastMarching2D();

/** \brief Constructor using grid and level set function.
 *  @param [in] g Instance of class Grid
 *  @param [in] ls Vector containing the level set function at grid nodes. The values are
 *  \a 0 on the interface (from which the distance is computed), positive on one side and 
 *  negative on the other side. They must contain the signed distance on the nodes surrounding
 *  the interface but can take any value on other nodes, provided they have the right sign.
 */
    FastMarching2D(const Grid&   g,
                   Vect<real_t>& ls);

/** \brief Constructor using grid, level set function and velocity to extend.
 *  @param [in] g Instance of class Grid
 *  @param [in] ls Vector containing the level set function at grid nodes. The values are
 *  0 on the interface (from which the distance is computed), positive on one side and 
 *  negative on the other side. They must contain the signed distance on the nodes surrounding
 *  the interface but can take any value on other nodes, provided their sign is right.
 *  @param [in] F Vector containing the front speed at grid nodes. Only values on nodes surrounding the
 *  interface are relevant.
 */
    FastMarching2D(const Grid&   g,
                   Vect<real_t>& ls,
                   Vect<real_t>& F);

/// \brief Destructor
    ~FastMarching2D();

/** \brief Execute Fast Marching Procedure.
 *  \details Once this function was called, the vector <tt>ls</tt> used in the constructor will
 *  contain the signed distance function and <tt>F</tt> will contain the extended speed.
 */
    void execute();

/// \brief Check distance function
    void Check();

 private:
    bool _ext;
    size_t _nx, _ny;
    int _bw, _current;
    real_t _hx, _hy;
    Grid _g;
    const Grid *_cg;
    Vect<size_t> _pos;
    Vect<real_t> *_A, *_v, _ls, _f, _u, _sol;
    Vect<Pt> _heap;

    void Init();

//  Function defining obstacle
    real_t Obstacle(Point2D<real_t> a);

//  Return distance between point x and the segment joining point y to z.
    real_t Dist(const Point2D<real_t>& x,
                const Point2D<real_t>& y,
                const Point2D<real_t>& z);

//  Add an element to the heap and reorganize it
    void Add(real_t Uij,
             size_t a,
             size_t b);

//  Modify an element of the heap
    void Modify(real_t Uij,
                size_t a,
                size_t b);

//  Delete the first element in the heap that has the minimal value
    Pt Delete();

//  Return the positive solution of a second degree equation
    inline real_t PositiveSol(real_t a,
                              real_t b,
                              real_t c)
    {
       real_t delta = b*b - 4*a*c;
       if (delta == 0)
          return(-0.5*b/a);
       else {
          if (delta < 0) {
             cerr << "Error in FastMarching2D: No positive solution found." << endl;
             exit(1);
          }
          else
             return (0.5*(sqrt(delta)-b)/a);
       }
    }

//  Distance between selected points and the boundary
    void ExtendLocalVelocity(Vect<real_t>& dis);

    void Int(size_t           k,
             size_t           l,
             Vect<real_t>&    a,
             Point2D<real_t>& x1,
             Point2D<real_t>& x2);

//  Interpolate a function on a segment
    real_t Interp(Point2D<real_t>& xx);

//  Set to 'gray' the neighbors of (i,j)
    void setGray(size_t i, size_t j);
    void setGrayWithObstacle(size_t i, size_t j);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
