/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

              Definition of class Grid for structured 2-D and 3-D grids

  ==============================================================================*/

#ifndef __GRID_H
#define __GRID_H

#include <stdlib.h>
#include <math.h>
#include <string>

#include <iostream>
using std::cout;
using std::ostream;
using std::endl;

#include <iomanip>
using std::setw;
using std::string;

#include "OFELI_Config.h"
#include "linear_algebra/Point.h"
#include "linear_algebra/Point2D.h"
#include "mesh/Node.h"
#include "mesh/Element.h"
#include "mesh/Side.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Grid.h
 *  \brief Definition file for class Grid.
 */

/*! \class Grid
 *  \ingroup Mesh
 *  \brief To manipulate structured grids.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Grid
{
  public:
    
/// \brief Construct a default grid with 10 intervals in each direction
    Grid();

/** \brief Construct a 1-D structured grid given its extremal coordinates and number of intervals
 *  @param [in] xm Minimal value for <tt>x</tt>
 *  @param [in] xM Maximal value for <tt>x</tt>
 *  @param [in] npx Number of grid intervals in the <tt>x</tt>-direction
 */
    Grid(real_t xm,
         real_t xM,
         size_t npx);

/** \brief Construct a 2-D structured grid given its extremal coordinates and number of intervals
 *  @param [in] xm Minimal value for <tt>x</tt>
 *  @param [in] xM Maximal value for <tt>x</tt>
 *  @param [in] ym Minimal value for <tt>y</tt>
 *  @param [in] yM Maximal value for <tt>y</tt>
 *  @param [in] npx Number of grid intervals in the <tt>x</tt>-direction
 *  @param [in] npy Number of grid intervals in the <tt>y</tt>-direction
 */
    Grid(real_t xm,
         real_t xM,
         real_t ym,
         real_t yM,
         size_t npx,
         size_t npy);

/** \brief Construct a 2-D structured grid given its extremal coordinates and number of intervals
 *  @param [in] m Minimal coordinate value
 *  @param [in] M Maximal coordinate value
 *  @param [in] npx Number of grid intervals in the <tt>x</tt>-direction
 *  @param [in] npy Number of grid intervals in the <tt>y</tt>-direction
 */
    Grid(Point<real_t> m,
         Point<real_t> M,
         size_t        npx,
         size_t        npy);

/** \brief Construct a 3-D structured grid given its extremal coordinates and number of intervals
 *  @param [in] xm Minimal value for <tt>x</tt>
 *  @param [in] xM Maximal value for <tt>x</tt>
 *  @param [in] ym Minimal value for <tt>y</tt>
 *  @param [in] yM Maximal value for <tt>y</tt>
 *  @param [in] zm Minimal value for <tt>z</tt>
 *  @param [in] zM Maximal value for <tt>z</tt>
 *  @param [in] npx Number of grid intervals in the <tt>x</tt>-direction
 *  @param [in] npy Number of grid intervals in the <tt>y</tt>-direction
 *  @param [in] npz Number of grid intervals in the <tt>z</tt>-direction
 */
    Grid(real_t xm,
         real_t xM,
         real_t ym,
         real_t yM,
         real_t zm,
         real_t zM, 
         size_t npx,
         size_t npy,
         size_t npz);

/** \brief Construct a 3-D structured grid given its extremal coordinates and number of intervals
 *  @param [in] m Minimal coordinate value
 *  @param [in] M Maximal coordinate value
 *  @param [in] npx Number of grid intervals in the <tt>x</tt>-direction
 *  @param [in] npy Number of grid intervals in the <tt>y</tt>-direction
 *  @param [in] npz Number of grid intervals in the <tt>z</tt>-direction
 */
    Grid(Point<real_t> m, Point<real_t> M, size_t npx, size_t npy, size_t npz);

/// \brief Destructor
    ~Grid() { }

/// \brief Set min. coordinates of the domain.
/// @param [in] x Minimal values of coordinates
    void setXMin(const Point<real_t>& x) { _xmin = x; }

/// Set max. coordinates of the domain.
/// @param [in] x Maximal values of coordinates
    void setXMax(const Point<real_t>& x) { _xmax = x; }

/** \brief Set Dimensions of the domain: 1-D case
 *  @param [in] xmin Minimal value of <tt>x</tt>-coordinate
 *  @param [in] xmax Maximal value of <tt>x</tt>-coordinate
 */
    void setDomain(real_t xmin,
                   real_t xmax);

/** \brief Set Dimensions of the domain: 2-D case
 *  @param [in] xmin Minimal value of <tt>x</tt>-coordinate
 *  @param [in] xmax Maximal value of <tt>x</tt>-coordinate
 *  @param [in] ymin Minimal value of <tt>y</tt>-coordinate
 *  @param [in] ymax Maximal value of <tt>y</tt>-coordinate
 */
    void setDomain(real_t xmin,
                   real_t xmax,
                   real_t ymin,
                   real_t ymax);

/** \brief Set Dimensions of the domain: 3-D case
 *  @param [in] xmin Minimal value of <tt>x</tt>-coordinate
 *  @param [in] xmax Maximal value of <tt>x</tt>-coordinate
 *  @param [in] ymin Minimal value of <tt>y</tt>-coordinate
 *  @param [in] ymax Maximal value of <tt>y</tt>-coordinate
 *  @param [in] zmin Minimal value of <tt>z</tt>-coordinate
 *  @param [in] zmax Maximal value of <tt>z</tt>-coordinate
 */
    void setDomain(real_t xmin,
                   real_t xmax,
                   real_t ymin,
                   real_t ymax,
                   real_t zmin,
                   real_t zmax);

/** \brief Set Dimensions of the domain: 3-D case
 *  @param [in] xmin Minimal coordinate value
 *  @param [in] xmax Maximal coordinate value
 */
    void setDomain(Point<real_t> xmin,
                   Point<real_t> xmax);

/// \brief Return min. Coordinates of the domain
    const Point<real_t>& getXMin() const { return _xmin; }

/// \brief Return max. Coordinates of the domain
    const Point<real_t>& getXMax() const { return _xmax; }

/** \brief Set number of grid intervals in the <tt>x</tt>, <tt>y</tt> and <tt>z</tt>-directions.
 *  \details Number of points is the number of intervals plus one in each direction
 *  @param [in] nx Number of grid intervals in the <tt>x</tt>-direction
 *  @param [in] ny Number of grid intervals in the <tt>y</tt>-direction (Default=<tt>0</tt>: 1-D grid)
 *  @param [in] nz Number of grid intervals in the <tt>z</tt>-direction (Default=<tt>0</tt>: 1-D or 2-D grid)
 *  @remark: The size of the grid (<tt>xmin</tt> and <tt>xmax</tt>) must have been defined before.
 */
    void setN(size_t nx,
              size_t ny=0,
              size_t nz=0);

/// \brief Set number of degrees of freedom for a node [Default: <tt>1</tt>]
    void setNbDOF(size_t n);

/// \brief Return number of grid intervals in the <tt>x</tt>-direction.
    size_t getNx() const { return _n.x; }

/// \brief Return number of grid intervals in the <tt>y</tt>-direction.
/// \details <tt>ny=0</tt> for 1-D domains (segments)
    size_t getNy() const { return _n.y; }

/// \brief Return number of grid intervals in the z-direction.
/// \details <tt>nz=0</tt> for 1-D (segments) and 2-D domains (rectangles)
    size_t getNz() const { return _n.z; }

/// \brief Return grid size in the x-direction
    real_t getHx() const { return _h.x; }

/// \brief Return grid size in the y-direction
    real_t getHy() const { return _h.y; }

/// \brief Return grid size in the z-direction
    real_t getHz() const { return _h.z; }
    
/// \brief Return coordinates a point with label <tt>i</tt> in a 1-D grid
    Point<real_t> getCoord(size_t i) const;
    
/// \brief Return coordinates a point with label <tt>(i,j)</tt> in a 2-D grid
    Point<real_t> getCoord(size_t i,
                           size_t j) const;

/// \brief Return coordinates a point with label <tt>(i,j,k)</tt> in a 3-D grid
    Point<real_t> getCoord(size_t i,
                           size_t j,
                           size_t k) const;

/// \brief Return total number of grid nodes
    size_t getNbNodes() const { return _nb_nodes; }

/// \brief Return total number of dof
    size_t getNbDOF() const { return _nb_nodes*_nb_dof; }

/// \brief Return x-coordinate of point with index <tt>i</tt>
    real_t getX(size_t i) const;

/// \brief Return y-coordinate of point with index <tt>j</tt>
    real_t getY(size_t j) const;

/// \brief Return z-coordinate of point with index <tt>k</tt>
    real_t getZ(size_t k) const;

/// \brief Return coordinates of point with indices <tt>(i,j)</tt>
    Point2D<real_t> getXY(size_t i,
                          size_t j) const;

/// \brief Return coordinates of point with indices <tt>(i,j,k)</tt>
    Point<real_t> getXYZ(size_t i,
                         size_t j,
                         size_t k) const;

/** \brief Return coordinates of center of a 1-D cell with indices <tt>i</tt>,
 *  <tt>i+1</tt> 
 */
    real_t getCenter(size_t i) const;

/** \brief Return coordinates of center of a 2-D cell with indices <tt>(i,j)</tt>,
 *  <tt>(i+1,j)</tt>, <tt>(i+1,j+1)</tt>, <tt>(i,j+1)</tt>
 */
    Point<real_t> getCenter(size_t i,
                            size_t j) const;

/** \brief Return coordinates of center of a 3-D cell with indices <tt>(i,j,k)</tt>,
 *  <tt>(i+1,j,k)</tt>, <tt>(i+1,j+1,k)</tt>, <tt>(i,j+1,k)</tt>, <tt>(i,j,k+1)</tt>, 
 *  <tt>(i+1,j,k+1)</tt>, <tt>(i+1,j+1,k+1)</tt>, <tt>(i,j+1,k+1)</tt> 
 */
    Point<real_t> getCenter(size_t i,
                            size_t j,
                            size_t k) const;

/** \brief Set a code for some grid points
 *  @param [in] exp Regular expression that determines the set of grid points
 *  on which the code is applied.
 *  @param [in] code Code to assign.
 */
    void setCode(string exp,
                 int    code);

/** \brief Set a code for grid points on sides
 *  @param [in] side Side for which code is assigned. Possible values are: <tt>MIN_X</tt>, 
 *  <tt>MAX_X</tt>, <tt>MIN_Y</tt>, <tt>MAX_Y</tt>, <tt>MIN_Z</tt>, <tt>MAX_Z</tt>
 *  @param [in] code Code to assign.
 */
    void setCode(int side,
                 int code);

/** \brief Return code for a side number
 *  @param [in] side Side for which code is returned. Possible values are: <tt>MIN_X</tt>, 
 *  <tt>MAX_X</tt>, <tt>MIN_Y</tt>, <tt>MAX_Y</tt>, <tt>MIN_Z</tt>, <tt>MAX_Z</tt>
 */
    int getCode(int side) const;

/** \brief Return code for a grid point
    @param [in] i <tt>i</tt>-th index for node for which code is to be returned.
    @param [in] j <tt>j</tt>-th index for node for which code is to be returned.
 */
    int getCode(size_t i,
                size_t j) const;

/** \brief Return code for a grid point
    @param [in] i <tt>i</tt>-th index for node for which code is to be returned.
    @param [in] j <tt>j</tt>-th index for node for which code is to be returned.
    @param [in] k <tt>k</tt>-th index for node for which code is to be returned.
 */
    int getCode(size_t i,
                size_t j,
                size_t k) const;

/// \brief Return space dimension
    size_t getDim() const { return _dim; }

/** \brief Change state of a cell from active to inactive (1-D grid)
 *  @param [in] i grid cell to remove
 */
    void Deactivate(size_t i);

/** \brief Change state of a cell from active to inactive (2-D grid)
 *  @param [in] i <tt>i</tt>-th index for grid cell to remove. 
 *  If this value is <tt>0</tt>, all cells <tt>(*,j)</tt> are deactivated 
 *  @param [in] j <tt>j</tt>-th index for grid cell to remove
 *  If this value is <tt>0</tt>, all cells <tt>(i,*)</tt> are deactivated
 *  @remark if <tt>i</tt> and <tt>j</tt> have value <tt>0</tt> all grid cells
 *  are deactivated !!
 */
    void Deactivate(size_t i,
                    size_t j);

/** \brief Change state of a cell from active to inactive (2-D grid)
 *  @param [in] i <tt>i</tt>-th index for grid cell to remove. 
 *  If this value is <tt>0</tt>, all cells <tt>(*,j,k)</tt> are deactivated 
 *  @param [in] j <tt>j</tt>-th index for grid cell to remove
 *  If this value is <tt>0</tt>, all cells <tt>(i,*,k)</tt> are deactivated
 *  @param [in] k <tt>k</tt>-th index for grid cell to remove
 *  If this value is <tt>0</tt>, all cells <tt>(i,j,*)</tt> are deactivated
 */
    void Deactivate(size_t i,
                    size_t j,
                    size_t k);

/** \brief Say if cell is active or not (1-D grid)
 *  @param [in] i Index of cell
 *  @return <tt>1</tt> if cell is active, <tt>0</tt> if not
 */
    int isActive(size_t i) const { return _active[i-1]; }

/** \brief Say if cell is active or not (2-D grid)
 *  @param [in] i <tt>i</tt>-th index of cell
 *  @param [in] j <tt>j</tt>-th index of cell
 *  @return <tt>1</tt> if cell is active, <tt>0</tt> if not
 */
    int isActive(size_t i,
                 size_t j) const
       { return _active[_n.y*(i-1)+j-1]; }

/** \brief Say if cell is active or not (3-D grid)
 *  @param [in] i <tt>i</tt>-th index of cell
 *  @param [in] j <tt>j</tt>-th index of cell
 *  @param [in] k <tt>k</tt>-th index of cell
 *  @return <tt>1</tt> if cell is active, <tt>0</tt> if not
 */
    int isActive(size_t i,
                 size_t j,
                 size_t k) const
       { return _active[_n.y*_n.z*(i-1)+_n.z*(j-1)+k-1]; }

  private:
    size_t _dim, _nb_nodes, _nb_dof;
    Point<real_t> _xmin, _xmax, _h;
    Point<size_t> _n;
    Point<int> _cm, _cM;
    vector<int> _active;
};

/** \fn ostream & operator<<(ostream& s, const Grid &g)
 *  \brief Output grid data.
 *  \ingroup Grid
 */
    ostream& operator<<(ostream&    s,
                        const Grid& g);

} /* namespace OFELI */

#endif
