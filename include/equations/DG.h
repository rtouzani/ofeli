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

                       Definition of parent class DG
            for initializations of the Discontinuous Galerkin method
 
  ==============================================================================*/


#ifndef __DG_H
#define __DG_H

#include "OFELI_Config.h"
#include "equations/Equation.h"
//#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/SpMatrix.h"

namespace OFELI {

/*! \file DG.h
 *  \brief Definition file for class DG.
 */

/*! \class DG
 *  \ingroup DG
 *  \brief Enables preliminary operations and utilities for the 
 *  Discontinous Galerkin method
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class DG {

 public:


/** \brief Constructor with mesh and degree of the method
 *  @param [in] ms Mesh instance
 *  @param [in] degree Polynomial degree of the DG method [Default: <tt>1</tt>]
 */
    DG(Mesh &  ms,
       size_t degree=1);

/// \brief Destructor
    ~DG();

/// \brief Set matrix graph
    int setGraph();


 protected:

   Mesh                              *_theMesh;
   Element                           *_theElement;
   Side                              *_theSide;
   size_t                            _nb_dof, _nb_sdof, _nb_eq, _ne, _nf;
   real_t                            _volume, _area, _length;
   size_t                            _degree, _nb_el_dof, _nb_sd_dof, _neq;
   SpMatrix<real_t>                  _A;
   Vect<real_t>                      _b, _x, *_u, *_f, *_dbc, *_nbc;
   Vect<Point<real_t> >              _N;
   vector<std::map<size_t,size_t> >  _g2l;
   size_t II(size_t i) const { return _nb_el_dof*(_ne-1)+i; }
   size_t IJ(size_t i) const { return _nb_el_dof*(_nf-1)+i; }

   void setDGLabel();

};

} /* namespace OFELI */

#endif
