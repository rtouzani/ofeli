/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

               Definition of some auxiliary optimization functions

  ==============================================================================*/

#ifndef __OPTIM_AUX_H
#define __OPTIM_AUX_H

#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file OptimAux.h
 *  \ingroup Solver
 *  \brief File that contains auxiliary functions for optimization.
 */

/** \fn int void BCAsConstraint(const Mesh &m, const Vect<double> &bc, Vect<double> &up, Vect<double> &low)
 *  \ingroup Solver
 *  \brief To impose Dirichlet boundary conditions in an optimization problem.
 *  If such conditions are to present, this function has to be invoked 
 *  by giving on input <tt>bc(i)</tt> as the value to impose for the
 *  <tt>i</tt>-th optimization variable.
 *  @param [in] m Mesh instance
 *  @param [in] bc Vect instance where <tt>bc(i)</tt> is the value to impose for dof <tt>i</tt>
 *  @param [out] up Vect instance that contains on output upper bounds for DOFs
 *  @param [out] low Vect instance that contains on output lower bounds for DOFs
 */
inline void BCAsConstraint(const Mesh&         m,
                           const Vect<real_t>& bc,
                                 Vect<real_t>& up,
                                 Vect<real_t>& low)
{
   for (size_t i=0; i<low.size(); i++) {
      low[i] = -1.e38;
      up[i] = 1.e38;
   }
   mesh_nodes(m) {
      for (size_t j=1; j<=The_node.getNbDOF(); j++) {
          size_t k = The_node.getDOF(j);
          if (The_node.getCode(j)>0)
             low(k) = up(k) = bc(k);
      }
   }
}

} /* namespace OFELI */

#endif
