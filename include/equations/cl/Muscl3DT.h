/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                            Definition of Class Muscl3DT
     Class for Muscl finite volumes for 3-D problems using tetrahedral meshes

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/

#ifndef __MUSCL3DT_H
#define __MUSCL3DT_H


#include <iostream>
using std::endl;

#include <valarray>
using std::valarray;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include "equations/cl/Muscl.h"
#include "mesh/MeshUtil.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Muscl3DT.h
 *  \brief Definition file for class Muscl3DT
 */

template<class T_,size_t N_> class LocalVect;
template<class T_> class Vect;


/*! \class Muscl3DT
 *  \ingroup ConservationLaws
 *  \brief Class for 3-D hyperbolic solvers with Muscl scheme using tetrahedra
 *
 * \author S. Clain, V. Clauzon
 * \copyright GNU Lesser Public License
 */

class Muscl3DT : public Muscl {

  public:

   using Muscl::_nb_sides;
   using Muscl::_nb_elements;

/// \brief Constructor using mesh
    Muscl3DT(Mesh& m);

/// \brief Destructor
    ~Muscl3DT();
   
/** \brief Function to reconstruct by the Muscl method
 *  @param [in] U Field to reconstruct
 *  @param [out] LU Left gradient vector
 *  @param [out] RU Right gradient vector
 *  @param [in] dof Label of dof to reconstruct
 */
    bool setReconstruction(const Vect<real_t>& U,
                                 Vect<real_t>& LU,
                                 Vect<real_t>& RU,
                                 size_t        dof);

/// \brief Return minimum area of faces in the mesh
    real_t getMinimumFaceArea() const { return _MinimumFaceArea; }

/// \brief Return minimum volume of elements in the mesh
    real_t getMinimumElementVolume() const { return _MinimumElementVolume; }

/// \brief Return maximum area of faces in the mesh
    real_t getMaximumFaceArea() const { return _MaximumFaceArea; }

/// \brief Return maximum volume of elements in the mesh
    real_t getMaximumElementVolume() const { return _MaximumElementVolume; }

/// \brief Return mean area of faces in the mesh
    real_t getMeanFaceArea() const { return _MeanFaceArea; }

/// \brief Return mean volume of elements in the mesh
    real_t getMeanElementVolume() const { return _MeanElementVolume; }

/// \brief Return minimum length of edges in the mesh
    real_t getMinimumEdgeLength() const { return _MinimumEdgeLength; }

/// \brief Return minimum volume by area in the mesh
    real_t getMinimumVolumebyArea() const { return _MinimumVolumebyArea; }

/// \brief Return maximum length of edges in the mesh
    real_t getMaximumEdgeLength() const { return _MaximumEdgeLength; }

/// \brief Return value of \a tau lim
    real_t getTauLim() const { return _taulim; }

/// \brief Return value of \a Comega
    real_t getComega() const { return _Comega; }

/// \brief Assign value of \a beta lim
    void setbetalim(real_t bl) { _betalim = bl; }

  protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// normals
   Vect<Point<real_t> > _n; 

// Ki/Sij
   Vect<real_t> _Lrate, _Rrate;

/// this is the main function here, compute all geometric stuff. Could take a while
/// implementation of the virtual function defined in Muscl.h
    void Initialize();

    void FirstOrder(const Vect<real_t>& U,
                          Vect<real_t>& LU,
                          Vect<real_t>& RU,
                          size_t        dof);
    void MultiSlopeM(const Vect<real_t>& U,
                           Vect<real_t>& LU,
                           Vect<real_t>& RU,
                           size_t        dof);
    void MultiSlopeQ(const Vect<real_t>& U,
                           Vect<real_t>& LU,
                           Vect<real_t>& RU,
                           size_t        dof);
    void GradientQ(const Vect<real_t>& U,
                         Vect<real_t>& LU,
                         Vect<real_t>& RU,
                         size_t        dof);
    void GradientM(const Vect<real_t>& U,
                         Vect<real_t>& LU,
                         Vect<real_t>& RU,
                         size_t        dof);
    void MultiSlopeMStable(const Vect<real_t>& U,
                                 Vect<real_t>& LU,
                                 Vect<real_t>& RU,
                                 size_t        dof);

    int trick_i(int i)
    {
       if (i<5)
          return i;
       else
          return (i%5)+1;
    }

//  get gradient using 4 data and geometry context of an Element 
    Point<real_t> getGrad(const LocalVect<real_t,4>& qU,
                                Element*             el);

//  algebra
    Point<real_t> getDetCrammer_33(const LocalMatrix<real_t,3,3>& A,
                                   const Point<real_t>&           B);

//  B_iB_j and B_iQ_j
    Vect<Point<real_t> > _BBv, _BQv;
//  norm of  BiMij,  B_iB_j and B_iQ_j 
    Vect<real_t> _BMl, _BBl, _BQl; 
//  B_iM_j
    Vect<Point<real_t> > _BMv; 
//  projection coeff for multislope methods
    Vect<real_t> _beta12, _beta13, _beta14, _betabis12, _betabis13, _betabis14;
//  projection coeff for multislpoe M method
    Vect<real_t> _mu1, _mu2, _mu3, _mu4;
//  If a cell is not multislope convenient, this variable is set to 0
    Vect<int> _order;

// mesh characteristics
   real_t _MinimumFaceArea, _MinimumElementVolume, _MaximumFaceArea, _MaximumElementVolume;
   real_t _MeanFaceArea, _MeanElementVolume, _MinimumEdgeLength, _MinimumVolumebyArea;
   real_t _MaximumEdgeLength, _taulim, _Comega, _betalim;

   real_t getMinSize32();
   real_t gettaulim();
   real_t getComega();
   void getgraphtau();
   void getgraphComega();

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};


//-----------------------------------------------------------------------------
// Associated functions
//-----------------------------------------------------------------------------

/// \fn ostream & operator<<(ostream& s, const Muscl3DT &m)
/// \brief Output mesh data as calculated in class Muscl3DT.
/// \ingroup Solver
    ostream & operator<<(ostream& s, const Muscl3DT &m);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
