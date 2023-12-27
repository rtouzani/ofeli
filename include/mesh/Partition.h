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

                           Definition of class Partition

  ==============================================================================*/

#ifndef __PARTITION_H
#define __PARTITION_H


#include "OFELI_Config.h"
#include "mesh/Mesh.h"
#include "io/output.h"

extern "C" {
#include "mesh/metis/metis.h"
}

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \class Partition
 *  \ingroup Mesh
 *  \brief To partition a finite element mesh into balanced submeshes.
 *
 * Class Partition enables partitioning a given mesh into a given number
 * of submeshes with a minimal connectivity. Partition uses the well known
 * <tt>metis</tt> library that is included in the OFELI library. A more detailed
 * description of metis can be found in the web site:\n
 * http://www.csit.fsu.edu/~burkardt/c_src/metis/metis.html
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class Partition
{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
 typedef struct { size_t side1, sd, side; } Interface;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 public:

//----------------------------   BASIC OPERATIONS   -----------------------------

/// \brief Default constructor
    Partition() : _theMesh(nullptr), _nb_submesh(1) { }

/** \brief Constructor to partition a mesh into submeshes
 *  @param [in] mesh Mesh instance
 *  @param [in] n Number of submeshes
 */
    Partition(Mesh&  mesh,
              size_t n);

/** \brief Constructor using already created submeshes
 *  @param [in] mesh Mesh instance
 *  @param [in] n Number of submeshes
 *  @param [in] epart Vector containing for each element its submesh label
 *  (Running from 0 to <tt>n-1</tt>
 */
    Partition(Mesh&        mesh,
              int          n,
              vector<int>& epart);

/// \brief Destructor
    ~Partition();

//-----------------------------   INSPECTORS  -----------------------------------

/// \brief Return number of submeshes
    size_t getNbSubMeshes() const { return _nb_submesh; }

/// \brief Return number of nodes in given submesh
    size_t getNbNodes(size_t i) const { return _theSubMesh[i]->getNbNodes(); }

/// \brief Return number of elements in given submesh
    size_t getNbElements(size_t i) const { return _theSubMesh[i]->getNbElements(); }

/// \brief Return the global Mesh instance
    Mesh *getMesh() { return _theMesh; }

/// \brief Return the submesh of label <tt>i</tt>
    Mesh *getMesh(size_t i) { return _theSubMesh[i]; }

/** \brief Return node label in subdomain by giving its label in initial mesh
 *  @param [in] sm Label of submesh
 *  @param [in] label Label of node in initial mesh
 */
    size_t getNodeLabelInSubMesh(size_t sm,
                                 size_t label) const
    { return _m2sm_node[sm][label-1]; }

/// \brief Return element label in subdomain by giving its label in initial mesh
    size_t getElementLabelInSubMesh(size_t sm,
                                    size_t label) const
    { return _m2sm_element[sm][label-1]; }

/** \brief Return node label in initial mesh by giving its label in submesh
 *  @param [in] sm Label of submesh
 *  @param [in] label Node label
 */
    size_t getNodeLabelInMesh(size_t sm,
                              size_t label) const
    { return _sm2m_node[sm][label-1]; }

/// \brief Return element label in initial mesh by giving its label in submesh
    size_t getElementLabelInMesh(size_t sm,
                                 size_t label) const
    { return _sm2m_element[sm][label-1]; }

/// \brief Return Number of interface sides for a given sub-mesh
    size_t getNbInterfaceSides(size_t sm) const { return _nb_interface_sides[sm]; }

/** \brief Return index of submesh that contains the <tt>i</tt>-th side label in sub-mesh <tt>sm</tt>
 *  @param [in] sm Submesh index
 *  @param [in] i Side label
 *  @return Index of submesh
 */
    size_t getSubMesh(size_t sm,
                      size_t i) const
    { return _interface_side[sm][i-1].sd; }

/** \brief Return reference to submesh
 *  @param [in] i Submesh index
 *  @return Reference to corresponding Mesh instance
 */
    Mesh& getSubMesh(size_t i) const { return *(_theSubMesh[i]); }

/** \brief Return i-th side label in a given submesh
 *  @param [in] sm Index of submesh
 *  @param [in] i Label of side
 */
    size_t getFirstSideLabel(size_t sm,
                             size_t i) const
    { return _interface_side[sm][i-1].side1; }

/** \brief Return side label in the neighbouring submesh corresponding to <tt>i</tt>-th side
 *  label in sub-mesh <tt>sm</tt>
 *  @param [in] sm Label of submesh
 *  @param [in] i Side label
 */
    size_t getSecondSideLabel(size_t sm,
                              size_t i) const
    { return _interface_side[sm][i-1].side; }

/** \brief Get number of connected nodes in a submesh
 *  @param [in] n Label of node for which connections are counted
 *  @param [in] s Label of submesh (starting from 0)
 */
    int getNbConnectInSubMesh(int n,
                              int s) const;

/** \brief Get number of connected nodes out of a submesh
 *  @param [in] n Label of node for which connections are counted
 *  @param [in] s Label of submesh (starting from 0)
 */
    int getNbConnectOutSubMesh(int n,
                               int s) const;

/** \brief Save a submesh in file
 *  @param [in] n Label of submesh
 *  @param [in] file Name of file in which submesh is saved
 */
    void put(size_t n,
             string file) const { _theSubMesh[n]->put(file); }

//-----------------------------   MODIFIERS  -----------------------------------

/// \brief Set Mesh instance
    void set(Mesh&  mesh,
             size_t n);

//------------------------   ASSOCIATED OPERATORS  -----------------------------

/// \brief Output class information
    friend ostream & operator<<(ostream&         s,
                                const Partition& p);

 private:

   Mesh           *_theMesh;
   vector<Mesh *> _theSubMesh;
   vector<size_t> _nb_interface_sides;
   size_t         **_sm2m_node, **_m2sm_node, **_sm2m_element, **_m2sm_element;
   int            _nb_submesh;
   Interface      **_interface_side;
   vector<int>    _npart, _epart;
   vector<vector<size_t> > _node_neig, _nnz;
   void Init();
   void Prepare();
   void NodeNeighborList();
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
