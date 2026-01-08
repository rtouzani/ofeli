/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

                           Definition of class 'MeshAdapt'

  ==============================================================================*/


#ifndef __MESH_ADAPT_H
#define __MESH_ADAPT_H

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <iostream>
using std::ostream;
using std::endl;
using std::cerr;
using std::vector;

#include <iomanip>
using std::setw;

#include <string>
using std::string;

#include "OFELI_Config.h"
#include "mesh/Domain.h"
#include "linear_algebra/Vect.h"
#include "mesh/MeshUtil.h"

#include "mesh/bamg/Meshio.h"
#include "mesh/bamg/Mesh2.h"
#include "mesh/bamg/QuadTree.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
void MeshErrorIO(ios&);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

using namespace bamg;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file MeshAdapt.h
 *  \brief Definition file for class MeshAdapt
 */

/*! \class MeshAdapt
 * \ingroup Mesh
 * \brief To adapt mesh in function of given solution
 * \details Class MeshAdapt enables modifying mesh according to a solution
 * vector defining at nodes. It concerns 2-D triangular meshes only. 
 * @remark Class MeshAdapt is mainly based on the software 'Bamg' developed
 * by F. Hecht, Universite Pierre et Marie Curie, Paris.
 * We warmly thank him for accepting incoporation of Bamg in the OFELI package 
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class MeshAdapt
{

 public:

/// \brief Default constructor
    MeshAdapt();

/// \brief Constructor using initial mesh
/// @param [in] ms Reference to initial mesh
    MeshAdapt(Mesh& ms);

/// \brief Constructor using a reference to class Domain
/// @param [in] dom Reference to Domain class
    MeshAdapt(Domain &dom);

/// \brief Destructor
    ~MeshAdapt() { }
   
/// \brief Get reference to Domain instance
    Domain &getDomain() const { return *_domain; }

/// \brief Get reference to current mesh
//    Mesh &getMesh() const { return *_ms[_iter-1]; }
    Mesh &getMesh() const { return *_theMesh; }

/// \brief Set reference to Domain instance
    void set(Domain &dom);

/// \brief Set reference to Mesh instance
    void set(Mesh &ms);

/// \brief Define label of node
    void setSolution(const Vect<real_t>& u);

/// \brief Set number of Jacobi iterations for smoothing
    void setJacobi(int n) { _nb_Jacobi = n; }

/// \brief Set number of smoothing iterations
    void setSmooth(int n) { _nb_smooth = n; }

/// \brief Metric is constructed with absolute error
    void AbsoluteError() { _abs_error = true; }

/// \brief Metric is constructed with relative error
    void RelativeError() { _abs_error = false; }

/// \brief Set error threshold for adaption
    void setError(real_t err) { _err = err; }

/// \brief Set minimal mesh size
    void setHMin(real_t h) { _hmin = h; }

/// \brief Set maximal mesh size
    void setHMax(real_t h) { _hmax = h; }
    
/// \brief Set minimal mesh size and set anisotropy
    void setHMinAnisotropy(real_t h) { _hmin_aniso = h; }

/// \brief Set relaxation parameter for smoothing
/// \details Default value for relaxation parameter is 1.8
    void setRelaxation(real_t omega) { _omega = omega; }

/// Set that adapted mesh construction is anisotropic
    void setAnisotropic() { _anisotropic = true; }

/// Set maximum ratio of anisotropy
    void MaxAnisotropy(real_t a) { _aniso_max = a; }

/**
 * \brief Change the metric such that the maximal subdivision of a background's edge
 * is bounded by the given number (always limited by 10)
 */
    void setMaxSubdiv(real_t s) { _max_subdiv = s; }

/** \brief Set maximum number of vertices
 *  \details Default value is 500000
 */
    void setMaxNbVertices(size_t n) { _max_nb_vertices = n; }

/** \brief Set ratio for a smoothing of the metric
 *  @param [in] r Ratio value.
 *  @note If <tt>r</tt> is 0 then no smoothing is performed, if <tt>r</tt> lies in
 *  <tt>[1.1,10]</tt> then the smoothing changes the metric such that the largest
 *  geometrical progression (speed of mesh size variation in mesh is bounded by
 *  <tt>r</tt>) (by default no smoothing)
 */
    void setRatio(real_t r) { _ratio = r; }

/// \brief Do not scale solution before metric computation
/// \details By default, solution is scaled (between 0 and 1)
    void setNoScaling() { _scaling = false; }

/// \brief Do not keep old vertices
/// \details By default, old vertices are kept
    void setNoKeep() { _keep_vertices = false; }

/// \brief set computation of the Hessian
    void setHessian() { _hessian = true; }

/// \brief Create mesh output file
    void setOutputMesh(string file) { _output_mesh_file = file; _set_outm = true; }

/// \brief Set Geometry file
    void setGeoFile(string file) { _geo_file = file; _set_geo = true; }

/// \brief Set error on geometry
    void setGeoError(real_t e) { _geo_err = e; }
   
/// \brief Set background mesh
    void setBackgroundMesh(string bgm) { _background_mesh_file = bgm; _set_bgm = true; }

/// \brief Split edges with two vertices on boundary
    void SplitBoundaryEdges() { _splitbedge = true; }

/// \brief Create a metric file
    void CreateMetricFile(string mf) { _ometric_file = mf; _set_ometric = true; }

/// \brief Set Metric file
    void setMetricFile(string mf) { _metric_file = mf; _set_metric = true; }

/// \brief Set solution defined on background mesh for metric construction
    void getSolutionMbb(string mbb) { _mbb_file = mbb; _set_mbb = _set_rbb = true; }

/// \brief Set solution defined on background mesh for metric construction
    void getSolutionMBB(string mBB) { _mBB_file = mBB; _set_mBB = true; }

/// \brief Read solution defined on the background mesh in <tt>bb</tt> file
/// \details Solution is interpolated on created mesh
    void getSolutionbb(string rbb) { _rbb_file = rbb; _set_rbb = true; }
    
/// \brief Read solution defined on the background mesh in <tt>BB</tt> file
/// \details Solution is interpolated on created mesh
    void getSolutionBB(string rBB) { _rBB_file = rBB; _set_rBB = true; }

/** \brief Get the interpolated solution on the new mesh
 *  \details The solution must have been saved on an output bb file
 *  @param [out] u Vector that contains on output the obtained solutions. This
 *                 vector is resized before being initialized
 *  @param [in] is [Default: <tt>1</tt>]
 */
    void getSolution(Vect<real_t>& u,
                     int           is=1);

/// \brief Write the file of interpolation of the solutions in <tt>bb</tt> file
    void getInterpolatedSolutionbb() { _wbb_file = "adapt.wbb"; _set_wbb = true; }

/// \brief Write the file of interpolation of the solutions in <tt>BB</tt> file
    void getInterpolatedSolutionBB() { _wBB_file = "adapt.wBB"; _set_wBB = true; }

/// \brief Set angular limit for a corner (in degrees)
/// \details The angle is defined from 2 normals of 2 consecutive edges
    void setTheta(real_t theta) { _cutoff_radian = theta*OFELI_PI/180.0; }

/// \brief Split triangles into 4 triangles
    void Split() { _allquad = 2; }

/** \brief Save a solution in metric file
 *  @param [in] file File name where the metric is stored
 *  @param [in] u Solution vector to store
 */
    void saveMbb(string              file,
                 const Vect<real_t>& u);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Interpolate a solution vector from a former mesh to the generated one
 *  @param [in] u Solution vector to store
 *  @param [in] v Interpolated solution vector
 *  @param [in] nx Number of x-grid cells (Interpolation uses an intermediate grid)
 *              [Default: 100]
 *  @param [in] ny Number of y-grid cells (Interpolation uses an intermediate grid)
 *              [Default: 100]
 */
/*   void Interpolate(const Vect<real_t>& u,
                          Vect<real_t>& v,
                          size_t        nx=100,
                          size_t        ny=100);*/

   void Interpolate(const Vect<real_t>& u,
                    Vect<real_t>&       v);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Run adaptation process
 *  @return Return code: 
 *  <ul>
 *     <li>= 0: Adaptation has been normally completed
 *     <li>= 1: An error occured
 *  </ul>
 */
    int run();

/** \brief Run adaptation process using a solution vector
 *  @param [in] u Solution vector defined on the input mesh
 *  @return Return code: 
 *  <ul>
 *     <li>= 0: Adaptation has been normally completed
 *     <li>= 1: An error occured
 *  </ul>
 */
    int run(const Vect<real_t>& u);


/** \brief Run adaptation process using a solution vector and interpolates solution on the
 *  adapted mesh
 *  @param [in] u Solution vector defined on the input mesh
 *  @param [in] v Solution vector defined on the (adapted) output mesh
 *  @return Return code: 
 *  <ul>
 *     <li>= 0: Adaptation has been normally completed
 *     <li>= 1: An error occured
 *  </ul>
 */
    int run(const Vect<real_t>& u,
            Vect<real_t>&       v);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  friend ostream & operator<<(ostream& s, const MeshAdapt& a);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 private:

   vector<Mesh *> _ms;
   Mesh *_theMesh;
   const Vect<real_t> *_u;
   int    _nb_Jacobi, _nb_smooth, _allquad;
   size_t _nb_nodes, _nb_elements, _max_nb_vertices, _iter, _nb_subdiv;
   bool   _abs_error, _anisotropic, _scaling, _keep_vertices;
   bool   _set_outm, _set_geo, _hessian, _set_bgm, _set_metric;
   bool   _splitbedge, _set_mbb, _set_mBB, _set_rbb, _set_rBB;
   bool   _set_wbb, _set_wBB, _set_meshr, _set_ometric;
   real_t _scale_fact, _err, _hmin, _hmax, _hmin_aniso, _aniso_max, _max_subdiv, _omega;
   real_t _ratio, _power, _coef, _geo_err, _cutoff_radian, _cut_off;
   string _output_mesh_file, _geo_file, _bb_file, _BB_file, _background_mesh_file;
   string _metric_file, _mbb_file, _mBB_file, _rbb_file, _rBB_file, _wbb_file;
   string _wBB_file, _ometric_file, _meshr_file;
   Domain *_domain;

   void setDefault();
};


//-----------------------------------------------------------------------------
// Associated functions
//-----------------------------------------------------------------------------

/// \fn ostream & operator<<(ostream& s, const MeshAdapt &a)
/// \brief Output MeshAdapt class data
/// \ingroup Mesh
    ostream & operator<<(ostream&         s,
                         const MeshAdapt& a);


/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
