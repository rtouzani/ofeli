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

                 Definition of Class Material to select materials

  ==============================================================================*/

#ifndef __MATERIAL_H
#define __MATERIAL_H

#include <string>
#include <iostream>
#include <iomanip>
using std::ostream;
using std::endl;
using std::cerr;
using std::setw;
using std::string;

#include "OFELI_Config.h"
#include "io/fparser/fparser.h"

extern FunctionParser theParser;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \defgroup Physics Physical properties of media
 *  \brief Physical properties of materials and media
 */

/*! \file Material.h
 *  \brief Definition file for class Material.
 */

/*! \class Material
 *  \ingroup Physics
 *  \brief To treat material data.
 * This class enables reading material data in material data files.
 * It also returns these informations by means of its members.
 */

enum {
   BY_VALUE,
   BY_POSITION,
   BY_FIELD
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct Prop {
   bool    exist;
   int     type;
   real_t  value;
   string  fp_xyzt;
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

class Material
{

 public:

/// \brief Default consructor.
/// \details It initializes the class and searches for the path
/// where are material data files.
    Material();

/// \brief Copy constructor
    Material(const Material& m);

/// \brief Destructor
    ~Material();

/** \brief Associate to material code number <tt>n</tt> the material named <tt>name</tt>
 *  \return Number of materials
 */
    int set(int           m,
            const string& name);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void scanXML();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return material name for material with code <tt>m</tt>
/// \details If such a material is not found, return a blank string.
    string getName(int m) const;

/// \brief Return material code for <tt>i</tt>-th material.
    int getCode(size_t i) const;

/// \brief Return Number of read materials.
    size_t getNbMat() const;

/// \brief Associate code <tt>m</tt> to current material.
    void setCode(int m);

/// Check if material code <tt>c</tt> is present.
/// \return <tt>0</tt> if succeeded, <tt>1</tt> if not.
    int check(int c);

/// \brief Return constant density.
    real_t Density();

/// \brief Return density at point <tt>x</tt> and time <tt>t</tt>
    real_t Density(const Point<real_t>& x,
                   real_t               t);

/// \brief Return constant specific heat.
    real_t SpecificHeat();

/// \brief Return specific heat at point <tt>x</tt> and time <tt>t</tt>
    real_t SpecificHeat(const Point<real_t>& x,
                        real_t               t);

/// \brief Return constant thermal conductivity.
    real_t ThermalConductivity();

/// \brief Return thermal conductivity at point <tt>x</tt> and time <tt>t</tt>
    real_t ThermalConductivity(const Point<real_t>& x,
                               real_t               t);

/// \brief Return constant melting temperature.
    real_t MeltingTemperature();

/// \brief Return melting temperature at point <tt>x</tt> and time <tt>t</tt>
    real_t MeltingTemperature(const Point<real_t>& x,
                              real_t               t);

/// \brief Return constant evaporation temperature.
    real_t EvaporationTemperature();

/// \brief Return evaporation temperature at point <tt>x</tt> and time <tt>t</tt>
    real_t EvaporationTemperature(const Point<real_t>& x,
                                  real_t               t);

/// \brief Return constant thermal expansion coefficient.
    real_t ThermalExpansion();

/// \brief Return thermal expansion coefficient at point <tt>x</tt> and time <tt>t</tt>
    real_t ThermalExpansion(const Point<real_t>& x,
                            real_t               t);

/// \brief Return constant latent heat for melting.
    real_t LatentHeatForMelting();

/// \brief Return latent heat for melting at point <tt>x</tt> and time <tt>t</tt>
    real_t LatentHeatForMelting(const Point<real_t>& x,
                                real_t               t);

/// \brief Return constant latent heat for evaporation.
    real_t LatentHeatForEvaporation();

/// \brief Return latent heat for evaporation at point <tt>x</tt> and time <tt>t</tt>
    real_t LatentHeatForEvaporation(const Point<real_t>& x,
                                    real_t               t);

/// \brief Return constant dielectric constant.
    real_t DielectricConstant();

/// \brief Return dielectric constant at point <tt>x</tt> and time <tt>t</tt>
    real_t DielectricConstant(const Point<real_t>& x,
                              real_t               t);

/// \brief Return constant electric conductivity.
    real_t ElectricConductivity();

/// \brief Return electric conductivity at point <tt>x</tt> and time <tt>t</tt>
    real_t ElectricConductivity(const Point<real_t>& x,
                                real_t               t);

/// \brief Return constant electric resistivity.
    real_t ElectricResistivity();

/// \brief Return electric resistivity at point <tt>x</tt> and time <tt>t</tt>
    real_t ElectricResistivity(const Point<real_t>& x,
                               real_t               t);

/// \brief Return constant magnetic permeability
    real_t MagneticPermeability();

/// \brief Return magnetic permeability at point <tt>x</tt> and time <tt>t</tt>
    real_t MagneticPermeability(const Point<real_t>& x,
                                real_t               t);

/// \brief Return constant viscosity
    real_t Viscosity();

/// \brief Return viscosity at point <tt>x</tt> and time <tt>t</tt>
    real_t Viscosity(const Point<real_t>& x,
                     real_t               t);

/// \brief Return constant Young modulus
    real_t YoungModulus();

/// \brief Return Young modulus at point <tt>x</tt> and time <tt>t</tt>
    real_t YoungModulus(const Point<real_t>& x,
                        real_t               t);

/// \brief Return constant Poisson ratio
    real_t PoissonRatio();

/// \brief Return Poisson ratio at point <tt>x</tt> and time <tt>t</tt>
    real_t PoissonRatio(const Point<real_t>& x,
                        real_t               t);

/// \brief Return constant <tt>i</tt>-th property
    real_t Property(int i);

/// \brief Return <tt>i</tt>-th property at point <tt>x</tt> and time <tt>t</tt>
    real_t Property(int                  i,
                    const Point<real_t>& x,
                    real_t               t);

/// \brief Operator =
    Material & operator=(const Material& m)
    {
       _nb_mat = m._nb_mat;
       _path = m._path;
       for (size_t i=0; i<_nb_mat; i++) {
          _code[i] = m._code[i];
          _mat[i] = m._mat[i];
       }
       //       _ff = m._ff;
       set_kw();
       return *this;
    }

    friend class XMLParser;

 private:

   int _state, _type;
   size_t _nb_mat, _index_mat;
   string _mat[MAX_NB_MATERIALS], _name, _path;
   int _code[MAX_NB_MATERIALS];
   void getValue(size_t i);
   void getFunction(size_t i);
   void getTable() { }
   real_t _data[10];

   void set_kw();
   real_t getProperty(Prop& prop);
   real_t getProperty(Prop& prop, const Point<real_t>& x, const real_t& t);
   Prop _density[MAX_NB_MATERIALS];
   Prop _specific_heat[MAX_NB_MATERIALS];
   Prop _thermal_conductivity[MAX_NB_MATERIALS];
   Prop _melting_temperature[MAX_NB_MATERIALS];
   Prop _evaporation_temperature[MAX_NB_MATERIALS];
   Prop _thermal_expansion[MAX_NB_MATERIALS];
   Prop _latent_heat_melting[MAX_NB_MATERIALS];
   Prop _latent_heat_evaporation[MAX_NB_MATERIALS];
   Prop _dielectric_constant[MAX_NB_MATERIALS];
   Prop _electric_conductivity[MAX_NB_MATERIALS];
   Prop _electric_resistivity[MAX_NB_MATERIALS];
   Prop _magnetic_permeability[MAX_NB_MATERIALS];
   Prop _viscosity[MAX_NB_MATERIALS];
   Prop _young_modulus[MAX_NB_MATERIALS];
   Prop _poisson_ratio[MAX_NB_MATERIALS];
   Prop _prop[MAX_NB_MATERIALS][10];
};


/** \fn ostream & operator<<(ostream& s, const Material& m)
 *  \brief Output material data.
 *  \ingroup Mesh
 */
    ostream& operator<<(      ostream&  s,
                        const Material& m);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
