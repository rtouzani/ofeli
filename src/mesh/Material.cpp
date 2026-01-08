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

                         Implementation of class 'Material'

  ==============================================================================*/

#include <iostream>
#include <fstream>
#include "mesh/Material.h"
#include "io/XMLParser.h"
#include "linear_algebra/Point.h"
#include "OFELIException.h"

using std::cout;
using std::ifstream;

namespace OFELI {

Material::Material()
{
   for (size_t i=0; i<MAX_NB_MATERIALS; i++)
      _code[i] = MY_RANDOM;
   _nb_mat = 0;
   if (getenv("OFELI_PATH_MATERIAL")==nullptr)
      _path = PATH_MATERIAL;
   else
      _path = getenv("OFELI_PATH_MATERIAL");
   string mat_file = _path + PATH_SEP;
   mat_file += "Generic.md";
   if (!ifstream(mat_file.c_str())) {
      cout << "Error in class Material: Material directory not found." << endl;
      cout << "Edit the file include/OFELI_Config.h and modify the line with PATH_MATERIAL," << endl;
      cout << "or Set the pathname for the directory Material in the environment variable OFELI_PATH_MATERIAL." << endl;
      exit(1);
   }
}


Material::Material(const Material& m)
{
   _nb_mat = m._nb_mat;
   _path = m._path;
   for (size_t i=0; i<_nb_mat; i++) {
      _code[i] = m._code[i];
      _mat[i] = m._mat[i];
   }
}


Material::~Material() { }


void Material::scanXML()
{
   for (size_t i=0; i<_nb_mat; i++) {
      string file = _path + PATH_SEP;
      file += _mat[i] + MATERIAL_EXT;
      XMLParser p(file,EType::MATERIAL);
      p.setMaterialNumber(i);
      p.getMaterial();
   }
}


int Material::set(int           m,
                  const string& name)
{
   string file = _path + PATH_SEP;
   string mf = file + "Generic.md";
   ifstream im(mf.c_str());
   if (im.fail()) {
      cout << "Error: File " << mf << " not found." << endl;
      cout << "Your OFELI installation settings show that material files must be in the directory: " << file << endl;
      exit(1);
   }
   if (_code[_nb_mat] != m) {
      _code[_nb_mat] = m;
      _mat[_nb_mat] = name;
      file += _mat[_nb_mat] + MATERIAL_EXT;
      {
         if (ifstream(file.c_str()).fail())
            file = _mat[_nb_mat] + MATERIAL_EXT;
      }
      XMLParser p(file,EType::MATERIAL);
      p.setMaterialNumber(_nb_mat);
      p.getMaterial();
   }
   _nb_mat++;
   return _nb_mat;
}


int Material::check(int c)
{
   for (size_t i=0; i<_nb_mat; i++)
      if (_code[i] == c)
         return 0;
   set(c,"Generic");
   return 1;
}


string Material::getName(int m) const
{
   for (size_t i=0; i<_nb_mat; i++)
      if (_code[i] == m)
         return _mat[i];
   return " ";
}


size_t Material::getNbMat() const
{
   return _nb_mat;
}


void Material::setCode(int m)
{
   _index_mat = 0;
   for (size_t i=0; i<_nb_mat; i++) {
      if (_code[i] == m)
         _index_mat = i;
   }
}


int Material::getCode(size_t i) const
{
   return _code[i-1];
}


real_t Material::getProperty(Prop& prop)
{
   real_t val;
   if (prop.type==BY_VALUE)
      val = prop.value;
   else {
      _theFct.set(prop.fp_xyzt);
      val = _theFct(0.,0.,0.,0.);
   }
   return val;
}


real_t Material::getProperty(Prop&                prop,
                             const Point<real_t>& x,
                             const real_t&        t)
{
   real_t val;
   if (prop.type==BY_VALUE)
      val = prop.value;
   else {
      _theFct.set(prop.fp_xyzt);
      val = _theFct(x,t);
   }
   return val;
}


real_t Material::Density()
{
   if (_density[_index_mat].exist)
      return getProperty(_density[_index_mat]);
   else
      throw OFELIException("Material::Density(): This property is not present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::Density(const Point<real_t>& x,
                         real_t               t)
{
   if (_density[_index_mat].exist)
      return getProperty(_density[_index_mat],x,t);
   else
      throw OFELIException("Material::Density(Point<real_t>,real_t): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::SpecificHeat()
{
   if (_specific_heat[_index_mat].exist)
      return getProperty(_specific_heat[_index_mat]);
   else
      throw OFELIException("Material::SpecificHeat(): This property is not present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::SpecificHeat(const Point<real_t>& x,
                              real_t               t)
{
   if (_specific_heat[_index_mat].exist)
      return getProperty(_specific_heat[_index_mat],x,t);
   else
      throw OFELIException("Material::SpecificHeat(Point<real_t>,real_t): This property is not present for material " +
   _mat[_index_mat]);
   return 0.;
}


real_t Material::ThermalConductivity()
{
   if (_thermal_conductivity[_index_mat].exist)
      return getProperty(_thermal_conductivity[_index_mat]);
   else
      throw OFELIException("Material::ThermalConductivity(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::ThermalConductivity(const Point<real_t>& x,
                                     real_t               t)
{
   if (_thermal_conductivity[_index_mat].exist)
      return getProperty(_thermal_conductivity[_index_mat],x,t);
   else
      throw OFELIException("Material::ThermalConductivity(Point<real_t>,real_t): This "
                           "property is not present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::MeltingTemperature()
{
   if (_melting_temperature[_index_mat].exist)
      return getProperty(_melting_temperature[_index_mat]);
   else
      throw OFELIException("Material::MeltingTemperature(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::MeltingTemperature(const Point<real_t>& x,
                                    real_t               t)
{
   if (_melting_temperature[_index_mat].exist)
      return getProperty(_melting_temperature[_index_mat],x,t);
   else
      throw OFELIException("Material::MeltingTemperature(Point<real_t>,real_t): This property "
                           "is not present for material "+_mat[_index_mat]);
   return 0.;
}


real_t Material::EvaporationTemperature()
{
   if (_evaporation_temperature[_index_mat].exist)
      return getProperty(_evaporation_temperature[_index_mat]);
   else
      throw OFELIException("Material::EvaporationTemperature(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::EvaporationTemperature(const Point<real_t>& x,
                                        real_t               t)
{
   if (_evaporation_temperature[_index_mat].exist)
      return getProperty(_evaporation_temperature[_index_mat],x,t);
   else
      throw OFELIException("Material::EvaporationTemperature(Point<real_t>,real_t): This property is not "
                           "present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::ThermalExpansion()
{
   if (_thermal_expansion[_index_mat].exist)
      return getProperty(_thermal_expansion[_index_mat]);
   else
      throw OFELIException("Material::ThermalExpansion(): This property is not present for material "+_mat[_index_mat]);
   return 0.;
}


real_t Material::ThermalExpansion(const Point<real_t>& x,
                                  real_t               t)
{
   if (_thermal_expansion[_index_mat].exist)
      return getProperty(_thermal_expansion[_index_mat],x,t);
   else
      throw OFELIException("Material::ThermalExpansion(Point<real_t>,real_t): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::LatentHeatForMelting()
{
   if (_latent_heat_melting[_index_mat].exist)
      return getProperty(_latent_heat_melting[_index_mat]);
   else
      throw OFELIException("Material::LatentHeatForMelting(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::LatentHeatForMelting(const Point<real_t>& x,
                                      real_t              t)
{
   if (_latent_heat_melting[_index_mat].exist)
      return getProperty(_latent_heat_melting[_index_mat],x,t);
   else
      throw OFELIException("Material::LatentHeatForMelting(Point<real_t>,real_t): This property is not present for material "
                           + _mat[_index_mat]);
   return 0.;
}


real_t Material::LatentHeatForEvaporation()
{
   if (_latent_heat_evaporation[_index_mat].exist)
      return getProperty(_latent_heat_evaporation[_index_mat]);
   else
      throw OFELIException("Material::LatentHeatForEvaporation(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::LatentHeatForEvaporation(const Point<real_t>& x,
                                          real_t               t)
{
   if (_latent_heat_evaporation[_index_mat].exist)
       return getProperty(_latent_heat_evaporation[_index_mat],x,t);
   else
      throw OFELIException("Material::LatentHeatForEvaporation(Point<real_t>,real_t): This property is "
                           "not present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::DielectricConstant()
{
   if (_dielectric_constant[_index_mat].exist)
      return getProperty(_dielectric_constant[_index_mat]);
   else
      throw OFELIException("Material::DielectricConstant(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::DielectricConstant(const Point<real_t>& x,
                                    real_t               t)
{
   if (_dielectric_constant[_index_mat].exist)
      return getProperty(_dielectric_constant[_index_mat],x,t);
   else
      throw OFELIException("Material::DielectricConstant(Point<real_t>,real_t): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::ElectricConductivity()
{
   if (_electric_conductivity[_index_mat].exist)
      return getProperty(_electric_conductivity[_index_mat]);
   else
      throw OFELIException("Material::ElectricConductivity(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::ElectricConductivity(const Point<real_t>& x,
                                      real_t               t)
{
   if (_electric_conductivity[_index_mat].exist)
      return getProperty(_electric_conductivity[_index_mat],x,t);
   else
      throw OFELIException("Material::ElectricConductivity(Point<real_t>,real_t): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::ElectricResistivity()
{
   if (_electric_resistivity[_index_mat].exist)
      return getProperty(_electric_resistivity[_index_mat]);
   else
      throw OFELIException("Material::ElectricResistivity(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::ElectricResistivity(const Point<real_t>& x,
                                     real_t               t)
{
   if (_electric_resistivity[_index_mat].exist)
      return getProperty(_electric_resistivity[_index_mat],x,t);
   else
      throw OFELIException("Material::ElectricResistivity(Point<real_t>,real_t): This property "
                           "is not present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::MagneticPermeability()
{
   if (_magnetic_permeability[_index_mat].exist)
      return getProperty(_magnetic_permeability[_index_mat]);
   else
      throw OFELIException("Material::Density(Point<real_t>,real_t): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::MagneticPermeability(const Point<real_t>& x,
                                      real_t               t)
{
   if (_magnetic_permeability[_index_mat].exist)
      return getProperty(_magnetic_permeability[_index_mat],x,t);
   else
      throw OFELIException("Material::MagneticPermeability(Point<real_t>,real_t): This property "
                           "is not present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::Viscosity()
{
   if (_viscosity[_index_mat].exist)
      return getProperty(_viscosity[_index_mat]);
   else
      throw OFELIException("Material::Viscosity(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::Viscosity(const Point<real_t>& x,
                           real_t               t)
{
   if (_viscosity[_index_mat].exist)
      return getProperty(_viscosity[_index_mat],x,t);
   else
      throw OFELIException("Material::Viscosity(Point<real_t>,real_t): This property is not "
                           "present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::YoungModulus()
{
   if (_young_modulus[_index_mat].exist)
      return getProperty(_young_modulus[_index_mat]);
   else
      throw OFELIException("Material::YoungModulus(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::YoungModulus(const Point<real_t>& x,
                              real_t               t)
{
   if (_young_modulus[_index_mat].exist)
      return getProperty(_young_modulus[_index_mat],x,t);
   else
      throw OFELIException("Material::YoungModulus(Point<real_t>,real_t): This property is not "
                           "present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::PoissonRatio()
{
   if (_poisson_ratio[_index_mat].exist)
      return getProperty(_poisson_ratio[_index_mat]);
   else
      throw OFELIException("Material::PoissonRatio(): This property is not present for material " +
                           _mat[_index_mat]);
   return 0.;
}


real_t Material::PoissonRatio(const Point<real_t>& x,
                              real_t               t)
{
   if (_poisson_ratio[_index_mat].exist)
      return getProperty(_poisson_ratio[_index_mat],x,t);
   else
      throw OFELIException("Material::PoissonRatio(Point<real_t>,real_t): This property is not "
                           "present for material " + _mat[_index_mat]);
   return 0.;
}


real_t Material::Property(int i)
{
   return getProperty(_prop[_index_mat][i-1]);
}


real_t Material::Property(int                  i,
                          const Point<real_t>& x,
                          real_t               t)
{
   return getProperty(_prop[_index_mat][i-1],x,t);
}


void Material::getValue(size_t i)
{
  /*   real_t val = _ff.getD();

   switch (_type) {

      case 0:
         _density[i].value = val;
         _density[i].type = BY_VALUE;
         break;

      case 1:
         _specific_heat[i].value = val;
         _specific_heat[i].type = BY_VALUE;
         break;

      case 2:
         _thermal_conductivity[i].value = val;
         _thermal_conductivity[i].type = BY_VALUE;
         break;

      case 3:
         _melting_temperature[i].value = val;
         _melting_temperature[i].type = BY_VALUE;
         break;

      case 4:
         _evaporation_temperature[i].value = val;
         _evaporation_temperature[i].type = BY_VALUE;
         break;

      case 5:
         _thermal_expansion[i].value = val;
         _thermal_expansion[i].type = BY_VALUE;
         break;

      case 6:
         _latent_heat_melting[i].value = val;
         _latent_heat_melting[i].type = BY_VALUE;
         break;

      case 7:
         _latent_heat_evaporation[i].value = val;
         _latent_heat_evaporation[i].type = BY_VALUE;
         break;

      case 8:
         _dielectric_constant[i].value = val;
         _dielectric_constant[i].type = BY_VALUE;
         break;

      case 9:
         _electric_conductivity[i].value = val;
         _electric_conductivity[i].type = BY_VALUE;
         break;

      case 10:
         _electric_resistivity[i].value = val;
         _electric_resistivity[i].type = BY_VALUE;
         break;

      case 11:
         _magnetic_permeability[i].value = val;
         _magnetic_permeability[i].type = BY_VALUE;
         break;

      case 12:
         _viscosity[i].value = val;
         _viscosity[i].type = BY_VALUE;
         break;

      case 13:
         _young_modulus[i].value = val;
         _young_modulus[i].type = BY_VALUE;
         break;

      case 14:
         _poisson_ratio[i].value = val;
         _poisson_ratio[i].type = BY_VALUE;
         break;
   }
}


void Material::getFunction(size_t i)
{
   int _typef;
   string tf = _ff.getS();
   if (tf == "xyzt")
     _typef = 0;
   else
     _typef = 1;
   string str = _ff.getE();

   switch (_type) {

      case 0:
         if (_typef==0) {
            _density[i].fp_xyzt = str;
            _density[i].type = BY_POSITION;
         }
         break;

      case 1:
         if (_typef==0) {
            _specific_heat[i].fp_xyzt = str;
            _specific_heat[i].type = BY_POSITION;
         }
         break;

      case 2:
         if (_typef==0) {
            _thermal_conductivity[i].fp_xyzt = str;
            _thermal_conductivity[i].type = BY_POSITION;
         }
         break;

      case 3:
         if (_typef==0) {
            _melting_temperature[i].fp_xyzt = str;
            _melting_temperature[i].type = BY_POSITION;
         }
         break;

      case 4:
         if (_typef==0) {
            _evaporation_temperature[i].fp_xyzt = str;
            _evaporation_temperature[i].type = BY_POSITION;
         }
         break;

      case 5:
         if (_typef==0) {
            _thermal_expansion[i].fp_xyzt = str;
            _thermal_expansion[i].type = BY_POSITION;
         }
         break;

      case 6:
         if (_typef==0) {
            _latent_heat_melting[i].fp_xyzt = str;
            _latent_heat_melting[i].type = BY_POSITION;
         }
         break;

      case 7:
         if (_typef==0) {
            _latent_heat_evaporation[i].fp_xyzt = str;
            _latent_heat_evaporation[i].type = BY_POSITION;
         }
         break;

      case 8:
         if (_typef==0) {
            _dielectric_constant[i].fp_xyzt = str;
            _dielectric_constant[i].type = BY_POSITION;
         }
         break;

      case 9:
         if (_typef==0) {
            _electric_conductivity[i].fp_xyzt = str;
            _electric_conductivity[i].type = BY_POSITION;
         }
         break;

      case 10:
         if (_typef==0) {
            _electric_resistivity[i].fp_xyzt = str;
            _electric_resistivity[i].type = BY_POSITION;
         }
         break;

      case 11:
         if (_typef==0) {
            _magnetic_permeability[i].fp_xyzt = str;
            _magnetic_permeability[i].type = BY_POSITION;
         }
         break;

      case 12:
         if (_typef==0) {
            _viscosity[i].fp_xyzt = str;
            _viscosity[i].type = BY_POSITION;
         }
         break;

      case 13:
         if (_typef==0) {
            _young_modulus[i].fp_xyzt = str;
            _young_modulus[i].type = BY_POSITION;
         }
         break;

      case 14:
         if (_typef==0) {
            _poisson_ratio[i].fp_xyzt = str;
            _poisson_ratio[i].type = BY_POSITION;
         }
         break;
   }*/
}


ostream& operator<<(ostream&        s,
                    const Material& m)
{
   s << "\nNumber of Materials: " << m.getNbMat() << endl;
   s << "  Code   Name" << endl;
   for (size_t i=1; i<=m.getNbMat(); i++)
     s << setw(6) << m.getCode(i) << "   " << m.getName(m.getCode(i)) << endl;
   s << endl;
   return s;
}

Material theMaterial;

} /* namespace OFELI */
