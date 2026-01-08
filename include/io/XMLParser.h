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

           Definition of class 'XMLParser' for parsing OFELI XML files

  ==============================================================================*/


#ifndef __XML_PARSER_H
#define __XML_PARSER_H

#include <map>
#include <fstream>
#include <string>
#include <vector>
using std::string;
using std::vector;
using std::map;

#include "OFELI_Config.h"
#include "io/xmlsp/xmlsp.h"
#include "io/Prescription.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/DMatrix.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

typedef std::pair<string,string> SString;
#define IPF_TOKEN(s)  if (it!=tokens.end()) s = *it; break;
#define IPF_TOKENI(s) if (it!=tokens.end()) s = stoi(*it); break;
#define IPF_TOKENF(s) if (it!=tokens.end()) s = stof(*it); break;
#define IPF_ISET(s)   if (!_scan) { if (a.first=="value") s = stoi(a.second); } break;
#define IPF_DSET(s)   if (!_scan) { if (a.first=="value") s = stof(a.second); } break;
#define IPF_SSET(s)   if (!_scan) { if (a.first=="value") s = a.second; } break;


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

class Mesh;
class Grid;
class IPF;
class Domain;
class Tabulation;

class XMLParser : public Parser
{

 public:

   enum {ASCII,BINARY};

   enum class Tag {PROJECT, TITLE, AUTHOR, VERBOSE, OUTPUT, SAVE, PLOT, DOMAIN_FILE, MESH_FILE, INIT_FILE,
                   RESTART_FILE, BC, BC_FILE, BF, BF_FILE, SF, SF_FILE, SAVE_FILE, PLOT_FILE, PRESCRIPTION_FILE,
                   AUX_FILE, DENSITY, SPECIFIC_HEAT, THERMAL_CONDUCTIVITY, MELTING_TEMPERATURE, EVAPORATION_TEMPERATURE,
                   THERMAL_EXPANSION, LATENT_HEAT_MELTING, LATENT_HEAT_EVAPORATION, DIELECTRIC_CONSTANT, 
                   ELECTRIC_CONDUCTIVITY, ELECTRIC_RESISTIVITY, ELECTRIC_PERMITTIVITY, MAGNETIC_PERMEABILITY,
                   POISSON_RATIO, RHO_CP, VISCOSITY, YOUNG_MODULUS, VECTOR, CODE, NB_DOF, ARRAY, DDOMAIN,
                   VERTEX, LINE, CIRCLE, CONTOUR, HOLE, REQUIRED_VERTEX, REQUIRED_EDGE, RECTANGLE, SUBDOMAIN,
                   MESH, NODES, ELEMENTS, SIDES, EDGES, TIME, STEP, NB_STEPS, TIME_STEP,
                   MAX_TIME, NB_ITER, TOLERANCE, LABEL, PARAMETER, INTEGER, DOUBLE, STRING, COMPLEX, VALUE, DOF,
                   GRID, CONSTANT, EXPRESSION, MATERIAL, MATERIAL_CODE, DATA, MATRIX, NAME,
                   PRESCRIPTION, SOLUTION, BOUNDARY_CONDITION, BODY_FORCE, POINT_FORCE, FLUX, INITIAL, FUNCTION, VARIABLE};

   XMLParser();
   XMLParser(string file, EType type=EType::MESH);
   XMLParser(string file, Mesh& ms, EType type=EType::MESH, int format=ASCII);
   XMLParser(string file, Grid& g, int format=ASCII);
   XMLParser(string file, Vect<real_t>& v, int format=ASCII);
   XMLParser(const XMLParser &p);
   ~XMLParser();
   void setMaterialNumber(int m);
   int scan(size_t k=1);
   int scan(vector<real_t>& t, int type);
   void setFile(string file);
   void set(Mesh& ms, int format=ASCII);
   void set(Grid& gr);
   int get(EType type, vector<Prescription::PPar>& p);
   int get(Mesh& ms, int format=ASCII);
   int get(Mesh& ms, vector<real_t>& v, string& name);
   int get(Mesh& ms, Vect<real_t>& v, real_t time=-1, string name="ANYTHING", int format=ASCII);
   int get(Grid& gr, vector<real_t>& v, string& name);
   int get(Grid& gr, Vect<real_t>& v, real_t time=-1, string name="ANYTHING", int format=ASCII);
   int get(Vect<real_t>& v, real_t time, string name="ANYTHING", int format=ASCII);
   int get(Vect<real_t>& v, const string& name="ANYTHING");
   int get(IPF& ipf);
   int get(Domain& dm);
   int get(Tabulation& t);
   int get(Matrix<real_t>* A);
   int get(DMatrix<real_t>& A);
   int getArray(Vect<real_t>& A);
   int getMaterial();
   int get(int type, Vect<real_t>& v, real_t time=0);
   size_t getNbDOF() const { return _nb_dof; }
   MatrixSize MSize() const;

 protected:
   std::ifstream _is;
   bool _is_opened, _is_closed, _set_mesh, _set_grid, _set_vector, _set_file;
   bool _set_prescription, _set_matrix, _set_domain, _prescription_opened, _compact, _value;
   real_t _x, _y, _z, _time, _sought_time, _val;
   int _access, _cm, _format, _var, _code, _rtype, _scan;
   EType _type, _prescription_type, _old_type;
   Prescription::PPar _par;
   string _file, _mesh_file, _vect_file, _el_shape, _sd_shape, _name, _sought_name, _tag_name, _xml, _mat, _slabel, _svalue;
   size_t _iter, _dof, _label, _nb_dof, _dim, _nb_nodes, _nb_elements, _nb_sides, _nb_edges, _tab_size, _vect_size;
   size_t _nb_el_nodes, _nb_sd_nodes, _nb_mat, _all_steps, _nb_funct;
   size_t _nb_var, _nx, _ny, _nz, _nt;
   MatrixType _storage;
   MatrixSize _msize;
   DOFSupport _dof_support;
   mutable Mesh *_theMesh;
   Grid *_theGrid;
   Vect<real_t> *_v, *_theArray;
   XMLParser *_parser;
   IPF *_ipf;
   Domain *_theDomain;
   Tabulation *_theTabulation;
   Matrix<real_t> *_theMatrix;
   vector<Prescription::PPar> *_vp;
   vector<real_t> *_ft, *_V;
   Fct _theFct;

   virtual bool on_tag_open(string tag_name, StringMap& attributes);
   virtual bool on_cdata(string cdata);
   virtual bool on_tag_close(string tag_name);
   virtual bool on_comment(string comment);
   virtual bool on_processing(string value);
   virtual bool on_doctype(string value);
   virtual bool on_document_begin();
   virtual bool on_document_end();
   virtual void on_error(int e, int line, int col, string message);
   void open();
   void read_project(const SString& a);
   void read_project_data(const vector<string> &tokens, std::vector<string>::iterator &it);
   void read_parameter_data(const vector<string> &tokens, std::vector<string>::iterator &it);
   void read_domain(const SString& a);
   void read_domain_data(const vector<string> &tokens, std::vector<string>::iterator &it);
   void read_mat(const SString& a);
   void read_mat_data(const vector<string> &tokens, vector<string>::iterator &it);
   void read_mesh(const SString& a);
   void read_mesh_data(const vector<string> &tokens, vector<string>::iterator &it);
   void read_vect(const SString& a);
   void read_vect_data(const vector<string> &tokens, vector<string>::iterator &it);
   void read_const_vect_data(const vector<string> &tokens, vector<string>::iterator &it);
   void read_const_vect_data();
   void read_exp_vect_data(const vector<string> &tokens, vector<string>::iterator &it);
   void parse_exp(size_t n, size_t k);
   void read_prescription(const SString& a);
   void read_prescribe_data(const vector<string> &tokens, vector<string>::iterator &it);
   void read_tab(const SString& a);
   void read_tab_data(const vector<string> &tokens, vector<string>::iterator &it);
   int read_matrix(const SString& a);
   int read_matrix_data(const vector<string> &tokens, vector<string>::iterator &it);
   MatrixSize getMatrixData();

   map<string,MatrixType> _STORAGE = {{"dense",DENSE},
                                      {"skyline",SKYLINE},
                                      {"sparse",SPARSE},
                                      {"diagonal",DIAGONAL},
                                      {"tridiagonal",TRIDIAGONAL},
                                      {"band",BAND},
                                      {"symmetric",SYMMETRIC},
                                      {"unsymmetric",UNSYMMETRIC},
                                      {"identity",IDENTITY}
                                     };

   map<string,Tag> _TS = {{"Project",Tag::PROJECT},
                          {"title",Tag::TITLE},
                          {"author",Tag::AUTHOR},
                          {"verbose",Tag::VERBOSE},
                          {"output",Tag::OUTPUT},
                          {"domain_file",Tag::DOMAIN_FILE},
                          {"mesh_file",Tag::MESH_FILE},
                          {"init_file",Tag::INIT_FILE},
                          {"restart_file",Tag::RESTART_FILE},
                          {"bc",Tag::BC},
                          {"bc_file",Tag::BC_FILE},
                          {"bf",Tag::BF},
                          {"bf_file",Tag::BF_FILE},
                          {"sf",Tag::SF},
                          {"sf_file",Tag::SF_FILE},
                          {"save",Tag::SAVE},
                          {"save_file",Tag::SAVE_FILE},
                          {"plot",Tag::PLOT},
                          {"plot_file",Tag::PLOT_FILE},
                          {"prescription_file",Tag::PRESCRIPTION_FILE},
                          {"aux_file",Tag::AUX_FILE},
                          {"Density",Tag::DENSITY},
                          {"SpecificHeat",Tag::SPECIFIC_HEAT},
                          {"ThermalConductivity",Tag::THERMAL_CONDUCTIVITY},
                          {"MeltingTemperature",Tag::MELTING_TEMPERATURE},
                          {"EvaporationTemperature",Tag::EVAPORATION_TEMPERATURE},
                          {"ThermalExpansion",Tag::THERMAL_EXPANSION},
                          {"LatentHeatMelting",Tag::LATENT_HEAT_MELTING},
                          {"LatentHeatEvaporation",Tag::LATENT_HEAT_EVAPORATION},
                          {"DielectricConstant",Tag::DIELECTRIC_CONSTANT},
                          {"ElectricConductivity",Tag::ELECTRIC_CONDUCTIVITY},
                          {"ElectricResistivity",Tag::ELECTRIC_RESISTIVITY},
                          {"ElectricPermittivity",Tag::ELECTRIC_PERMITTIVITY},
                          {"MagneticPermeability",Tag::MAGNETIC_PERMEABILITY},
                          {"PoissonRatio",Tag::POISSON_RATIO},
                          {"RhoCp",Tag::RHO_CP},
                          {"Viscosity",Tag::VISCOSITY},
                          {"YoungModulus",Tag::YOUNG_MODULUS},
                          {"Vector",Tag::VECTOR},
                          {"Field",Tag::VECTOR},
                          {"code",Tag::CODE},
                          {"nb_dof",Tag::NB_DOF},
                          {"array",Tag::ARRAY},
                          {"Domain",Tag::DDOMAIN},
                          {"vertex",Tag::VERTEX},
                          {"line",Tag::LINE},
                          {"circle",Tag::CIRCLE},
                          {"contour",Tag::CONTOUR},
                          {"hole",Tag::HOLE},
                          {"RequiredVertex",Tag::REQUIRED_VERTEX},
                          {"RequiredEdge",Tag::REQUIRED_EDGE},
                          {"rectangle",Tag::RECTANGLE},
                          {"subdomain",Tag::SUBDOMAIN},
                          {"Mesh",Tag::MESH},
                          {"Nodes",Tag::NODES}, {"nodes",Tag::NODES},
                          {"Elements",Tag::ELEMENTS}, {"elements",Tag::ELEMENTS},
                          {"Sides",Tag::SIDES}, {"sides",Tag::SIDES},
                          {"Edges",Tag::EDGES}, {"edges",Tag::EDGES},
                          {"time",Tag::TIME},
                          {"step",Tag::STEP},
                          {"nb_steps",Tag::NB_STEPS},
                          {"time_step",Tag::TIME_STEP},
                          {"max_time",Tag::MAX_TIME},
                          {"nb_iter",Tag::NB_ITER},
                          {"tolerance",Tag::TOLERANCE},
                          {"parameter",Tag::PARAMETER},
                          {"int",Tag::INTEGER},
                          {"double",Tag::DOUBLE},
                          {"string",Tag::STRING},
                          {"complex",Tag::COMPLEX},
                          {"value",Tag::VALUE},
                          {"dof",Tag::DOF},
                          {"Grid",Tag::GRID},
                          {"constant",Tag::CONSTANT},
                          {"expression",Tag::EXPRESSION},
                          {"Material",Tag::MATERIAL},
                          {"material_code",Tag::MATERIAL_CODE},
                          {"Data",Tag::DATA},
                          {"Matrix",Tag::MATRIX},
                          {"name",Tag::NAME},
                          {"Prescription",Tag::PRESCRIPTION},
                          {"Solution",Tag::SOLUTION},
                          {"BoundaryCondition",Tag::BOUNDARY_CONDITION},
                          {"BodyForce",Tag::BODY_FORCE}, {"Source",Tag::BODY_FORCE},
                          {"PointForce",Tag::POINT_FORCE},
                          {"Flux",Tag::FLUX}, {"Traction",Tag::FLUX}, {"BoundaryForce",Tag::FLUX},
                          {"Initial",Tag::INITIAL},
                          {"Function",Tag::FUNCTION},
                          {"variable",Tag::VARIABLE}, {"Variable",Tag::VARIABLE}
                         };
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif
