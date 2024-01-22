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

                      Implementation of class 'XMLParser'

  ==============================================================================*/

#include "io/XMLParser.h"
#include "util/macros.h"
#include "mesh/Domain.h"
#include "mesh/Mesh.h"
#include "mesh/MeshUtil.h"
#include "mesh/Grid.h"
#include "io/IPF.h"
#include "mesh/Material.h"
#include "shape_functions/Line2.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Quad4.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Hexa8.h"
#include "shape_functions/Penta6.h"
#include "io/Tabulation.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/DMatrix_impl.h"
#include "OFELIException.h"

namespace OFELI {

extern Material theMaterial;

XMLParser::XMLParser()
          : _is_opened(false), _set_mesh(false), _set_grid(false), _set_vector(false),
            _set_file(false), _set_prescription(false), _set_matrix(false), _set_domain(false),
            _prescription_opened(false), _scan(0), _nb_dof(1), _dim(2), _nb_nodes(0),
            _nb_elements(0), _nb_sides(0), _nb_edges(0), _nb_mat(0),
            _dof_support(NODE_DOF), _theMesh(nullptr), _theGrid(nullptr), _v(nullptr),
            _parser(nullptr), _ipf(nullptr)
{ }


XMLParser::XMLParser(string file,
                     EType  type)
          : _is_opened(false), _set_mesh(true), _set_grid(false), _set_vector(false), _set_file(true),
            _set_prescription(false), _set_matrix(false), _set_domain(false), _prescription_opened(false),
            _scan(0), _file(file), _nb_dof(1), _dim(2), _nb_nodes(0), _nb_elements(0), _nb_sides(0),
            _nb_edges(0), _nb_mat(0), _dof_support(NODE_DOF), _theMesh(nullptr), _theGrid(nullptr),
            _v(nullptr), _parser(nullptr), _ipf(nullptr)
{
   open();
   _type = _old_type = type;
   _x = _y = _z = _time = 0.;
}


XMLParser::XMLParser(string file,
                     Mesh&  ms,
                     EType  type,
                     int    format)
          : _is_opened(false), _set_mesh(true), _set_grid(false), _set_vector(false), _set_file(true),
            _set_prescription(false), _set_matrix(false), _set_domain(false), _prescription_opened(false),
            _format(format), _file(file), _nb_dof(1), _nb_mat(0), _dof_support(NODE_DOF),
            _theMesh(&ms), _theGrid(nullptr), _v(nullptr), _parser(nullptr), _ipf(nullptr)
{
   _type = _old_type = type;
   _scan = 0;
   _nb_nodes = _theMesh->getNbNodes();
   _nb_elements = _theMesh->getNbElements();
   _nb_sides = _theMesh->getNbSides();
   _nb_edges = _theMesh->getNbEdges();
   _dim = _theMesh->getDim();
   open();
   if (type==EType::MESH)
      get(ms);
}


XMLParser::XMLParser(string file,
                     Grid&  g,
                     int    format)
          : _is_opened(false), _set_mesh(false), _set_grid(true), _set_vector(false), _set_file(true),
            _set_prescription(false), _set_matrix(false), _set_domain(false), _prescription_opened(false),
            _format(format), _file(file), _nb_dof(1), _nb_mat(0), _dof_support(NODE_DOF),
            _theMesh(nullptr), _theGrid(&g), _v(nullptr), _parser(nullptr), _ipf(nullptr)
{
   _scan = 0;
   _type = _old_type = EType::GRID;
   _nb_nodes = _theMesh->getNbNodes();
   _nb_elements = _theMesh->getNbElements();
   _nb_sides = _theMesh->getNbSides();
   _nb_edges = _theMesh->getNbEdges();
   _dim = _theMesh->getDim();
   open();
}


XMLParser::XMLParser(const XMLParser& p)
          : _is_opened(p._is_opened), _is_closed(p._is_closed), _set_mesh(p._set_mesh), _set_grid(p._set_grid), 
            _set_vector(p._set_vector), _set_file(p._set_file), _set_prescription(p._set_prescription),
            _set_matrix(p._set_matrix), _set_domain(p._set_domain),
            _time(p._time), _sought_time(p._sought_time),
            _format(p._format), _file(p._file), _mesh_file(p._mesh_file),
            _sought_name(p._sought_name), _tag_name(p._tag_name), _xml(p._xml), _mat(p._mat),
            _nb_dof(p._nb_dof), _dim(p._dim), _nb_nodes(p._nb_nodes), _nb_elements(p._nb_elements),
            _nb_sides(p._nb_sides), _nb_edges(p._nb_edges), _nb_el_nodes(p._nb_el_nodes),
            _nb_sd_nodes(p._nb_sd_nodes), _dof_support(p._dof_support),_theMesh(p._theMesh), _theGrid(nullptr),
            _v(p._v), _parser(p._parser), _ipf(p._ipf)
{
   _scan = p._scan;
   _type = p._type;
   _old_type = p._old_type;
}


XMLParser::~XMLParser() { }


void XMLParser::set(Mesh& ms,
                    int   format)
{
   _theMesh = &ms;
   _format = format;
}


void XMLParser::set(Grid& gr)
{
   _theGrid = &gr;
}


void XMLParser::setMaterialNumber(int m)
{
   _nb_mat = m;
}


void XMLParser::open()
{
   if (_is_opened)
      return;

// read chunks, parse them, and save for one chunk parse
   set_skip_whitespaces(true);
   _is.open(_file.c_str());
   if (_is.fail())
      throw OFELIException("In XMLParser::open(): File " + _file + " cannot be opened.");

// begin parsing
   if (!begin())
      throw OFELIException("In XMLParser::open: Failed to initialize parser for file " + _file);
   string chunk;
   char buff[64];
   while (!_is.eof() && _is.good()) {
      _is.read(buff,64);
      int len = int(_is.gcount());
      _xml.append(buff,len);
      chunk.assign(buff,len);
   }
   if (_is.bad()) {
      _is.close();
      throw OFELIException("In XMLParser::open(): Failed to read file.");
   }
   _is.close();
   if (!end())
      throw OFELIException("In XMLParser::open(): Failed to finalize parsing.");
   _is_opened = true;
}


int XMLParser::scan(size_t k)
{
   if (Verbosity>10 || _scan>1) {
      cout << "Scanning xml file: " << _file << endl;
      cout << "----------------------------------------------------------------------" << endl;
   }
   _scan = k;
   _set_mesh = _set_vector = true;
   if (parse(_xml)) {
      if (Verbosity>10 || _scan>1)
         cout << "Parse done" << endl;
      _scan = 0;
      cout << "----------------------------------------------------------------------" << endl;
      cout << "Scanning complete." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::scan(size_t ): Failed to parse XML file.");
   return -1;
}


int XMLParser::scan(vector<real_t>& t,
                    int             type)
{
   _scan = true;
   _set_mesh = _set_prescription = false;
   _set_vector = true;
   _old_type = _type;
   _type = EType::VECTOR;
   _ft = &t;
   _ft->clear();
   _rtype = type;
   _compact = true;
   if (parse(_xml)) {
      if (Verbosity>10 || _scan>1) {
         cout << "Parse done" << endl;
         cout << "----------------------------------------------------------------------" << endl;
         cout << "Scanning complete." << endl;
      }
      return 0;
   }
   else
      throw OFELIException("In XMLParser::scan(vector<real_t>,int): Failed to parse XML file.");
   return -1;
}


int XMLParser::get(Domain& dm)
{
   _theDomain = &dm;
   _set_mesh = _set_vector = false;
   _set_domain = true;
   _scan = 0;
   _theMesh = nullptr;
   _old_type = _type;
   _type = EType::DDOMAIN;
   _theDomain->_nb_dof = 1;
   _theDomain->_dim = 2;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Domain): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Tabulation& t)
{
   _theTabulation = &t;
   _set_mesh = _set_vector = _set_domain = false;
   _scan = 0;
   _theMesh = nullptr;
   _old_type = _type;
   _type = EType::FUNCTION;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Tabulation): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(DMatrix<real_t>& A)
{
   _old_type = _type;
   _type = EType::MATRIX;
   _theMatrix = &A;
   _set_matrix = true;
   _scan = 0;
   _name = "";
   _theMesh = nullptr;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(DMatrix): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Matrix<real_t>* A)
{
   _old_type = _type;
   _type = EType::MATRIX;
   _theMatrix = A;
   _set_matrix = true;
   _scan = 0;
   _theMesh = nullptr;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Matrix): Failed to parse XML file.");
   return 0;
}


int XMLParser::getArray(Vect<real_t>& A)
{
   _theArray = &A;
   _set_mesh = _set_vector = _set_domain = false;
   _scan = 0;
   _theMesh = nullptr;
   _old_type = _type;
   _type = EType::FUNCTION;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::getArray(Vect): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(IPF& ipf)
{
   _ik1 = _ik2 = _dk1 = _dk2 = _ck = _mk = _pk = _dk = 0;
   _ipf = &ipf;
   _scan = 0;
   _set_mesh = _set_vector = false;
   _theMesh = nullptr;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(IPF): Failed to parse XML file.");
   return 0;
}


int XMLParser::getMaterial()
{
   _set_mesh = _set_vector = false;
   _scan = false;
   _theMesh = nullptr;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::getMaterial(): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(EType                       type,
                   vector<Prescription::PPar>& p)
{
   _prescription_type = type;
   _old_type = _type;
   _type = EType::PRESCRIBE;
   _scan = 0;
   _vp = &p;
   _vp->clear();
   _set_mesh = _set_vector = false;
   _compact = true;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(EType,vector<EType>): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Mesh& ms,
                   int   format)
{
   _set_mesh = true;
   _set_vector = false;
   _scan = 0;
   _theMesh = &ms;
   _nb_nodes = _theMesh->getNbNodes();
   _nb_elements = _theMesh->getNbElements();
   _nb_sides = _theMesh->getNbSides();
   _nb_edges = _theMesh->getNbEdges();
   _format = format;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Mesh,int): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Vect<real_t>& v,
                   const string& name)
{
   _set_mesh = false;
   _set_vector = true;
   _scan = 0;
   _sought_name = name;
   _sought_time = -1.0;
   _nx = _ny = _nz = _nt = 0;
   _compact = true;
   _v = &v;
   _all_steps = 0;
   _name = _sought_name;
   _theMesh = nullptr;
   if (_v->WithMesh())
      _theMesh = &(_v->getMesh());
   _nb_dof = 1;
   _v->setName(_name);
   _old_type = _type;
   _type = EType::VECTOR;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done" << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Vect<real_t>,string): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Vect<real_t>& v,
                   real_t        time,
                   string        name,
                   int           format)
{
   _set_mesh = false;
   _set_vector = true;
   _scan = 0;
   _sought_time = time;
   _sought_name = name;
   _v = &v;
   _theMesh = nullptr;
   if (_v->WithMesh())
      _theMesh = &(_v->getMesh());
   _dof_support = v.getDOFType();
   _format = format;
   _nb_dof = v.getNbDOF();
   _v->setName(_name);
   _compact = true;
   _all_steps = 0;
   _nx = 0, _ny = _nz = 1;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done" << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Vect<real_t>,real_t,string,int): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Mesh&         ms,
                   Vect<real_t>& v,
                   real_t        time,
                   string        name,
                   int           format)
{
   _set_mesh = _set_vector = true;
   _scan = 0;
   _sought_time = time;
   _sought_name = name;
   _v = &v;
   _theMesh = &ms;
   _nb_nodes = _theMesh->getNbNodes();
   _nb_elements = _theMesh->getNbElements();
   _nb_sides = _theMesh->getNbSides();
   _nb_edges = _theMesh->getNbEdges();
   _dof_support = v.getDOFType();
   _nx = 0, _ny = _nz = 1;
   _format = format;
   _compact = true;
   _all_steps = 0;
   _nb_dof = v.getNbDOF();
   _name = v.getName();
   _v->setMesh(*_theMesh,_dof_support,_nb_dof);
   _v->setName(_name);
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done" << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Vect<real_t>,real_t,string,int): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Mesh&                    ms,
                   vector<vector<real_t> >& v,
                   string&                  name)
{
   _set_mesh = _set_vector = true;
   _scan = 0;
   _theMesh = &ms;
   _nb_nodes = _theMesh->getNbNodes();
   _nb_elements = _theMesh->getNbElements();
   _nb_sides = _theMesh->getNbSides();
   _nb_edges = _theMesh->getNbEdges();
   _V = &v;
   _all_steps = 1;
   _compact = true;
   _old_type = _type;
   _type = EType::VECTOR;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done" << endl;
      name = _name;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Mesh,vector<vector<real_t> >,string): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Grid&         gr,
                   Vect<real_t>& v,
                   real_t        time,
                   string        name,
                   int           format)
{
   _set_mesh = _set_grid = _set_vector = true;
   _scan = 0;
   _sought_time = time;
   _sought_name = name;
   _v = &v;
   _theGrid = &gr;
   _nx = _theGrid->getNx(), _ny = _theGrid->getNy(), _nz = _theGrid->getNz();
   _format = format;
   _compact = true;
   _all_steps = 0;
   _nb_dof = v.getNbDOF();
   _name = v.getName();
   _v->setGrid(*_theGrid);
   _v->setName(_name);
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done" << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Grid,Vect<real_t>,real_t,string,int): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Grid&                    gr,
                   vector<vector<real_t> >& v,
                   string&                  name)
{
   _set_grid = _set_vector = true;
   _set_mesh = false;
   _scan = 0;
   _theGrid = &gr;
   _nx = _theGrid->getNx(), _ny = _theGrid->getNy(), _nz = _theGrid->getNz();
   _V = &v;
   _all_steps = 1;
   _compact = true;
   _old_type = _type;
   _type = EType::VECTOR;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done" << endl;
      name = _name;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Grid,vector<vector<real_t> >,string): Failed to parse XML file.");
   return 0;
}


void XMLParser::setFile(string file)
{
   _file = file;
   _set_file = true;
}


bool XMLParser::on_tag_open(string     tag_name,
                            StringMap& attributes)
{
   _nb_el_nodes = 0;
   _tag_name = tag_name;
   switch (_TS[_tag_name]) {

      case Tag::DDOMAIN:
         _set_domain = true;
         _type = EType::DDOMAIN;
         if (_scan>1)
            cout << "Found a domain description." << endl;
         break;

      case Tag::MATERIAL:
         _type = EType::MATERIAL;
         if (_scan>0)
            cout << "Found a material description." << endl;
         break;

      case Tag::PRESCRIPTION:
         _set_prescription = true;
         _old_type = _type;
         _type = EType::PRESCRIBE;
         if (_scan>1)
            cout << "Found a prescription." << endl;
         else {
            _par.code = 0;
            _par.dof = 1;
            _par.bx = _par.by = _par.bz = _par.bt = false;
         }
         break;

      case Tag::SOLUTION:
         if (_scan>2)
            cout << "-> Prescription of exact solution" << endl;
         else
            _par.type = EType::SOLUTION;
         break;

      case Tag::BOUNDARY_CONDITION:
         if (_scan>2)
            cout << "-> Prescription of a boundary condition" << endl;
         else
            _par.type = EType::BOUNDARY_CONDITION;
         break;

      case Tag::BODY_FORCE:
         if (_scan>2)
            cout << "-> Prescription of a body force" << endl;
         else
            _par.type = EType::BODY_FORCE;
         break;

      case Tag::POINT_FORCE:
         if (_scan>2)
            cout << "-> Prescription of a pointwise force" << endl;
         else
            _par.type = EType::POINT_FORCE;
         break;

      case Tag::FLUX:
         if (_scan>2)
            cout << "-> Prescription of a boundary force" << endl;
         else
            _par.type = EType::BOUNDARY_FORCE;
         break;

      case Tag::INITIAL:
         if (_scan>2)
            cout << "-> Prescription of an initial field" << endl;
         else
            _par.type = EType::INITIAL;
         break;

      case Tag::FUNCTION:
         if (_scan>2)
            cout << "-> Definition of a tabulation" << endl;
         break;

      default:
         break;
   }
   for (const auto &a: attributes) {

      if (_ipf && _type==EType::PROJECT)
         read_project(a);

      if (_set_domain && _type==EType::DDOMAIN)
         read_domain(a);

      if (_set_mesh && _type==EType::MESH)
         read_mesh(a);

      if (_set_vector)
         read_vect(a);

      if (_set_matrix)
         read_matrix(a);

      //      if (_material)
      //         read_mat(a);

      if (_set_prescription)
         read_prescription(a);

      if (_type==EType::FUNCTION)
         read_tab(a);
   }
   return true;
}


void XMLParser::read_project(const SString& a)
{
   string mat;
   if (_tag_name=="Project") {
      if (a.first=="name")
         _ipf->_project = a.second;
   }
   if (a.first=="value") {

      switch (_TS[_tag_name]) {
      
         case Tag::VERBOSE:
            if (_scan>1) 
               cout << "Verbosity parameter set." << endl;
            else {
               if (a.first=="value")
                  _ipf->_verbose = stoi(a.second);
            }
            break;

         case Tag::MATERIAL_CODE:
            if (_scan>1) 
               cout << "Material code and name given." << endl;
            else {
               if (a.first=="code")
                  _code = stoi(a.second);
               if (a.first=="material" || a.first=="value")
                  mat = a.second;
               theMaterial.set(_code,mat);
            }
            break;

         case Tag::OUTPUT:
            if (_scan>1) 
               cout << "Output parameter selected." << endl;
            else {
               if (a.first=="value")
                  _ipf->_output = stoi(a.second);
            }
            break;

         case Tag::SAVE:
            if (_scan>1) 
               cout << "Save parameter selected." << endl;
            else {
               if (a.first=="value")
                  _ipf->_save = stoi(a.second);
            }
            break;
   
         case Tag::PLOT:
            if (_scan>1) 
               cout << "Frequency for plotting save parameter selected." << endl;
            else {
               if (a.first=="value")
                  _ipf->_plot = stoi(a.second);
            }
            break;

         case Tag::BC:
            if (_scan>1)
               cout << "Boundary condition toggle selected." << endl;
            else
               _ipf->_bc = stoi(a.second);
            break;

         case Tag::BF:
            if (_scan>1)
               cout << "Body force toggle selected." << endl;
            else
               _ipf->_bf = stoi(a.second);
            break;

         case Tag::SF:
            if (_scan>1)
               cout << "Boundary force toggle selected." << endl;
            else
               _ipf->_sf = stoi(a.second);
            break;

         case Tag::INITIAL:
            if (_scan>1)
               cout << "Initial solution toggle selected." << endl;
            else
               _ipf->_ini = stoi(a.second);
            break;

         case Tag::PRESCRIPTION:
            if (_scan>1)
               cout << "Prescription toggle selected." << endl;
            else
               _ipf->_data = stoi(a.second);
            break;

         case Tag::NB_STEPS:
            if (_scan>1)
               cout << "Number of time steps parameter selected." << endl;
            else
               _ipf->_nb_steps = stoi(a.second);
            break;

         case Tag::NB_ITER:
            if (_scan>1)
               cout << "Number of iterations parameter selected." << endl;
            else
               _ipf->_nb_iter = stoi(a.second);
            break;

         case Tag::TIME_STEP:
            if (_scan>1)
               cout << "Time step parameter selected." << endl;
            else
               _ipf->_time_step = stof(a.second);
            break;

         case Tag::MAX_TIME:
            if (_scan>1)
               cout << "Maximal time parameter selected." << endl;
            else
               _ipf->_max_time = stof(a.second);
            break;

         case Tag::TOLERANCE:
            if (_scan>1)
               cout << "Tolerance parameter selected." << endl;
            else
               _ipf->_tolerance = stof(a.second);
            break;

         case Tag::INTEGER:
            if (_scan>1)
               cout << "Extra integer parameter selected." << endl;
            else
               _ipf->_int_par[_ik1++] = stoi(a.second);
            break;

         case Tag::DOUBLE:
            if (_scan>1)
               cout << "Extra double parameter selected." << endl;
            else
               _ipf->_real_par[_dk1++] = stof(a.second);
            break;

         case Tag::STRING:
            if (_scan>1)
               cout << "Extra string parameter selected." << endl;
            else
               if (a.first=="value")
                  _ipf->_string_par[_dk2++] = a.second;
            break;

         case Tag::COMPLEX:
            if (_scan>1)
               cout << "Extra complex parameter selected." << endl;
            else {
               if (a.first=="real")
                  _ipf->_complex_par[_ck] += complex_t(stof(a.second));
               if (a.first=="imag")
                  _ipf->_complex_par[_ck++] += complex_t(0.,stof(a.second));
            }
            break;

         case Tag::DOMAIN_FILE:
            if (_scan>1)
               cout << "Domain file name given." << endl;
            else
               _ipf->_domain_file = a.second;
            break;

         case Tag::MESH_FILE:
            if (_scan>1)
               cout << "Mesh file name given." << endl;
            else
               _ipf->_mesh_file[_mk++] = a.second;
            break;

         case Tag::INIT_FILE:
            if (_scan>1)
               cout << "Initial field file name given." << endl;
            else
               _ipf->_init_file = a.second;
            break;

         case Tag::RESTART_FILE:
            if (_scan>1) 
               cout << "Restart file name given." << endl;
            else
               _ipf->_restart_file = a.second;
            break;

         case Tag::BC_FILE:
            if (_scan>1) 
               cout << "Boundary condition file name given." << endl;
            else
               _ipf->_bc_file = a.second;
            break;

         case Tag::BF_FILE:
            if (_scan>1)
               cout << "Body force file name given." << endl;
            else
               _ipf->_bf_file = a.second;
            break;

         case Tag::SF_FILE:
            if (_scan>1) 
               cout << "Boundary force file name given." << endl;
            else
               _ipf->_sf_file = a.second;
            break;

         case Tag::SAVE_FILE:
            if (_scan>1) 
               cout << "Save file name given." << endl;
            else
               _ipf->_save_file = a.second;
            break;

         case Tag::PLOT_FILE:
            if (_scan>1) 
               cout << "Plot file name given." << endl;
            else
               _ipf->_plot_file[_pk++] = a.second;
            break;

         case Tag::PRESCRIPTION_FILE:
            if (_scan>1)
               cout << "Prescription file name given." << endl;
            else
               _ipf->_data_file[_dk++] = a.second;
            break;

         case Tag::AUX_FILE:
            if (_scan>1)
               cout << "Auxiliary file name given." << endl;
            else
               _ipf->_aux_file[_ik2++] = a.second;
            break;

         case Tag::DENSITY:
            if (_scan>1)
               cout << "Density value given." << endl;
            else
               _ipf->_mp._density = a.second;
            break;

         case Tag::ELECTRIC_CONDUCTIVITY:
            if (_scan>1) 
               cout << "Electric conductivity value given." << endl;
            else
               _ipf->_mp._electric_cond = a.second;
            break;

         case Tag::ELECTRIC_PERMITTIVITY:
            if (_scan>1)
               cout << "Electric permittivity value given." << endl;
            else
               _ipf->_mp._electric_perm = a.second;
            break;

         case Tag::MAGNETIC_PERMEABILITY:
            if (_scan>1)
               cout << "Magnetic permeability value given." << endl;
            else
               _ipf->_mp._magnetic_perm = a.second;
            break;

         case Tag::POISSON_RATIO:
            if (_scan>1)
               cout << "Poisson ratio value given." << endl;
            else
               _ipf->_mp._poisson = a.second;
            break;

         case Tag::THERMAL_CONDUCTIVITY:
            if (_scan>1)
               cout << "Thermal conductivity value given." << endl;
            else
               _ipf->_mp._thermal_cond = a.second;
            break;

         case Tag::RHO_CP:
            if (_scan>1) 
               cout << "Value of density x specific heat given." << endl;
            else
               _ipf->_mp._rho_cp = a.second;
            break;

         case Tag::VISCOSITY:
            if (_scan>1)
               cout << "Viscosity value given." << endl;
            else
               _ipf->_mp._visc = a.second;
            break;

         case Tag::YOUNG_MODULUS:
            if (_scan>1)
               cout << "Young modulus value given." << endl;
            else
               _ipf->_mp._young = a.second;
            break;

         default:
            break;
      }
   }

   switch (_TS[_tag_name]) {
      
         case Tag::PARAMETER:
            if (_scan>1)
               cout << "keyword parameter selected." << endl;
            else {
               if (a.first=="label")
                  _ipf->_param_label.push_back(a.second);
               if (a.first=="value")
                  _ipf->_param_value.push_back(a.second);
               if (a.first=="ext" || a.first=="extension")
                  _ipf->_param_ext.push_back(a.second);
            }
            break;

         case Tag::ARRAY:
         if (_scan>1)
            cout << "keyword array selected." << endl;
         else {
            if (a.first=="label")
               _ipf->_array_label.push_back(a.second);
            if (a.first=="size")
               _ipf->_array_size.push_back(stoi(a.second));
            if (a.first=="ext" || a.first=="extension")
               _ipf->_array_ext.push_back(a.second);
         }
         break;

         default:
            break;
   }
}


void XMLParser::read_domain(const SString& a)
{
   real_t x=0., y=0., h=0.;
   int code=0;
   switch (_TS[_tag_name]) {

      case Tag::DDOMAIN:
         if (a.first=="dim") {
            _dim = stoi(a.second);
            if (_scan==0)
               _theDomain->_dim = _dim;
         }
         else if (a.first=="nb_dof")
            _nb_dof = stoi(a.second);
         if (_scan==0)
            _theDomain->_nb_dof = _nb_dof;
         break;

      case Tag::VERTEX:
         if (a.first=="x")
            x = stof(a.second);
         if (a.first=="y")
            y = stof(a.second);
         if (a.first=="code")
            code = stoi(a.second);
         if (a.first=="h") {
            h = stof(a.second);
            _theDomain->insertVertex(x,y,h,code);
         }
         break;

      case Tag::LINE:
         break;

      case Tag::CONTOUR:
         break;

      case Tag::HOLE:
         break;

      case Tag::SUBDOMAIN:
         break;

      default:
         break;
   }
}


void XMLParser::read_tab(const SString& a)
{
   switch (_TS[_tag_name]) {

      case Tag::FUNCTION:
         if (_scan)
            cout << "Found a tabulated function:" << endl;
         if (a.first=="name") {
            _theTabulation->setFunction(a.second);
            _nb_funct = _theTabulation->_nb_funct;
            _nb_var = 0;
         }
         break;

      case Tag::VARIABLE:
         if (a.first=="label") {
            _theTabulation->setVariable(a.second);
            _nb_var++; 
         }
         else if (a.first=="nb_points") {
            _tab_size = stoi(a.second);
            if (_scan==0)
               _theTabulation->Funct[_nb_funct-1].Np[_nb_var-1] = _tab_size;
         }
         else if (a.first=="min")
            _theTabulation->Funct[_nb_funct-1].Min[_nb_var-1] = stof(a.second);    
         else if (a.first=="max")
            _theTabulation->Funct[_nb_funct-1].Max[_nb_var-1] = stof(a.second);    
         break;

      case Tag::DATA:
         break;

      default:
         break;
   }
}


void XMLParser::read_matrix(const SString& a)
{
   _type = EType::MATRIX;
   if (_TS[_tag_name]==Tag::MATRIX) {
      if (a.first=="nb_rows")
         _nb_rows = stoi(a.second);
      if (a.first=="nb_cols" || a.first=="nb_columns")
         _nb_cols = stoi(a.second);
      if (a.first=="name") {
         _name = a.second;
         _theMatrix->setName(_name);
      }
   }
}


void XMLParser::read_mat(const SString& a)
{
   if (_scan)
      return;

   switch (_TS[_tag_name]) {

      case Tag::MATERIAL:
         if (a.first=="name") {
         theMaterial._name = a.second;
         theMaterial._density[_nb_mat].exist = false;
         theMaterial._specific_heat[_nb_mat].exist = false;
         theMaterial._thermal_conductivity[_nb_mat].exist = false;
         theMaterial._melting_temperature[_nb_mat].exist = false;
         theMaterial._evaporation_temperature[_nb_mat].exist = false;
         theMaterial._thermal_expansion[_nb_mat].exist = false;
         theMaterial._latent_heat_melting[_nb_mat].exist = false;
         theMaterial._latent_heat_evaporation[_nb_mat].exist = false;
         theMaterial._dielectric_constant[_nb_mat].exist = false;
         theMaterial._electric_conductivity[_nb_mat].exist = false;
         theMaterial._electric_resistivity[_nb_mat].exist = false;
         theMaterial._magnetic_permeability[_nb_mat].exist = false;
         theMaterial._viscosity[_nb_mat].exist = false;
         theMaterial._young_modulus[_nb_mat].exist = false;
         theMaterial._poisson_ratio[_nb_mat].exist = false;
         }
         break;

      case Tag::DENSITY:
         if (_scan>1)
            cout << "   Found density value." << endl;
         if (_scan==0) {
            theMaterial._density[_nb_mat].exist = true;
            theMaterial._density[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._density[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::SPECIFIC_HEAT:
         if (_scan>1)
            cout << "   Found specific heat value." << endl;
         if (_scan==0) {
            theMaterial._specific_heat[_nb_mat].exist = true;
            theMaterial._specific_heat[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._specific_heat[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::THERMAL_CONDUCTIVITY:
         if (_scan>1)
            cout << "   Found thermal conductivity value." << endl;
         if (_scan==0) {
            theMaterial._thermal_conductivity[_nb_mat].exist = true;
            theMaterial._thermal_conductivity[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._thermal_conductivity[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::MELTING_TEMPERATURE:
         if (_scan>1)
            cout << "   Found melting temperature value." << endl;
         if (_scan==0) {
            theMaterial._melting_temperature[_nb_mat].exist = true;
            theMaterial._melting_temperature[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._melting_temperature[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::EVAPORATION_TEMPERATURE:
         if (_scan>1)
            cout << "   Found evaporation temperature value." << endl;
         if (_scan==0) {
            theMaterial._evaporation_temperature[_nb_mat].exist = true;
            theMaterial._evaporation_temperature[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._evaporation_temperature[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::THERMAL_EXPANSION:
         if (_scan>1)
            cout << "   Found thermal expansion value." << endl;
         if (_scan==0) {
            theMaterial._thermal_expansion[_nb_mat].exist = true;
            theMaterial._thermal_expansion[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._thermal_expansion[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::LATENT_HEAT_MELTING:
         if (_scan>1)
            cout << "   Found latent heat for melting value." << endl;
         if (_scan==0) {
            theMaterial._latent_heat_melting[_nb_mat].exist = true;
            theMaterial._latent_heat_melting[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._latent_heat_melting[_nb_mat-1].value = stof(a.second);
         }
         break;

      case Tag::LATENT_HEAT_EVAPORATION:
         if (_scan>1)
            cout << "   Found latent heat for evaporation value." << endl;
         if (_scan==0) {
            theMaterial._latent_heat_evaporation[_nb_mat].exist = true;
            theMaterial._latent_heat_evaporation[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._latent_heat_evaporation[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::DIELECTRIC_CONSTANT:
         if (_scan>1)
            cout << "   Found dielectric constant value." << endl;
         if (_scan==0) {
            theMaterial._dielectric_constant[_nb_mat].exist = true;
            theMaterial._dielectric_constant[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._dielectric_constant[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::ELECTRIC_CONDUCTIVITY:
         if (_scan>1)
            cout << "   Found electric conductivity value." << endl;
         if (_scan==0) {
            theMaterial._electric_conductivity[_nb_mat].exist = true;
            theMaterial._electric_conductivity[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._electric_conductivity[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::ELECTRIC_RESISTIVITY:
         if (_scan>1)
            cout << "   Found electric resistivity value." << endl;
         if (_scan==0) {
            theMaterial._electric_resistivity[_nb_mat].exist = true;
            theMaterial._electric_resistivity[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._electric_resistivity[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::MAGNETIC_PERMEABILITY:
         if (_scan>1)
            cout << "   Found magnetic permeability value." << endl;
         if (_scan==0) {
            theMaterial._magnetic_permeability[_nb_mat].exist = true;
            theMaterial._magnetic_permeability[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._magnetic_permeability[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::VISCOSITY:
         if (_scan>1)
            cout << "   Found viscosity value." << endl;
         if (_scan==0) {
            theMaterial._viscosity[_nb_mat].exist = true;
            theMaterial._viscosity[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._viscosity[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::YOUNG_MODULUS:
         if (_scan>1)
            cout << "   Found Young modulus value." << endl;
         if (_scan==0) {
            theMaterial._young_modulus[_nb_mat].exist = true;
            theMaterial._young_modulus[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._young_modulus[_nb_mat].value = stof(a.second);
         }
         break;

      case Tag::POISSON_RATIO:
         if (_scan>1)
            cout << "   Found Poisson ratio value." << endl;
         if (_scan==0) {
            theMaterial._poisson_ratio[_nb_mat].exist = true;
            theMaterial._poisson_ratio[_nb_mat].type = BY_VALUE;
            if (a.first=="value")
               theMaterial._poisson_ratio[_nb_mat].value = stof(a.second);
         }
         break;
   
      default:
         break;
   }
}


void XMLParser::read_mesh(const SString& a)
{
   _type = EType::MESH;
   switch (_TS[_tag_name]) {

      case Tag::MESH:
         if (_dim==1)
            _el_shape = "line", _nb_el_nodes = 2;
         else if (_dim==2)
            _el_shape = "triangle", _nb_el_nodes = 3, _sd_shape = "line", _nb_sd_nodes = 2;
         else if (_dim==3)
            _el_shape = "tetrahedron", _nb_el_nodes = 4, _sd_shape = "triangle", _nb_sd_nodes = 3;
         if (a.first=="file" && _scan==0)
            new XMLParser(a.second,*_theMesh);
         else if (a.first=="dim") {
            _dim = stoi(a.second);
            if (_scan==0)
               _theMesh->setDim(_dim);
         }
         else if (a.first=="nb_dof")
            _nb_dof = stoi(a.second);
         if (_scan)
            cout << "Found a finite element mesh:" << endl;
         break;

      case Tag::ELEMENTS:
         if (a.first=="shape")
            _el_shape = a.second;
         if (a.first=="nodes")
            _nb_el_nodes = stoi(a.second);
         break;

      case Tag::SIDES:
         if (a.first=="shape")
            _sd_shape = a.second;
         if (a.first=="nodes")
            _nb_sd_nodes = stoi(a.second);
         break;

      case Tag::MATERIAL:
         break;

      default:
         break;
   }
}


void XMLParser::read_prescription(const pair<string,string>& a)
{
   static vector<string> bccs {"PERIODIC_A","PERIODIC_B","CONTACT","CONTACT_M","CONTACT_S","SLIP"};
   int bcc[6] = {9999,-9999,9998,9997,-9997,9996};
 
   if (a.first=="file" && _scan==0)
      new XMLParser(a.second,*_theMesh);
   if (a.first=="code")
      _par.code = BoundaryConditionCode(bccs,bcc,a.second);
   else if (a.first=="dof")
      _par.dof = stoi(a.second);
   else if (a.first=="x")
      _par.x = stof(a.second), _par.bx = true;
   else if (a.first=="y")
      _par.y = stof(a.second), _par.by = true;
   else if (a.first=="z")
      _par.z = stof(a.second), _par.bz = true;
   else if (a.first=="time")
      _par.t = stof(a.second), _par.bt = true;
   else if (a.first=="value") {
   }
   else
      ;
}


void XMLParser::read_vect(const SString& a)
{
   _vect_size = 0;
   if ((_tag_name=="Vector"||_tag_name=="Field") && 
       (_type==EType::MESH||_type==EType::GRID||_type==EType::VECTOR)) {
      if (a.first=="file" && _scan==0) {
         _parser = new XMLParser(a.second);
         _parser->set(*_theMesh);
         _parser->_v = _v;
      }
      else if (a.first=="name")
         _name = a.second;
      else if (a.first=="type") {
         if (a.second=="Node")
            _dof_support = NODE_DOF;
         else if (a.second=="Element")
            _dof_support = ELEMENT_DOF;
         else if (a.second=="Side")
            _dof_support = SIDE_DOF;
         else if (a.second=="Edge")
            _dof_support = EDGE_DOF;
         else
            ;
      }
      else if (a.first=="nb_dof")
         _nb_dof = stoi(a.second);
      else if (a.first=="nx")
         _nx = stoi(a.second);
      else if (a.first=="ny")
         _ny = stoi(a.second);
      else if (a.first=="nz")
         _nz = stoi(a.second);
      else if (a.first=="nt")
         _nt = stoi(a.second);
      else if (a.first=="size")
         _vect_size = stoi(a.second);
   }
   if (_tag_name=="Step") {
      _compact = false;
      if (a.first=="time")
         _time = stof(a.second);
   }
   else if (_tag_name=="constant" || _tag_name=="const") {
      _dof = 0;
      if (a.first=="value")
         _val = stof(a.second);
      else if (a.first=="dof")
         _dof = stoi(a.second);
   }
   else if (_tag_name=="expression" || _tag_name=="expr") {
      _dof = 0;
      if (a.first=="dof")
         _dof = stoi(a.second);
   }
}


bool XMLParser::on_cdata(string cdata)
{
   _name = "";
   vector<string> tokens;
   {
      string buf;
      std::stringstream ss(cdata);
      while (ss >> buf)
         tokens.push_back(buf);
   }
   vector<string>::iterator it=tokens.begin();
   if (_type==EType::PROJECT)
      read_project_data(tokens,it);

   else if (_type==EType::MESH)
      read_mesh_data(tokens,it);

   else if (_type==EType::DDOMAIN)
      read_domain_data(tokens,it);

// Material in mesh data
   else if (_tag_name=="Material" && _set_mesh) {
      while (it!=tokens.end()) {
         _cm = stoi(*it++);
         _mat = *it++;
         _nb_mat = theMaterial.set(_cm,_mat);
      }
   }

// Material
   else if (_type==EType::MATERIAL)
      read_mat_data(tokens,it);

// Prescribe
   else if (_set_prescription && _type==EType::PRESCRIBE && _prescription_type==_par.type)
//   else if (_set_prescription && _prescription_type==_par.type)
      read_prescribe_data(tokens,it);

// Vector
   else if ((_tag_name=="Step"||(_tag_name=="Field"&&_compact)||(_tag_name=="Vector"&&_compact)) && _set_vector)
      read_vect_data(tokens,it);
   else if ((_tag_name=="constant"||_tag_name=="const") && _set_vector)
      read_const_vect_data(tokens,it);
   else if ((_tag_name=="expression"||_tag_name=="expr") && _set_vector)
      read_exp_vect_data(tokens,it);

// Tabulated function
   else if (_type==EType::FUNCTION) {
      if (_tag_name=="Data")
         read_tab_data(tokens,it);
   }

// Matrix
   else if (_type==EType::MATRIX)
      read_matrix_data(tokens,it);

   else
      ;
   return true;
}


void XMLParser::read_prescribe_data(const vector<string>&     tokens,
                                    vector<string>::iterator& it)
{
   _par.fct = tokens[0];
   for (size_t i=1; i<tokens.size(); ++i)
      _par.fct += tokens[i];
   _vp->push_back(_par);
}


void XMLParser::read_vect_data(const vector<string>&     tokens,
                               vector<string>::iterator& it)
{
   size_t nb = 0;
   if (_theMesh!=nullptr) {
      if (_dof_support==NODE_DOF)
         nb = _nb_nodes*_nb_dof;
      else if (_dof_support==ELEMENT_DOF)
         nb = _nb_elements*_nb_dof;
      else if (_dof_support==SIDE_DOF)
         nb = _nb_sides*_nb_dof;
      else if (_dof_support==EDGE_DOF)
         nb = _nb_edges*_nb_dof;
      else
         nb = _nx*_ny*_nz;
   }

   if (_scan) {
      if (_compact) {
         while (it!=tokens.end()) {
            _time = stof(*it++);
            _ft->push_back(_time);
            for (size_t n=0; n<nb; n++)
               *it++;
         }
      }
      else
         _ft->push_back(_time);
      return;
   }
   if (_all_steps>0 && _compact) {
      size_t i = 0;
      while (it!=tokens.end()) {
         _time = stof(*it++);
         for (size_t j=0; j<nb; j++)
            (*_V)[i].push_back(stof(*it++));
         i++;
      }
   }
   else {
      if (_compact) {
         if ((_name==_sought_name || _sought_name=="ANYTHING") && _set_vector) {
            if (_nx>0) {
               if (_ny==0) {
                  _vect_size = _nx;
                  _v->setSize(_nx);
                  for (size_t i=0; i<_vect_size; ++i)
                     (*_v)[i] = stof(*it++);
               }
               else if (_nz==0) {
                  _vect_size = _nx*_ny;
                  _v->setSize(_nx,_ny);
                  for (size_t i=0; i<_vect_size; ++i)
                     (*_v)[i] = stof(*it++);
               }
               else if (_nt==0) {
                  _vect_size = _nx*_ny*_nz;
                  _v->setSize(_nx,_ny,_nz);
                  for (size_t i=0; i<_vect_size; ++i)
                     (*_v)[i] = stof(*it++);
               }
               else {
                  _vect_size = _nx*_ny*_nz*_nt;
                  _v->setSize(_nx,_ny,_nz,_nt);
                  for (size_t i=0; i<_vect_size; ++i)
                     (*_v)[i] = stof(*it++);
               }
            }
            else {
               if (nb > 0) {
                  while (it!=tokens.end()) {
                     _time = stof(*it++);
                     if (_time==_sought_time || _sought_time==-1.0) {
                        _v->setMesh(*_theMesh,_dof_support,_nb_dof);
                        _v->setName(_name);
                        _v->setTime(_time);
                        for (size_t n=1; n<=_v->getNb(); ++n)
                           for (size_t k=1; k<=_nb_dof; ++k)
                              (*_v)(n,k) = stof(*it++);
                     }
                     else {
                        for (size_t n=1; n<=_v->getNb()*_nb_dof; ++n)
                           *it++;
                     }
                  }
               }
               else {
                  vector<real_t> V;
                  _time = stof(*it++);
                  while (it!=tokens.end())
                     V.push_back(stof(*it++));
                  _v->setSize(V.size());
                  for (size_t i=0; i<V.size(); ++i)
                     (*_v)[i] = V[i];
               }
            }
         }
      }
      else {
         if ((_name==_sought_name || _sought_name=="ANYTHING") &&
             (_time==_sought_time || _sought_time==-1.0) && _set_vector) {
            _v->setMesh(*_theMesh,_dof_support,_nb_dof);
            _v->setName(_name);
            _v->setTime(_time);
            for (size_t n=1; n<=_v->getNb(); ++n)
               for (size_t k=1; k<=_nb_dof; ++k)
                  (*_v)(n,k) = stof(*it++);
         }
      }
   }
}


void XMLParser::read_const_vect_data(const vector<string>&     tokens,
                                     vector<string>::iterator& it)
{
   if (_scan) {
      _ft->push_back(_time);
      if (_scan>1) {
         if (_dof_support==NODE_DOF)
            cout << "Found a nodewise field vector, Name: " << _name << ", Time = " 
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==ELEMENT_DOF)
            cout << "Found an elementwise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==SIDE_DOF)
            cout << "Found a sidewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==EDGE_DOF)
            cout << "Found an edgewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else
            ;
      }
   }
   else {
      if ((_name==_sought_name || _sought_name=="ANYTHING") &&
          (_time==_sought_time || _sought_time==-1.0) && _set_vector) {
         _v->setTime(_time);
         _val = stof(*it++);
         for (size_t k=1; k<=_v->getNb(); k++) {
            if (_dof==0)
               for (size_t l=1; l<=_v->getNbDOF(); l++)
                  (*_v)(k,l) = _val;
            else
               (*_v)(k,_dof) = _val;
         }
      }
   }
}


void XMLParser::read_const_vect_data()
{
   if (_scan) {
      _ft->push_back(_time);
      if (_scan>1) {
         if (_dof_support==NODE_DOF)
            cout << "Found a nodewise field vector, Name: " << _name << ", Time = " 
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==ELEMENT_DOF)
            cout << "Found an elementwise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==SIDE_DOF)
            cout << "Found a sidewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==EDGE_DOF)
            cout << "Found an edgewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else
            ;
      }
      return;
   }

   vector<string>::iterator it;
   if ((_name==_sought_name || _sought_name=="ANYTHING") &&
       (_time==_sought_time || _sought_time==-1.0) && _set_vector) {
      _v->setTime(_time);
      _val = stof(*it++);
      for (size_t k=1; k<=_v->getNb(); k++) {
         if (_dof==0) {
            for (size_t l=1; l<=_v->getNbDOF(); l++)
               (*_v)(k,l) = _val;
         }
         else
            (*_v)(k,_dof) = _val;
      }
   }
}


void XMLParser::read_exp_vect_data(const vector<string>&     tokens,
                                   vector<string>::iterator& it)
{
   if (_scan) {
      _ft->push_back(_time);
      if (_scan>1) {
         if (_dof_support==NODE_DOF)
            cout << "Found a nodewise field vector, Name: " << _name << ", Time = " 
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==ELEMENT_DOF)
            cout << "Found an elementwise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==SIDE_DOF)
            cout << "Found a sidewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==EDGE_DOF)
            cout << "Found an edgewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else
            ;
      }
      return;
   }

   if ((_name==_sought_name || _sought_name=="ANYTHING") &&
       (_time==_sought_time || _sought_time==-1.0) && _set_vector) {
      _v->setTime(_time);
      _theFct.set(*it++);
      for (size_t n=1; n<=_v->getNb(); n++) {
         if (_dof==0)
            for (size_t k=1; k<=_v->getNbDOF(); k++)
               parse_exp(n,k);
         else
            parse_exp(n,_dof);
      }
   }
}


void XMLParser::read_tab_data(const vector<string>&     tokens,
                              vector<string>::iterator& it)
{
   _theTabulation->setSizes();
   size_t i=0;
   while (it!=tokens.end())
      _theTabulation->Funct[_nb_funct-1].Val(++i) = stof(*it++);
}


void XMLParser::read_matrix_data(const vector<string>&     tokens,
                                 vector<string>::iterator& it)
{
   _theMatrix->setSize(_nb_rows,_nb_cols);
   size_t i=1, j=0;
   while (it!=tokens.end()) {
      if (++j>_nb_cols && i<_nb_rows)
         i++, j = 1;
      if (i>_nb_rows)
         throw OFELIException("In XMLParser::read_matrix_data(..): Insufficient number of matrix entries.");
      (*_theMatrix)(i,j) = stof(*it++);
   }
}


void XMLParser::parse_exp(size_t n,
                          size_t k)
{
   Point<real_t> c;
   if (_dof_support==NODE_DOF)
      c = (*_theMesh)[n]->getCoord();
   else if (_dof_support==SIDE_DOF) {
      theSide = _theMesh->getPtrSide(n);
      if (theSide->getShape()==LINE)
	 c = Line2(theSide).getCenter();
      else if (theSide->getShape()==TRIANGLE)
         c = Triang3(theSide).getCenter();
      else if (theSide->getShape()==QUADRILATERAL)
         c = Quad4(theSide).getCenter();
   }
   else if (_dof_support==ELEMENT_DOF) {
      theElement = (*_theMesh)(n);
      if (theElement->getShape()==LINE)
         c = Line2(theElement).getCenter();
      else if (theElement->getShape()==TRIANGLE)
         c = Triang3(theElement).getCenter();
      else if (theElement->getShape()==QUADRILATERAL)
         c = Quad4(theElement).getCenter();
      else if (theElement->getShape()==TETRAHEDRON)
         c = Tetra4(theElement).getCenter();
      else if (theElement->getShape()==HEXAHEDRON)
         c = Hexa8(theElement).getCenter();
      else if (theElement->getShape()==PENTAHEDRON)
         c = Penta6(theElement).getCenter();
   }
   (*_v)(n,k) = _theFct(c);
}


void XMLParser::read_mat_data(const vector<string>&     tokens,
                              vector<string>::iterator& it)
{
   switch (_TS[_tag_name]) {

      case Tag::DENSITY:
         theMaterial._density[_nb_mat].exist = true;
         theMaterial._density[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._density[_nb_mat].value = stoi(*it);
         break;

      case Tag::SPECIFIC_HEAT:
         theMaterial._specific_heat[_nb_mat].exist = true;
         theMaterial._specific_heat[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._specific_heat[_nb_mat].value = stof(*it);
         break;

      case Tag::THERMAL_CONDUCTIVITY:
         theMaterial._thermal_conductivity[_nb_mat].exist = true;
         theMaterial._thermal_conductivity[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._thermal_conductivity[_nb_mat].value = stof(*it);
         break;

      case Tag::MELTING_TEMPERATURE:
         theMaterial._melting_temperature[_nb_mat].exist = true;
         theMaterial._melting_temperature[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._melting_temperature[_nb_mat].value = stof(*it);
         break;

      case Tag::EVAPORATION_TEMPERATURE:
         theMaterial._evaporation_temperature[_nb_mat].exist = true;
         theMaterial._evaporation_temperature[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._evaporation_temperature[_nb_mat].value = stof(*it);
         break;

      case Tag::THERMAL_EXPANSION:
         theMaterial._thermal_expansion[_nb_mat].exist = true;
         theMaterial._thermal_expansion[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._thermal_expansion[_nb_mat].value = stof(*it);
         break;

      case Tag::LATENT_HEAT_MELTING:
         theMaterial._latent_heat_melting[_nb_mat].exist = true;
         theMaterial._latent_heat_melting[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._latent_heat_melting[_nb_mat].value = stof(*it);
         break;

      case Tag::LATENT_HEAT_EVAPORATION:
         theMaterial._latent_heat_evaporation[_nb_mat].exist = true;
         theMaterial._latent_heat_evaporation[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._latent_heat_evaporation[_nb_mat].value = stof(*it);
         break;

      case Tag::DIELECTRIC_CONSTANT:
         theMaterial._dielectric_constant[_nb_mat].exist = true;
         theMaterial._dielectric_constant[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._dielectric_constant[_nb_mat].value = stof(*it);
         break;

      case Tag::ELECTRIC_CONDUCTIVITY:
         theMaterial._electric_conductivity[_nb_mat].exist = true;
         theMaterial._electric_conductivity[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._electric_conductivity[_nb_mat].value = stof(*it);
         break;

      case Tag::ELECTRIC_RESISTIVITY:
         theMaterial._electric_resistivity[_nb_mat].exist = true;
         theMaterial._electric_resistivity[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._electric_resistivity[_nb_mat].value = stof(*it);
         break;

      case Tag::MAGNETIC_PERMEABILITY:
         theMaterial._magnetic_permeability[_nb_mat].exist = true;
         theMaterial._magnetic_permeability[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._magnetic_permeability[_nb_mat].value = stof(*it);
         break;

      case Tag::VISCOSITY:
         theMaterial._viscosity[_nb_mat].exist = true;
         theMaterial._viscosity[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._viscosity[_nb_mat].value = stof(*it);
         break;

      case Tag::YOUNG_MODULUS:
         theMaterial._young_modulus[_nb_mat].exist = true;
         theMaterial._young_modulus[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._young_modulus[_nb_mat].value = stof(*it);
         break;

      case Tag::POISSON_RATIO:
         theMaterial._poisson_ratio[_nb_mat].exist = true;
         theMaterial._poisson_ratio[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._poisson_ratio[_nb_mat].value = stof(*it);
         break;

      default:
         break;
   }
}


void XMLParser::read_project_data(const vector<string>&     tokens,
                                  vector<string>::iterator& it)
{
   switch (_TS[_tag_name]) {

      case Tag::VERBOSE:
         IPF_TOKENI(_ipf->_verbose)

      case Tag::MATERIAL_CODE:
         if (it!=tokens.end())
            theMaterial.set(_code,it->c_str());
         break;

      case Tag::OUTPUT:
         IPF_TOKENI(_ipf->_output)

      case Tag::SAVE:
         IPF_TOKENI(_ipf->_save)

      case Tag::PLOT:
         IPF_TOKENI(_ipf->_plot)

      case Tag::BC:
         IPF_TOKENI(_ipf->_bc)

      case Tag::BODY_FORCE:
         IPF_TOKENI(_ipf->_bf)

      case Tag::FLUX:
         IPF_TOKENI(_ipf->_sf)

      case Tag::INITIAL:
         IPF_TOKENI(_ipf->_ini)

      case Tag::DATA:
         IPF_TOKENI(_ipf->_data)
      
      case Tag::NB_STEPS:
         IPF_TOKENI(_ipf->_nb_steps)

      case Tag::TIME_STEP:
         IPF_TOKENF(_ipf->_time_step)

      case Tag::MAX_TIME:
         IPF_TOKENF(_ipf->_max_time)

      case Tag::NB_ITER:
         IPF_TOKENI(_ipf->_data)

      case Tag::TOLERANCE:
         IPF_TOKENF(_ipf->_tolerance)

      case Tag::INTEGER:
         IPF_TOKENI(_ipf->_int_par[_ik1++])

      case Tag::DOUBLE:
         IPF_TOKENF(_ipf->_real_par[_dk1++])

      case Tag::STRING:
         IPF_TOKEN(_ipf->_string_par[_dk2++])

      case Tag::COMPLEX:
         while (it!=tokens.end()) {
            _ipf->_complex_par[_ck]  = complex_t(stof(*it++));
            _ipf->_complex_par[_ck] += complex_t(0.,stof(*it++));
         }
         _ck++;
         break;

      case Tag::PARAMETER:
         if (it!=tokens.end())
            _ipf->_param_value.push_back(*it);
         break;

      case Tag::DOMAIN_FILE:
         IPF_TOKEN(_ipf->_domain_file)

      case Tag::MESH_FILE:
         IPF_TOKEN(_ipf->_mesh_file[_mk++])

      case Tag::INIT_FILE:
         IPF_TOKEN(_ipf->_init_file)

      case Tag::RESTART_FILE:
         IPF_TOKEN(_ipf->_restart_file)

      case Tag::BC_FILE:
         IPF_TOKEN(_ipf->_bc_file)

      case Tag::BF_FILE:
         IPF_TOKEN(_ipf->_bf_file)

      case Tag::SF_FILE:
         IPF_TOKEN(_ipf->_sf_file)

      case Tag::SAVE_FILE:
         IPF_TOKEN(_ipf->_save_file)

      case Tag::PLOT_FILE:
         IPF_TOKEN(_ipf->_plot_file[_pk++])

      case Tag::PRESCRIPTION_FILE:
         IPF_TOKEN(_ipf->_data_file[_dk++])

      case Tag::AUX_FILE:
         IPF_TOKEN(_ipf->_aux_file[_ik2++])

      case Tag::DENSITY:
         IPF_TOKEN(_ipf->_mp._density)

      case Tag::ELECTRIC_CONDUCTIVITY:
         IPF_TOKEN(_ipf->_mp._electric_cond)

      case Tag::ELECTRIC_PERMITTIVITY:
         IPF_TOKEN(_ipf->_mp._electric_perm)

      case Tag::MAGNETIC_PERMEABILITY:
         IPF_TOKEN(_ipf->_mp._magnetic_perm)

      case Tag::POISSON_RATIO:
         IPF_TOKEN(_ipf->_mp._poisson)

      case Tag::THERMAL_CONDUCTIVITY:
         IPF_TOKEN(_ipf->_mp._thermal_cond)

      case Tag::RHO_CP:
         IPF_TOKEN(_ipf->_mp._rho_cp)

      case Tag::VISCOSITY:
         IPF_TOKEN(_ipf->_mp._visc)

      case Tag::YOUNG_MODULUS:
         IPF_TOKEN(_ipf->_mp._young)

      case Tag::ARRAY:
         {
            size_t k = 0;
            while (it!=tokens.end() && _scan==false)
               (_ipf->_array_value[k++]).push_back(stof(*it++));
            _ipf->_array_size.push_back(k);
            break;
         }

      default:
         break;
   }
}


void XMLParser::read_domain_data(const vector<string>&     tokens,
                                 vector<string>::iterator& it)
{
   if (!_set_domain)
      return;
   vector<size_t> c;
   switch (_TS[_tag_name]) {

      case Tag::VERTEX:
         while (it!=tokens.end()) {
               real_t x[3];
            for (size_t j=0; j<_dim; j++)
               x[j] = stof(*it++);
            int c1 = stoi(*it++);
            real_t h = stof(*it++);
            _theDomain->insertVertex(x[0],x[1],h,c1);
         }
         break;

      case Tag::LINE:
         while (it!=tokens.end()) {
            size_t v1 = stoi(*it++);
            size_t v2 = stoi(*it++);
            int dc = stoi(*it++), nc = 0;
            if (dc<0)
               nc = -dc, dc = 0;
            _theDomain->insertLine(v1,v2,dc,nc);
         }
         break;

      case Tag::CIRCLE:
         while (it!=tokens.end()) {
            size_t n1 = stoi(*it++);
            size_t n2 = stoi(*it++);
            size_t n3 = stoi(*it++);
            int dc = stoi(*it++), nc = dc;
            _theDomain->insertCircle(n1,n2,n3,dc,nc);
         }
         break;

      case Tag::CONTOUR:
         c.clear();
         while (it!=tokens.end())
            c.push_back(stoi(*it++));
         _theDomain->insertContour(c);
         break;

      case Tag::HOLE:
         c.clear();
         while (it!=tokens.end())
            c.push_back(stoi(*it++));
         _theDomain->insertHole(c);
         break;

      case Tag::REQUIRED_VERTEX:
         while (it!=tokens.end())
            _theDomain->insertRequiredVertex(stoi(*it++));
         break;

      case Tag::REQUIRED_EDGE:
         while (it!=tokens.end())
            _theDomain->insertRequiredEdge(stoi(*it++));
         break;

      case Tag::SUBDOMAIN:
         while (it!=tokens.end()) {
            size_t ln = stoi(*it++);
            int orient = stoi(*it++);
            int ref = stoi(*it++);
            _theDomain->insertSubDomain(ln,orient,ref);
         }
         break;

      case Tag::RECTANGLE:
         real_t x[4];
         while (it!=tokens.end()) {
            x[0] = stof(*it++);
            x[1] = stof(*it++);
            x[2] = stof(*it++);
            x[3] = stof(*it++);
            size_t n1 = stoi(*it++);
            size_t n2 = stoi(*it++);
            real_t r = stof(*it++);
            int c1 = stoi(*it++);
            int c2 = stoi(*it++);
            int c3 = stoi(*it++);
            int c4 = stoi(*it++);
            string file = *it++;
            _theDomain->Rectangle(x,n1,n2,r,c1,c2,c3,c4);
            _theMesh->put(file);
         }
         break;

      default:
         break;
   }
}


void XMLParser::read_mesh_data(const vector<string>&     tokens,
                               vector<string>::iterator& it)
{
   int code[MAX_NBDOF_NODE];
   size_t first_dof = 1;
   if (!_set_mesh)
      return;

   switch (_TS[_tag_name]) {

      case Tag::NODES:
         first_dof = 1;
         _nb_nodes = 0;
         while (it!=tokens.end()) {
            Point<real_t> a;
            a.x = stof(*it++);
            if (_dim > 1)
               a.y = stof(*it++);
            if (_dim > 2)
               a.z = stof(*it++);
            int mark = stoi(*it++);
            _nb_nodes++;
            if (_scan==0) {
               Node *nd = new Node(_nb_nodes,a);
               nd->setNbDOF(_nb_dof);
               nd->setDOF(first_dof,_nb_dof);
               DOFCode(mark,_nb_dof,code);
               nd->setCode(code);
               _theMesh->Add(nd);
            }
         }
         break;

      case Tag::ELEMENTS:
         _nb_elements = 0;
         while (it!=tokens.end()) {
            int nnd[20];
            for (size_t j=0; j<_nb_el_nodes; j++)
               nnd[j] = stoi(*it++);
            int code = stoi(*it++);
            _nb_elements++;
            if (_scan==0) {
               Element *el = new Element(_nb_elements,_el_shape,code);
               for (size_t k=0; k<_nb_el_nodes; k++)
                  el->Add((*_theMesh)[nnd[k]]);
               _theMesh->Add(el);
            }
         }
         break;

      case Tag::SIDES:
         first_dof = 1;
         _nb_sides = 0;
         while (it!=tokens.end()) {
            int nnd[8];
            for (size_t j=0; j<_nb_sd_nodes; j++)
               nnd[j] = stoi(*it++);
            int mark = stoi(*it++);
            _nb_sides++;
            if (_scan==0) {
               Side *sd = new Side(_nb_sides,_sd_shape);
               for (size_t k=0; k<_nb_sd_nodes; k++)
                  sd->Add((*_theMesh)[nnd[k]]);
               sd->setNbDOF(_nb_dof);
               sd->setDOF(first_dof,_nb_dof);
               DOFCode(mark,_nb_dof,code);
               for (size_t l=0; l<_nb_dof; l++)
                  sd->setCode(l+1,code[l]);
               _theMesh->Add(sd);
            }
         }
         break;

      case Tag::MATERIAL:
         while (it!=tokens.end()) {
            int n = stoi(*it++);
            string name = *it++;
            _nb_mat = theMaterial.set(n,name);
         }
         break;

      default:
         break;
   }
}


bool XMLParser::on_tag_close(string tag_name)
{
   if (Verbosity>10)
      cout << "CLOSING TAG: " << tag_name << endl;
   if (tag_name=="Constant" || tag_name=="Const")
      read_const_vect_data();
   switch (_TS[tag_name]) {

      case Tag::DDOMAIN:
         _set_domain = false;
         _type = _old_type;
         break;

      case Tag::MESH:
         _set_mesh = false;
         _type = _old_type;
         break;

      case Tag::PRESCRIPTION:
         _set_prescription = false;
         _type = _old_type;
         break;

      case Tag::FUNCTION:
         _type = _old_type;
         break;

      default:
         break;
   }

   return true;
}


bool XMLParser::on_comment(string comment)
{
   if (Verbosity>10)
      cout << "Parsing document: COMMENT: " << comment << endl;
   return true;
}


bool XMLParser::on_processing(string value)
{
   return true;
}


bool XMLParser::on_doctype(string value)
{
   if (Verbosity > 10)
      cout << "Parsing document: Doctype: " << value << endl;
   return true;
}


bool XMLParser::on_document_begin()
{
   if (Verbosity > 10)
      cout << "Parsing document: Begin" << endl;
   return true;
}


bool XMLParser::on_document_end()
{
   if (Verbosity > 10)
      cout << "Parsing document: End" << endl;
   return true;
}


void XMLParser::on_error(int    e,
                         int    line,
                         int    col,
                         string message)
{
   cout << "ERROR(" << e << "): " << message << ", at " << line << ": " << col << endl;
}

} /* namespace OFELI */
