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

                      Implementation of class 'XMLParser'

  ==============================================================================*/


#include "io/XMLParser.h"
#include "util/macros.h"
#include "mesh/Domain.h"
#include "mesh/Mesh.h"
#include "mesh/MeshUtil.h"
#include "io/IPF.h"
#include "mesh/Material.h"
#include "shape_functions/Line2.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Quad4.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Hexa8.h"
#include "shape_functions/Penta6.h"
#include "io/Tabulation.h"
#include "equations/AbsEqua.h"
#include "linear_algebra/Vect_impl.h"
#include "OFELIException.h"

namespace OFELI {

extern Material theMaterial;


XMLParser::XMLParser()
          : _is_opened(false), _set_mesh(false), _set_field(false),
            _set_file(false), _set_prescription(false), _set_domain(false),
            _prescription_opened(false), _nb_dof(1), _dim(2), _nb_nodes(0),
            _nb_elements(0), _nb_sides(0), _nb_edges(0), _scan(0), _dof_support(NODE_FIELD),
            _nb_mat(0), _theMesh(nullptr), _v(nullptr), _parser(nullptr), _ipf(nullptr)
{ }


XMLParser::XMLParser(string file,
                     int    type)
          : _is_opened(false), _set_mesh(true), _set_field(false), _set_file(true),
            _set_domain(false), _prescription_opened(false), _type(type),
            _file(file), _nb_dof(1), _dim(2), _nb_nodes(0), _nb_elements(0), _nb_sides(0),
            _nb_edges(0), _scan(0), _dof_support(NODE_FIELD), _nb_mat(0), _theMesh(nullptr),
            _v(nullptr), _parser(nullptr), _ipf(nullptr)
{
   open();
}


XMLParser::XMLParser(string file,
                     Mesh&  ms,
                     int    type,
                     int    format)
          : _is_opened(false), _set_mesh(true), _set_field(false), _set_file(true),
            _set_domain(false), _prescription_opened(false), _type(type),
            _format(format), _file(file), _nb_dof(1), _scan(0), _dof_support(NODE_FIELD),
            _nb_mat(0), _theMesh(&ms), _v(nullptr), _parser(nullptr), _ipf(nullptr)
{
   _nb_nodes = _theMesh->getNbNodes();
   _nb_elements = _theMesh->getNbElements();
   _nb_sides = _theMesh->getNbSides();
   _nb_edges = _theMesh->getNbEdges();
   _dim = _theMesh->getDim();
   open();
   if (type==MESH)
      get(ms);
}


XMLParser::XMLParser(const XMLParser& p)
          : _is_opened(p._is_opened), _is_closed(p._is_closed), _set_mesh(p._set_mesh),
            _set_field(p._set_field), _set_file(p._set_file), _set_domain(p._set_domain),
            _time(p._time), _sought_time(p._sought_time), _type(p._type),
            _format(p._format), _file(p._file), _mesh_file(p._mesh_file),
            _sought_name(p._sought_name), _tag_name(p._tag_name), _xml(p._xml), _mat(p._mat),
            _nb_dof(p._nb_dof), _dim(p._dim), _nb_nodes(p._nb_nodes), _nb_elements(p._nb_elements),
            _nb_sides(p._nb_sides), _nb_edges(p._nb_edges), _scan(p._scan), _nb_el_nodes(p._nb_el_nodes),
            _nb_sd_nodes(p._nb_sd_nodes), _dof_support(p._dof_support),_theMesh(p._theMesh),
            _v(p._v), _parser(p._parser), _ipf(p._ipf)
{ }


XMLParser::~XMLParser() { }


void XMLParser::set(Mesh& ms,
                    int   format)
{
   _theMesh = &ms;
   _format = format;
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


int XMLParser::scan(size_t ind)
{
   if (Verbosity>10 || _scan>1) {
      cout << "Scanning xml file: " << _file << endl;
      cout << "----------------------------------------------------------------------" << endl;
   }
   _scan = ind;
   _set_mesh = _set_field = true;
   if (parse(_xml)) {
      if (Verbosity>10 || _scan>1)
         cout << "Parse done" << endl;
      _scan = false;
      cout << "----------------------------------------------------------------------" << endl;
      cout << "Scanning complete." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::scan(size_t): Failed to parse XML file.");
   return -1;
}


int XMLParser::scan(vector<real_t>& t,
                    int             type,
                    size_t          ind)
{
   _scan = ind;
   _set_mesh = _set_prescription = false;
   _set_field = true;
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
      throw OFELIException("In XMLParser::scan(vector<real_t>,int,size_t): Failed to parse XML file.");
   return -1;
}


int XMLParser::get(Domain& dm)
{
   _theDomain = &dm;
   _set_mesh = _set_field = false;
   _set_domain = true;
   _scan = 0;
   _theMesh = nullptr;
   _type = DOMAIN_;
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
   _set_mesh = _set_field = _set_domain = false;
   _scan = 0;
   _theMesh = nullptr;
   _type = FUNCTION;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Tabulation): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(IPF& ipf)
{
   _ik1 = _ik2 = _dk1 = _dk2 = _ck = _mk = _pk = _dk = 0;
   _ipf = &ipf;
   _set_mesh = _set_field = false;
   _scan = 0;
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
   _set_mesh = _set_field = false;
   _scan = 0;
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


int XMLParser::get(int                      type,
                   vector<PrescriptionPar>& p)
{
   _prescription_type = type;
   _vp = &p;
   _vp->clear();
   _set_mesh = _set_field = false;
   _scan = 0;
   _compact = true;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done." << endl;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(int,vector<PrescriptionPar>): Failed to parse XML file.");
   return 0;
}


int XMLParser::get(Mesh& ms,
                   int   format)
{
   _set_mesh = true;
   _set_field = false;
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
   _set_field = true;
   _sought_name = name;
   _sought_time = -1.0;
   _scan = 0;
   _compact = true;
   _v = &v;
   _all_steps = 0;
   _name = _sought_name;
   _theMesh = nullptr;
   if (_v->WithMesh())
      _theMesh = &(_v->getMesh());
   _nb_dof = 1;
   _v->setName(_name);
   _nx = v.getNx(), _ny = v.getNy(), _nz = v.getNz();
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
   _set_field = true;
   _sought_time = time;
   _sought_name = name;
   _scan = 0;
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
   _set_mesh = true;
   _set_field = true;
   _sought_time = time;
   _sought_name = name;
   _scan = 0;
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
   _v->setMesh(*_theMesh,_nb_dof,_dof_support);
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
   _set_mesh = _set_field = true;
   _theMesh = &ms;
   _nb_nodes = _theMesh->getNbNodes();
   _nb_elements = _theMesh->getNbElements();
   _nb_sides = _theMesh->getNbSides();
   _nb_edges = _theMesh->getNbEdges();
   _scan = 0;
   _V = &v;
   _all_steps = 1;
   _compact = true;
   _type = FIELD;
   if (parse(_xml)) {
      if (Verbosity>10)
         cout << "Parse done" << endl;
      name = _name;
      return 0;
   }
   else
      throw OFELIException("In XMLParser::get(Mesh,vector<vector<real_t> >): Failed to parse XML file.");
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
   if (_tag_name=="Domain") {
      _set_domain = true;
      if (_scan>1)
         cout << "Found a domain description." << endl;
   }

   else if (_tag_name=="Material") {
      if (_scan>0)
         cout << "Found a material description." << endl;
   }

   else if (_tag_name=="Prescription") {
      _set_prescription = true;
      if (_scan>1)
         cout << "Found a prescription." << endl;
      else {
         _par.code = 0;
         _par.dof = 1;
         _par.bx = _par.by = _par.bz = _par.bt = false;
      }
   }

   else if (_tag_name=="Solution") {
      if (_scan>2)
         cout << "-> Prescription of exact solution" << endl;
      else
         _par.type = SOLUTION;
   }

   else if (_tag_name=="BoundaryCondition") {
      if (_scan>2)
         cout << "-> Prescription of a boundary condition" << endl;
      else
         _par.type = BOUNDARY_CONDITION;
   }

   else if (_tag_name=="BodyForce" || _tag_name=="Source") {
      if (_scan>2)
         cout << "-> Prescription of a body force" << endl;
      else
         _par.type = BODY_FORCE;
   }

   else if (_tag_name=="PointForce") {
      if (_scan>2)
         cout << "-> Prescription of a pointwise force" << endl;
      else
         _par.type = POINT_FORCE;
   }

   else if (_tag_name=="Flux" || _tag_name=="Traction"|| _tag_name=="BoundaryForce") {
      if (_scan>2)
         cout << "-> Prescription of a boundary force" << endl;
      else
         _par.type = BOUNDARY_FORCE;
   }

   else if (_tag_name=="Initial") {
      if (_scan>2)
         cout << "-> Prescription of an initial field" << endl;
      else
         _par.type = INITIAL_FIELD;
   }

   else if (_tag_name=="Function") {
      if (_scan>2)
         cout << "-> Definition of a tabulation" << endl;
   }

   else
      _par.type = 0;

   for (StringMap::iterator i=attributes.begin(); i!=attributes.end(); i++) {
      if (_ipf && _scan==0) {
         if (_type==PROJECT)
            read_project(i);
      }

      if (_set_domain && _type==DOMAIN_)
         read_domain(i);

      if (_set_mesh && _type==MESH)
         read_mesh(i);

      if (_set_field) {
         _type = FIELD;
         read_field(i);
      }

      //      if (_material)
      //         read_mat(i);

      if (_type==PRESCRIBE && _set_prescription)
         read_prescription(i);

      if (_type==FUNCTION)
         read_tab(i);
   }
   return true;
}


void XMLParser::read_project(const StringMap::iterator& i)
{
   string mat;
   if (_tag_name=="Project") {
      if (i->first=="name")
         _ipf->_project = i->second;
   }
   if (i->first=="value") {
      if (_tag_name=="verbose") {
         if (_scan>1) 
            cout << "Verbosity parameter set." << endl;
         else {
            if (i->first=="value")
               _ipf->_verbose = atoi((i->second).c_str());
         }
      }
      else if (_tag_name=="material_code") {
         if (_scan>1) 
            cout << "Material code and name given." << endl;
         else {
            if (i->first=="code")
               _code = atoi((i->second).c_str());
            if (i->first=="material" || i->first=="value")
               mat = i->second;
            theMaterial.set(_code,mat);
         }
      }
      else if (_tag_name=="output") {
         if (_scan>1) 
            cout << "Output parameter selected." << endl;
         else {
            if (i->first=="value")
               _ipf->_output = atoi((i->second).c_str());
         }
      }
      else if (_tag_name=="save") {
         if (_scan>1) 
            cout << "Save parameter selected." << endl;
         else {
            if (i->first=="value")
               _ipf->_save = atoi((i->second).c_str());
         }
      }
      else if (_tag_name=="plot") {
         if (_scan>1) 
            cout << "Frequency for plotting save parameter selected." << endl;
         else {
            if (i->first=="value")
               _ipf->_plot = atoi((i->second).c_str());
         }
      }
      else if (_tag_name=="bc") {
         if (_scan>1) 
            cout << "Boundary condition toggle selected." << endl;
         else
            _ipf->_bc = atoi((i->second).c_str());
      }
      else if (_tag_name=="bf") {
         if (_scan>1) 
            cout << "Body force toggle selected." << endl;
         else
            _ipf->_bf = atoi((i->second).c_str());
      }
      else if (_tag_name=="sf") {
         if (_scan>1) 
            cout << "Boundary force toggle selected." << endl;
         else
            _ipf->_sf = atoi((i->second).c_str());
      }
      else if (_tag_name=="init") {
         if (_scan>1) 
            cout << "Initial solution toggle selected." << endl;
         else
            _ipf->_ini = atoi((i->second).c_str());
      }
      else if (_tag_name=="prescription") {
         if (_scan>1) 
            cout << "Prescription toggle selected." << endl;
         else
            _ipf->_data = atoi((i->second).c_str());
      }
      else if (_tag_name=="nb_steps") {
         if (_scan>1) 
            cout << "Number of time steps parameter selected." << endl;
         else
            _ipf->_nb_steps = atoi((i->second).c_str());
      }
      else if (_tag_name=="nb_iter") {
         if (_scan>1) 
            cout << "Number of iterations parameter selected." << endl;
         else
            _ipf->_nb_iter = atoi((i->second).c_str());
      }
      else if (_tag_name=="time_step") {
         if (_scan>1) 
            cout << "Time step parameter selected." << endl;
         else
            _ipf->_time_step = atof((i->second).c_str());
      }
      else if (_tag_name=="max_time") {
         if (_scan>1) 
            cout << "Maximal time parameter selected." << endl;
         else
            _ipf->_max_time = atof((i->second).c_str());
      }
      else if (_tag_name=="tolerance") {
         if (_scan>1) 
            cout << "Tolerance parameter selected." << endl;
         else
            _ipf->_tolerance = atof((i->second).c_str());
      }
      else if (_tag_name=="int" || _tag_name=="integer") {
         if (_scan>1) 
            cout << "Extra integer parameter selected." << endl;
         else
            _ipf->_int_par[_ik1++] = atoi((i->second).c_str());
      }
      else if (_tag_name=="double") {
         if (_scan>1) 
            cout << "Extra double parameter selected." << endl;
         else
            _ipf->_real_par[_dk1++] = atof((i->second).c_str());
      }
      else if (_tag_name=="string") {
         if (_scan>1) 
            cout << "Extra string parameter selected." << endl;
         else
            if (i->first=="value")
               _ipf->_string_par[_dk2++] = i->second;
      }
      else if (_tag_name=="complex") {
         if (_scan>1)
            cout << "Extra complex parameter selected." << endl;
         else {
            if (i->first=="real")
               _ipf->_complex_par[_ck] += complex_t(atof((i->second).c_str()));
            if (i->first=="imag")
               _ipf->_complex_par[_ck++] += complex_t(0.,atof((i->second).c_str()));
         }
      }
      else if (_tag_name=="domain_file") {
         if (_scan>1)
            cout << "Domain file name given." << endl;
         else
            _ipf->_domain_file = i->second;
      }
      else if (_tag_name=="mesh_file") {
         if (_scan>1)
            cout << "Mesh file name given." << endl;
         else
            _ipf->_mesh_file[_mk++] = i->second;
      }
      else if (_tag_name=="init_file") {
         if (_scan>1)
            cout << "Initial field file name given." << endl;
         else
            _ipf->_init_file = i->second;
      }
      else if (_tag_name=="restart_file") {
         if (_scan>1) 
            cout << "Restart file name given." << endl;
         else
            _ipf->_restart_file = i->second;
      }
      else if (_tag_name=="bc_file") {
         if (_scan>1) 
            cout << "Boundary condition file name given." << endl;
         else
            _ipf->_bc_file = i->second;
      }
      else if (_tag_name=="bf_file") {
         if (_scan>1) 
            cout << "Body force file name given." << endl;
         else
            _ipf->_bf_file = i->second;
      }
      else if (_tag_name=="sf_file") {
         if (_scan>1) 
            cout << "Boundary force file name given." << endl;
         else
            _ipf->_sf_file = i->second;
      }
      else if (_tag_name=="save_file") {
         if (_scan>1) 
            cout << "Save file name given." << endl;
         else
            _ipf->_save_file = i->second;
      }
      else if (_tag_name=="plot_file") {
         if (_scan>1) 
            cout << "Plot file name given." << endl;
         else
            _ipf->_plot_file[_pk++] = i->second;
      }
      else if (_tag_name=="prescription_file") {
         if (_scan>1)
            cout << "Prescription file name given." << endl;
         else
            _ipf->_data_file[_dk++] = i->second;
      }
      else if (_tag_name=="aux_file") {
         if (_scan>1) 
            cout << "Auxiliary file name given." << endl;
         else
            _ipf->_aux_file[_ik2++] = i->second;
      }
      else if (_tag_name=="density") {
         if (_scan>1) 
            cout << "Density value given." << endl;
         else
            _ipf->_mp._density = i->second;
      }
      else if (_tag_name=="electric_conductivity") {
         if (_scan>1) 
            cout << "Electric conductivity value given." << endl;
         else
            _ipf->_mp._electric_cond = i->second;
      }
      else if (_tag_name=="electric_permittivity") { 
         if (_scan>1) 
            cout << "Electric permittivity value given." << endl;
         else
            _ipf->_mp._electric_perm = i->second;
      }
      else if (_tag_name=="magnetic_permeability") {
         if (_scan>1) 
            cout << "Magnetic permeability value given." << endl;
         else
            _ipf->_mp._magnetic_perm = i->second;
      }
      else if (_tag_name=="poisson_ratio") {
         if (_scan>1) 
            cout << "Poisson ratio value given." << endl;
         else
            _ipf->_mp._poisson = i->second;
      }
      else if (_tag_name=="thermal_conductivity") {
         if (_scan>1) 
            cout << "Thermal conductivity value given." << endl;
         else
            _ipf->_mp._thermal_cond = i->second;
      }
      else if (_tag_name=="rho_cp") {
         if (_scan>1) 
            cout << "Value of density x specific heat given." << endl;
         else
            _ipf->_mp._rho_cp = i->second;
      }
      else if (_tag_name=="viscosity") {
         if (_scan>1) 
            cout << "Viscosity value given." << endl;
         else
            _ipf->_mp._visc = i->second;
      }
      else if (_tag_name=="young_modulus") {
         if (_scan>1) 
            cout << "Young modulus value given." << endl;
         else
            _ipf->_mp._young = i->second;
      }
      else
         ;
   }

   if (_tag_name=="parameter" || _tag_name=="param") {
      if (_scan>1)
         cout << "keyword parameter selected." << endl;
      else {
         if (i->first=="label")
            _ipf->_param_label.push_back(i->second);
         if (i->first=="value")
            _ipf->_param_value.push_back(i->second);
         if (i->first=="ext" || i->first=="extension")
            _ipf->_param_ext.push_back(i->second);
      }
   }

   if (_tag_name=="array") {
      if (_scan>1)
         cout << "keyword array selected." << endl;
      else {
         if (i->first=="label")
            _ipf->_array_label.push_back(i->second);
         if (i->first=="size")
            _ipf->_array_size.push_back(atoi((i->second).c_str()));
         if (i->first=="ext" || i->first=="extension")
            _ipf->_array_ext.push_back(i->second);
      }
   }
}


void XMLParser::read_domain(const StringMap::iterator& i)
{
   if (_tag_name=="Domain") {
      if (i->first=="dim") {
         _dim = atoi((i->second).c_str());
         if (_scan==false)
            _theDomain->_dim = _dim;
      }
      else if (i->first=="nb_dof") {
         _nb_dof = atoi((i->second).c_str());
         if (_scan==false)
            _theDomain->_nb_dof = _nb_dof;
      }
      if (_scan)
         cout << "Found a domain:" << endl;
   }

   else if (_tag_name=="vertex") {
      real_t x[3]={0,0,0}, h;
      int code=0;
      if (i->first=="x")
         x[0] = atof((i->second).c_str());
      if (i->first=="y")
         x[1] = atof((i->second).c_str());
      if (i->first=="z")
         x[2] = atof((i->second).c_str());
      if (i->first=="code")
         code = atoi((i->second).c_str());
      if (i->first=="h") {
         h = atof((i->second).c_str());
         _theDomain->insertVertex(x[0],x[1],h,code);
      }
   }

   else if (_tag_name=="line") {
   }

   else if (_tag_name=="contour") {
   }

   else if (_tag_name=="hole") {
   }

   else if (_tag_name=="subdomain") {
   }

   else
      ;
}


void XMLParser::read_tab(const StringMap::iterator& i)
{
   if (_tag_name=="Function") {
      if (_scan)
         cout << "Found a tabulated function:" << endl;
      if (i->first=="name") {
         _theTabulation->setFunction(i->second);
         _nb_funct = _theTabulation->_nb_funct;
         _nb_var = 0;
      }
   }

   else if (_tag_name=="Variable") {
      if (i->first=="label") {
         _theTabulation->setVariable(i->second);
         _nb_var++; 
      }
      else if (i->first=="nb_pts") {
         size_t k = atoi((i->second).c_str());
         if (_scan==false)
            _theTabulation->_funct(_nb_funct).Np[_nb_var-1] = k;
      }
      else if (i->first=="min")
         _theTabulation->_funct(_nb_funct).Min[_nb_var-1] = atof((i->second).c_str());    
      else if (i->first=="max")
         _theTabulation->_funct(_nb_funct).Max[_nb_var-1] = atof((i->second).c_str());    
   }

   else if (_tag_name=="Data") {
   }

   else
      ;
}


void XMLParser::read_mat(const StringMap::iterator& i)
{
   if (_scan)
      return;
   if (_tag_name=="Material") {
      if (i->first=="name") {
         theMaterial._name = i->second;
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
   }
   if (_tag_name=="Density") {
      if (_scan>1)
         cout << "   Found density value." << endl;
      if (_scan==0) {
         theMaterial._density[_nb_mat].exist = true;
         theMaterial._density[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._density[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="SpecificHeat") {
      if (_scan>1)
         cout << "   Found specific heat value." << endl;
      if (_scan==0) {
         theMaterial._specific_heat[_nb_mat].exist = true;
         theMaterial._specific_heat[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._specific_heat[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="ThermalConductivity") {
      if (_scan>1)
         cout << "   Found thermal conductivity value." << endl;
      if (_scan==0) {
         theMaterial._thermal_conductivity[_nb_mat].exist = true;
         theMaterial._thermal_conductivity[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._thermal_conductivity[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="MeltingTemperature") {
      if (_scan>1)
         cout << "   Found melting temperature value." << endl;
      if (_scan==0) {
         theMaterial._melting_temperature[_nb_mat].exist = true;
         theMaterial._melting_temperature[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._melting_temperature[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="EvaporationTemperature") {
      if (_scan>1)
         cout << "   Found evaporation temperature value." << endl;
      if (_scan==0) {
         theMaterial._evaporation_temperature[_nb_mat].exist = true;
         theMaterial._evaporation_temperature[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._evaporation_temperature[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="ThermalExpansion") {
      if (_scan>1)
         cout << "   Found thermal expansion value." << endl;
      if (_scan==0) {
         theMaterial._thermal_expansion[_nb_mat].exist = true;
         theMaterial._thermal_expansion[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._thermal_expansion[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="LatentHeatMelting") {
      if (_scan>1)
         cout << "   Found latent heat for melting value." << endl;
      if (_scan==0) {
         theMaterial._latent_heat_melting[_nb_mat].exist = true;
         theMaterial._latent_heat_melting[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._latent_heat_melting[_nb_mat-1].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="LatentHeatEvaporation") {
      if (_scan>1)
         cout << "   Found latent heat for evaporation value." << endl;
      if (_scan==0) {
         theMaterial._latent_heat_evaporation[_nb_mat].exist = true;
         theMaterial._latent_heat_evaporation[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._latent_heat_evaporation[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="DielectricConstant") {
      if (_scan>1)
         cout << "   Found dielectric constant value." << endl;
      if (_scan==0) {
         theMaterial._dielectric_constant[_nb_mat].exist = true;
         theMaterial._dielectric_constant[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._dielectric_constant[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="ElectricConductivity") {
      if (_scan>1)
         cout << "   Found electric conductivity value." << endl;
      if (_scan==0) {
         theMaterial._electric_conductivity[_nb_mat].exist = true;
         theMaterial._electric_conductivity[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._electric_conductivity[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="ElectricResistivity") {
      if (_scan>1)
         cout << "   Found electric resistivity value." << endl;
      if (_scan==0) {
         theMaterial._electric_resistivity[_nb_mat].exist = true;
         theMaterial._electric_resistivity[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._electric_resistivity[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="MagneticPermeability") {
      if (_scan>1)
         cout << "   Found magnetic permeability value." << endl;
      if (_scan==0) {
         theMaterial._magnetic_permeability[_nb_mat].exist = true;
         theMaterial._magnetic_permeability[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._magnetic_permeability[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="Viscosity") {
      if (_scan>1)
         cout << "   Found viscosity value." << endl;
      if (_scan==0) {
         theMaterial._viscosity[_nb_mat].exist = true;
         theMaterial._viscosity[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._viscosity[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="YoungModulus") {
      if (_scan>1)
         cout << "   Found Young modulus value." << endl;
      if (_scan==0) {
         theMaterial._young_modulus[_nb_mat].exist = true;
         theMaterial._young_modulus[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._young_modulus[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else if (_tag_name=="PoissonRatio") {
      if (_scan>1)
         cout << "   Found Poisson ratio value." << endl;
      if (_scan==0) {
         theMaterial._poisson_ratio[_nb_mat].exist = true;
         theMaterial._poisson_ratio[_nb_mat].type = BY_VALUE;
         if (i->first=="value")
            theMaterial._poisson_ratio[_nb_mat].value = atof((i->second).c_str());
      }
   }
   else
      ;
}


void XMLParser::read_mesh(const StringMap::iterator& i)
{
   if (_tag_name=="Mesh") {
      if (_dim==1) {
         _el_shape = "line";
         _nb_el_nodes = 2;
      }
      else if (_dim==2) {
         _el_shape = "triangle";
         _nb_el_nodes = 3;
         _sd_shape = "line";
         _nb_sd_nodes = 2;
      }
      else if (_dim==3) {
         _el_shape = "tetrahedron";
         _nb_el_nodes = 4;
         _sd_shape = "triangle";
         _nb_sd_nodes = 3;
      }
      else
         ;
      if (i->first=="file" && _scan==0)
         new XMLParser(i->second,*_theMesh);
      else if (i->first=="dim") {
         _dim = atoi((i->second).c_str());
         if (_scan==false)
            _theMesh->setDim(_dim);
      }
      else if (i->first=="nb_dof")
         _nb_dof = atoi((i->second).c_str());
      if (_scan)
         cout << "Found a finite element mesh:" << endl;
   }

   else if (_tag_name=="Elements") {
      if (i->first=="shape")
         _el_shape = i->second;
      if (i->first=="nodes")
         _nb_el_nodes = atoi((i->second).c_str());
   }

   else if (_tag_name=="Sides") {
      if (i->first=="shape")
         _sd_shape = i->second;
      if (i->first=="nodes")
         _nb_sd_nodes = atoi((i->second).c_str());
   }

   else if (_tag_name=="Material") {
   }

   else
      ;
}


void XMLParser::read_prescription(const StringMap::iterator& i)
{
   vector<string> bccs(6);
   bccs[0] = "PERIODIC_A"; bccs[1] = "PERIODIC_B"; bccs[2] = "CONTACT";
   bccs[3] = "CONTACT_M"; bccs[4] = "CONTACT_S"; bccs[5] = "SLIP";
   int bcc[6] = {9999,-9999,9998,9997,-9997,9996};
 
   if (i->first=="file" && _scan==0)
      new XMLParser(i->second,*_theMesh);
   if (_scan)
      cout << "Found a prescription:" << endl;

   if (i->first=="code")
      _par.code = BoundaryConditionCode(bccs,bcc,i->second);
   else if (i->first=="dof")
      _par.dof = atoi((i->second).c_str());
   else if (i->first=="x")
      _par.x = atof((i->second).c_str()), _par.bx = true;
   else if (i->first=="y")
      _par.y = atof((i->second).c_str()), _par.by = true;
   else if (i->first=="z")
      _par.z = atof((i->second).c_str()), _par.bz = true;
   else if (i->first=="time")
      _par.t = atof((i->second).c_str()), _par.bt = true;
   else if (i->first=="value") {
   }
   else
      ;
}


void XMLParser::read_field(const StringMap::iterator& i)
{
   if (_tag_name=="Field" && (_type==MESH||_type==FIELD)) {
      if (i->first=="file" && _scan==0) {
         _parser = new XMLParser(i->second);
         _parser->set(*_theMesh);
         _parser->_v = _v;
      }
      else if (i->first=="name")
         _name = i->second;
      else if (i->first=="type") {
         if (i->second=="Node")
            _dof_support = NODE_FIELD;
         else if (i->second=="Element")
            _dof_support = ELEMENT_FIELD;
         else if (i->second=="Side")
            _dof_support = SIDE_FIELD;
         else if (i->second=="Edge")
            _dof_support = EDGE_FIELD;
         else
            ;
      }
      else if (i->first=="nb_dof")
         _nb_dof = atoi((i->second).c_str());
      if (i->first=="nx")
         _nx = atoi((i->second).c_str());
      if (i->first=="ny")
         _ny = atoi((i->second).c_str());
      if (i->first=="nz")
         _nz = atoi((i->second).c_str());
   }
   if (_tag_name=="Step") {
      _compact = false;
      if (i->first=="time")
         _time = atof((i->second).c_str());
   }
   else if (_tag_name=="constant" || _tag_name=="const") {
      _dof = 0;
      if (i->first=="value")
         _val = atof((i->second).c_str());
      else if (i->first=="dof")
         _dof = atoi((i->second).c_str());
   }
   else if (_tag_name=="expression" || _tag_name=="expr") {
      _dof = 0;
      if (i->first=="dof")
         _dof = atoi((i->second).c_str());
   }
}


bool XMLParser::on_cdata(string cdata)
{
   vector<string> tokens;
   {
      string buf;
      std::stringstream ss(cdata);
      while (ss >> buf)
         tokens.push_back(buf);
   }
   vector<string>::iterator it=tokens.begin();

   if (_type==PROJECT)
      read_project_data(tokens,it);

   else if (_type==MESH)
      read_mesh_data(tokens,it);

   else if (_type==DOMAIN_)
      read_domain_data(tokens,it);

// Material in mesh data
   else if (_tag_name=="Material" && _set_mesh) {
      while (it!=tokens.end()) {
         _cm = atoi((*it++).c_str());
         _mat = *it++;
         _nb_mat = theMaterial.set(_cm,_mat);
      }
   }

// Material
   else if (_type==MATERIAL)
      read_mat_data(tokens,it);

// Prescribe
   else if (_type==PRESCRIBE && _set_prescription && _prescription_type==_par.type)
      read_prescribe_data(tokens,it);

// Field
   else if ((_tag_name=="Step"||(_tag_name=="Field"&&_compact)) && _set_field)
      read_field_data(tokens,it);
   else if ((_tag_name=="constant"||_tag_name=="const") && _set_field)
      read_const_field_data(tokens,it);
   else if ((_tag_name=="expression"||_tag_name=="expr") && _set_field)
      read_exp_field_data(tokens,it);

// Tabulated functions
   else if (_type==FUNCTION) {
      if (_tag_name=="Data")
         read_tab_data(tokens,it);
   }

   else
      ;
   return true;
}


void XMLParser::read_prescribe_data(const vector<string>&     tokens,
                                    vector<string>::iterator& it)
{
   if (it!=tokens.end() && _scan==0) {
      _par.fct = *it++;
      _vp->push_back(_par);
   }
}


void XMLParser::read_field_data(const vector<string>&     tokens,
                                vector<string>::iterator& it)
{
   size_t nb = 0;
   if (_theMesh) {
      if (_dof_support==NODE_FIELD)
         nb = _nb_nodes*_nb_dof;
      else if (_dof_support==ELEMENT_FIELD)
         nb = _nb_elements*_nb_dof;
      else if (_dof_support==SIDE_FIELD)
         nb = _nb_sides*_nb_dof;
      else if (_dof_support==EDGE_FIELD)
         nb = _nb_edges*_nb_dof;
      else
         nb = _nx*_ny*_nz;
   }

   if (_scan) {
      if (_compact) {
         while (it!=tokens.end()) {
            _time = atof((*it++).c_str());
            _ft->push_back(_time);
            for (size_t n=0; n<nb; n++)
               *it++;
         }
      }
      else
         _ft->push_back(_time);
      if (_dof_support==NODE_FIELD)
         cout << "Found a nodewise field vector, Name: " << _name << ", Time = " 
              << _time << ", Nb. of DOF: " << _nb_dof << endl;
      else if (_dof_support==ELEMENT_FIELD)
         cout << "Found an elementwise field vector, Name: " << _name << ", Time = "
              << _time << ", Nb. of DOF: " << _nb_dof << endl;
      else if (_dof_support==SIDE_FIELD)
         cout << "Found a sidewise field vector, Name: " << _name << ", Time = "
              << _time << ", Nb. of DOF: " << _nb_dof << endl;
      else if (_dof_support==EDGE_FIELD)
         cout << "Found an edgewise field vector, Name: " << _name << ", Time = "
              << _time << ", Nb. of DOF: " << _nb_dof << endl;
      else
         ;
   }
   else {
      if (_all_steps>0 && _compact) {
         size_t i = 0;
         while (it!=tokens.end()) {
            _time = atof((*it++).c_str());
            for (size_t j=0; j<nb; j++)
               (*_V)[i].push_back(atof((*it++).c_str()));
            i++;
         }
      }
      else {
         if (_compact) {
            if ((_name==_sought_name || _sought_name=="ANYTHING") && _set_field) {
               if (_nx > 0) {
                  _v->setSize(_nx,_ny,_nz);
                  while (it!=tokens.end()) {
                     if (_sought_time >= 0)
                        _time = atof((*it++).c_str());
                     for (size_t i=1; i<=_nx; i++) {
                        for (size_t j=1; j<=_ny; j++) {
                           for (size_t k=1; k<=_nz; k++)
                              (*_v)(i,j,k) = atof((*it++).c_str());
                        }
                     }
                  }
               }
               else {
                  if (nb > 0) {
                     while (it!=tokens.end()) {
                        _time = atof((*it++).c_str());
                        if (_time==_sought_time || _sought_time==-1.0) {
                           _v->setMesh(*_theMesh,_nb_dof,_dof_support);
                           _v->setName(_name);
                           _v->setTime(_time);
                           for (size_t n=1; n<=_v->getNb(); n++)
                              for (size_t k=1; k<=_nb_dof; k++)
                                 (*_v)(n,k) = atof((*it++).c_str());
                        }
                        else {
                           for (size_t n=1; n<=_v->getNb()*_nb_dof; n++)
                              *it++;
                        }
                     }
                  }
                  else {
                     vector<real_t> V;
                     _time = atof((*it++).c_str());
                     while (it!=tokens.end())
                        V.push_back(atof((*it++).c_str()));
                     _v->setSize(V.size());
                     for (size_t i=0; i<V.size(); i++)
                        (*_v)[i] = V[i];
                  }
               }
            }
         }
         else {
            if ((_name==_sought_name || _sought_name=="ANYTHING") &&
                (_time==_sought_time || _sought_time==-1.0) && _set_field) {
               _v->setMesh(*_theMesh,_nb_dof,_dof_support);
               _v->setName(_name);
               _v->setTime(_time);
               for (size_t n=1; n<=_v->getNb(); n++)
                  for (size_t k=1; k<=_nb_dof; k++)
                     (*_v)(n,k) = atof((*it++).c_str());
            }
         }
      }
   }
}


void XMLParser::read_const_field_data(const vector<string>&     tokens,
                                      vector<string>::iterator& it)
{
   if (_scan) {
      _ft->push_back(_time);
      if (_scan>1) {
         if (_dof_support==NODE_FIELD)
            cout << "Found a nodewise field vector, Name: " << _name << ", Time = " 
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==ELEMENT_FIELD)
            cout << "Found an elementwise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==SIDE_FIELD)
            cout << "Found a sidewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==EDGE_FIELD)
            cout << "Found an edgewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else
            ;
      }
   }
   else {
      if ((_name==_sought_name || _sought_name=="ANYTHING") &&
          (_time==_sought_time || _sought_time==-1.0) && _set_field) {
         _v->setTime(_time);
         _val = atof((*it++).c_str());
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


void XMLParser::read_const_field_data()
{
   vector<string> tokens;
   vector<string>::iterator it;
   if (_scan) {
      _ft->push_back(_time);
      if (_scan>1) {
         if (_dof_support==NODE_FIELD)
            cout << "Found a nodewise field vector, Name: " << _name << ", Time = " 
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==ELEMENT_FIELD)
            cout << "Found an elementwise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==SIDE_FIELD)
            cout << "Found a sidewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==EDGE_FIELD)
            cout << "Found an edgewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else
            ;
      }
   }
   else {
      if ((_name==_sought_name || _sought_name=="ANYTHING") &&
          (_time==_sought_time || _sought_time==-1.0) && _set_field) {
         _v->setTime(_time);
         _val = atof((*it++).c_str());
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
}


void XMLParser::read_exp_field_data(const vector<string>&     tokens,
                                    vector<string>::iterator& it)
{
   if (_scan) {
      _ft->push_back(_time);
      if (_scan>1) {
         if (_dof_support==NODE_FIELD)
            cout << "Found a nodewise field vector, Name: " << _name << ", Time = " 
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==ELEMENT_FIELD)
            cout << "Found an elementwise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==SIDE_FIELD)
            cout << "Found a sidewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else if (_dof_support==EDGE_FIELD)
            cout << "Found an edgewise field vector, Name: " << _name << ", Time = "
                 << _time << ", Nb. of DOF: " << _nb_dof << endl;
         else
            ;
      }
   }
   else {
      if ((_name==_sought_name || _sought_name=="ANYTHING") &&
          (_time==_sought_time || _sought_time==-1.0) && _set_field) {
         _v->setTime(_time);
         theParser.Parse((*it++).c_str(),"x,y,z,t");	    
         for (size_t n=1; n<=_v->getNb(); n++) {
            if (_dof==0)
               for (size_t k=1; k<=_v->getNbDOF(); k++)
                  parse_exp(n,k);
            else
               parse_exp(n,_dof);
         }
      }
   }
}


void XMLParser::read_tab_data(const vector<string>&     tokens,
                              vector<string>::iterator& it)
{
   _theTabulation->setSizes();
   size_t i=1;
   while (it!=tokens.end())
      _theTabulation->_funct(_nb_funct).Val(i++) = atof((*it++).c_str());
}


void XMLParser::parse_exp(size_t n,
                          size_t k)
{
   Point<real_t> c;
   if (_dof_support==NODE_FIELD)
      c = _theMesh->getPtrNode(n)->getCoord();
   else if (_dof_support==SIDE_FIELD) {
      theSide = _theMesh->getPtrSide(n);
      if (theSide->getShape()==LINE)
	 c = Line2(theSide).getCenter();
      else if (theSide->getShape()==TRIANGLE)
         c = Triang3(theSide).getCenter();
      else if (theSide->getShape()==QUADRILATERAL)
         c = Quad4(theSide).getCenter();
   }
   else if (_dof_support==ELEMENT_FIELD) {
      theElement = _theMesh->getPtrElement(n);
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
   real_t z=theParser.Eval(c,_time);
   (*_v)(n,k) = z;
}


void XMLParser::read_mat_data(const vector<string>&     tokens,
                              vector<string>::iterator& it)
{
   if (_tag_name=="Density") {
      if (_scan)
         cout << "   Density is given." << endl;
      else {
         theMaterial._density[_nb_mat].exist = true;
         theMaterial._density[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._density[_nb_mat].value = atoi(it->c_str());
      }
   }
   else if (_tag_name=="SpecificHeat") {
      if (_scan)
         cout << "   Specific heat is given." << endl;
      else {
         theMaterial._specific_heat[_nb_mat].exist = true;
         theMaterial._specific_heat[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._specific_heat[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="ThermalConductivity") {
      if (_scan)
         cout << "   Thermal conductivity is given." << endl;
      else {
         theMaterial._thermal_conductivity[_nb_mat].exist = true;
         theMaterial._thermal_conductivity[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._thermal_conductivity[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="MeltingTemperature") {
      if (_scan)
         cout << "   Melting temperature is given." << endl;
      else {
         theMaterial._melting_temperature[_nb_mat].exist = true;
         theMaterial._melting_temperature[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._melting_temperature[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="EvaporationTemperature") {
      if (_scan)
         cout << "   Evaporation temperature is given." << endl;
      else {
         theMaterial._evaporation_temperature[_nb_mat].exist = true;
         theMaterial._evaporation_temperature[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._evaporation_temperature[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="ThermalExpansion") {
      if (_scan)
         cout << "   Thermal expansion is given." << endl;
      else {
         theMaterial._thermal_expansion[_nb_mat].exist = true;
         theMaterial._thermal_expansion[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._thermal_expansion[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="LatentHeatMelting") {
      if (_scan)
         cout << "   Latent heat for melting is given." << endl;
      else {
         theMaterial._latent_heat_melting[_nb_mat].exist = true;
         theMaterial._latent_heat_melting[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._latent_heat_melting[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="LatentHeatEvaporation") {
      if (_scan)
         cout << "   Latent heat for evaporation is given." << endl;
      else {
         theMaterial._latent_heat_evaporation[_nb_mat].exist = true;
         theMaterial._latent_heat_evaporation[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._latent_heat_evaporation[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="DielectricConstant") {
      if (_scan)
         cout << "   Dielectric constant is given." << endl;
      else {
         theMaterial._dielectric_constant[_nb_mat].exist = true;
         theMaterial._dielectric_constant[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._dielectric_constant[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="ElectricConductivity") {
      if (_scan)
         cout << "   Electric conductivity is given." << endl;
      else {
         theMaterial._electric_conductivity[_nb_mat].exist = true;
         theMaterial._electric_conductivity[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._electric_conductivity[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="ElectricResistivity") {
      if (_scan)
         cout << "   Electric resistivity is given." << endl;
      else {
         theMaterial._electric_resistivity[_nb_mat].exist = true;
         theMaterial._electric_resistivity[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._electric_resistivity[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="MagneticPermeability") {
      if (_scan)
         cout << "   Magnetic permeability is given." << endl;
      else {
         theMaterial._magnetic_permeability[_nb_mat].exist = true;
         theMaterial._magnetic_permeability[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._magnetic_permeability[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="Viscosity") {
      if (_scan)
         cout << "   Viscosity is given." << endl;
      else {
         theMaterial._viscosity[_nb_mat].exist = true;
         theMaterial._viscosity[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._viscosity[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="YoungModulus") {
      if (_scan)
         cout << "   Young modulus is given." << endl;
      else {
         theMaterial._young_modulus[_nb_mat].exist = true;
         theMaterial._young_modulus[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._young_modulus[_nb_mat].value = atof(it->c_str());
      }
   }
   else if (_tag_name=="PoissonRatio") {
      if (_scan)
         cout << "   Poisson ratio is given." << endl;
      else {
         theMaterial._poisson_ratio[_nb_mat].exist = true;
         theMaterial._poisson_ratio[_nb_mat].type = BY_VALUE;
         if (it!=tokens.end())
            theMaterial._poisson_ratio[_nb_mat].value = atof(it->c_str());
      }
   }
   else
      ;
}


void XMLParser::read_project_data(const vector<string>&     tokens,
                                  vector<string>::iterator& it)
{
   if (_tag_name=="verbose") {
      if (it!=tokens.end())
         _ipf->_verbose = atoi(it->c_str());
   }
   else if (_tag_name=="material_code") {
      if (it!=tokens.end())
         theMaterial.set(_code,it->c_str());
   }
   else if (_tag_name=="output") {
      if (it!=tokens.end())
         _ipf->_output = atoi(it->c_str());
   }
   else if (_tag_name=="save") {
      if (it!=tokens.end())
         _ipf->_save = atoi(it->c_str());
   }
   else if (_tag_name=="plot") {
      if (it!=tokens.end())
         _ipf->_plot = atoi(it->c_str());
   }
   else if (_tag_name=="bc") {
      if (it!=tokens.end())
         _ipf->_bc = atoi(it->c_str());
   }
   else if (_tag_name=="bf") {
      if (it!=tokens.end())
         _ipf->_bf = atoi(it->c_str());
   }
   else if (_tag_name=="sf") {
      if (it!=tokens.end())
         _ipf->_sf = atoi(it->c_str());
   }
   else if (_tag_name=="init") {
      if (it!=tokens.end())
         _ipf->_ini = atoi(it->c_str());
   }
   else if (_tag_name=="data") {
      if (it!=tokens.end())
         _ipf->_data = atoi(it->c_str());
   }
   else if (_tag_name=="nb_steps") {
      if (it!=tokens.end())
         _ipf->_nb_steps = atoi(it->c_str());
   }
   else if (_tag_name=="nb_iter") {
      if (it!=tokens.end())
         _ipf->_nb_iter = atoi(it->c_str());
   }
   else if (_tag_name=="time_step") {
      if (it!=tokens.end())
         _ipf->_time_step = atof(it->c_str());
   }
   else if (_tag_name=="max_time") {
      if (it!=tokens.end())
         _ipf->_max_time = atof(it->c_str());
   }
   else if (_tag_name=="tolerance") {
      if (it!=tokens.end())
         _ipf->_tolerance = atof(it->c_str());
   }
   else if (_tag_name=="int") {
      if (it!=tokens.end())
         _ipf->_int_par[_ik1++] = atoi(it->c_str());
   }
   else if (_tag_name=="double") {
      if (it!=tokens.end())
         _ipf->_real_par[_dk1++] = atoi(it->c_str());
   }
   else if (_tag_name=="string") {
      if (it!=tokens.end())
         _ipf->_string_par[_dk2++] = *it;
   }
   else if (_tag_name=="complex") {
      while (it!=tokens.end()) {
         _ipf->_complex_par[_ck]  = complex_t(atof((*it++).c_str()));
         _ipf->_complex_par[_ck] += complex_t(0.,atof((*it++).c_str()));
      }
      _ck++;
   }
   else if (_tag_name=="parameter" || _tag_name=="param") {
      if (it!=tokens.end())
         _ipf->_param_value.push_back(*it);
   }
   else if (_tag_name=="domain_file") {
      if (it!=tokens.end())
         _ipf->_domain_file = *it;
   } 
   else if (_tag_name=="mesh_file") {
      if (it!=tokens.end())
         _ipf->_mesh_file[_mk++] = *it;
   } 
   else if (_tag_name=="init_file") {
      if (it!=tokens.end())
         _ipf->_init_file = *it;
   } 
   else if (_tag_name=="restart_file") {
      if (it!=tokens.end())
         _ipf->_restart_file = *it;
   }
   else if (_tag_name=="bc_file") {
      if (it!=tokens.end())
         _ipf->_bc_file = *it;
   } 
   else if (_tag_name=="bf_file") {
      if (it!=tokens.end())
         _ipf->_bf_file = *it;
   }
   else if (_tag_name=="sf_file") {
      if (it!=tokens.end())
         _ipf->_sf_file = *it;
   }
   else if (_tag_name=="save_file") {
      if (it!=tokens.end())
         _ipf->_save_file = *it;
   }
   else if (_tag_name=="plot_file") {
      if (it!=tokens.end())
         _ipf->_plot_file[_pk++] = *it;
   } 
   else if (_tag_name=="prescription_file") {
      if (it!=tokens.end())
         _ipf->_data_file[_dk++] = *it;
   } 
   else if (_tag_name=="aux_file") {
      if (it!=tokens.end())
         _ipf->_aux_file[_ik2++] = *it;
   } 
   else if (_tag_name=="density") {
      if (it!=tokens.end())
         _ipf->_mp._density = *it;
   } 
   else if (_tag_name=="electric_conductivity") {
      if (it!=tokens.end())
         _ipf->_mp._electric_cond = *it;
   } 
   else if (_tag_name=="electric_permittivity") {
      if (it!=tokens.end())
         _ipf->_mp._electric_perm = *it;
   } 
   else if (_tag_name=="magnetic_permeability") {
      if (it!=tokens.end())
         _ipf->_mp._magnetic_perm = *it;
   } 
   else if (_tag_name=="poisson_ratio") {
      if (it!=tokens.end())
         _ipf->_mp._poisson = *it;
   } 
   else if (_tag_name=="thermal_conductivity") {
      if (it!=tokens.end())
         _ipf->_mp._thermal_cond = *it;
   } 
   else if (_tag_name=="rho_cp") {
      if (it!=tokens.end())
         _ipf->_mp._rho_cp = *it;
   } 
   else if (_tag_name=="viscosity") {
      if (it!=tokens.end())
         _ipf->_mp._visc = *it;
   }
   else if (_tag_name=="young_modulus") {
      if (it!=tokens.end())
         _ipf->_mp._young = *it;
   }
   else if (_tag_name=="array") {
      size_t k = 0;
      while (it!=tokens.end() && _scan==0)
         (_ipf->_array_value[k++]).push_back(atof((*it++).c_str()));
      _ipf->_array_size.push_back(k);
   }
   else ;
}


void XMLParser::read_domain_data(const vector<string>&     tokens,
                                 vector<string>::iterator& it)
{
   vector<size_t> c;
   if (_tag_name=="vertex" && _set_domain) {
      while (it!=tokens.end()) {
         real_t x[3];
         for (size_t j=0; j<_dim; j++)
            x[j] = atof((*it++).c_str());
         int c1 = atoi((*it++).c_str());
         real_t h = atof((*it++).c_str());
         _theDomain->insertVertex(x[0],x[1],h,c1);
      }
   }
   else if (_tag_name=="line" && _set_domain) {
      while (it!=tokens.end()) {
         size_t v1 = atoi((*it++).c_str());
         size_t v2 = atoi((*it++).c_str());
         int dc = atoi((*it++).c_str()), nc = dc;
         _theDomain->insertLine(v1,v2,dc,nc);
      }
   }
   else if (_tag_name=="circle" && _set_domain) {
      while (it!=tokens.end()) {
         size_t n1 = atoi((*it++).c_str());
         size_t n2 = atoi((*it++).c_str());
         size_t n3 = atoi((*it++).c_str());
         int dc = atoi((*it++).c_str()), nc = dc;
         _theDomain->insertCircle(n1,n2,n3,dc,nc);
      }
   }
   else if (_tag_name=="contour" && _set_domain) {
      c.clear();
      while (it!=tokens.end())
         c.push_back(atoi((*it++).c_str()));
      _theDomain->insertContour(c);
   }
   else if (_tag_name=="hole" && _set_domain) {
      c.clear();
      while (it!=tokens.end())
         c.push_back(atoi((*it++).c_str()));
      _theDomain->insertHole(c);
   }
   else if (_tag_name=="RequiredVertex" && _set_domain) {
      while (it!=tokens.end())
         _theDomain->insertRequiredVertex(atoi((*it++).c_str()));
   }
   else if (_tag_name=="RequiredEdge" && _set_domain) {
      while (it!=tokens.end())
         _theDomain->insertRequiredEdge(atoi((*it++).c_str()));
   }
   else if ((_tag_name=="subdomain"||_tag_name=="SubDomain") && _set_domain) {
      while (it!=tokens.end()) {
         size_t ln = atoi((*it++).c_str());
         int orient = atoi((*it++).c_str());
         int ref = atoi((*it++).c_str());
         _theDomain->insertSubDomain(ln,orient,ref);
      }
   }
   else if (_tag_name=="rectangle" && _set_domain) {
      real_t x[4];
      while (it!=tokens.end()) {
         x[0] = atof((*it++).c_str());
         x[1] = atof((*it++).c_str());
         x[2] = atof((*it++).c_str());
         x[3] = atof((*it++).c_str());
         size_t n1 = atoi((*it++).c_str());
         size_t n2 = atoi((*it++).c_str());
         real_t r = atof((*it++).c_str());
         int c1 = atoi((*it++).c_str());
         int c2 = atoi((*it++).c_str());
         int c3 = atoi((*it++).c_str());
         int c4 = atoi((*it++).c_str());
         string file = *it++;
         _theDomain->Rectangle(x,n1,n2,r,c1,c2,c3,c4);
         _theMesh->put(file);
      }
   }
}


void XMLParser::read_mesh_data(const vector<string>&     tokens,
                               vector<string>::iterator& it)
{
   int code[MAX_NBDOF_NODE];
   size_t first_dof = 1;
   if (_tag_name=="Nodes" || _tag_name=="Sides")
      first_dof = 1;

// Material
   if (_tag_name=="Material") {
      while (it!=tokens.end()) {
         int n = atoi((*it++).c_str());
         string name = *it++;
         _nb_mat = theMaterial.set(n,name);
      }
   }

// Nodes
   if (_tag_name=="Nodes" && _set_mesh) {
      _nb_nodes = 0;
      while (it!=tokens.end()) {
         Point<real_t> a;
         a.x = atof((*it++).c_str());
         if (_dim > 1)
            a.y = atof((*it++).c_str());
         if (_dim > 2)
            a.z = atof((*it++).c_str());
         int mark = atoi((*it++).c_str());
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
      if (_scan)
         cout << setw(8) << _nb_nodes << " nodes" << endl;
   }

// Elements
   if (_tag_name=="Elements" && _set_mesh) {
      _nb_elements = 0;
      while (it!=tokens.end()) {
         int nnd[20];
         for (size_t j=0; j<_nb_el_nodes; j++)
            nnd[j] = atoi((*it++).c_str());
         int code = atoi((*it++).c_str());
         _nb_elements++;
         if (_scan==0) {
            Element *el = new Element(_nb_elements,_el_shape,code);
            for (size_t k=0; k<_nb_el_nodes; k++)
               el->Add((*_theMesh)[nnd[k]]);
            _theMesh->Add(el);
         }
      }
      if (_scan)
         cout << setw(8) << _nb_elements << " elements" << endl;
   }

// Sides
   if (_tag_name=="Sides" && _set_mesh) {
      _nb_sides = 0;
      while (it!=tokens.end()) {
         int nnd[8];
         for (size_t j=0; j<_nb_sd_nodes; j++)
            nnd[j] = atoi((*it++).c_str());
         int mark = atoi((*it++).c_str());
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
      if (_scan)
         cout << setw(8) << _nb_sides << " sides" << endl;
   }
}


bool XMLParser::on_tag_close(string tag_name)
{
   if (Verbosity>10)
      cout << "CLOSING TAG: " << tag_name << endl;
   if (tag_name=="Constant" || tag_name=="Const") {
      read_const_field_data();
   }
   if (tag_name=="Prescription")
      _set_prescription = false;
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
