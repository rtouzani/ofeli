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

           Definition of class 'XMLParser' for parsing OFELI XML files

  ==============================================================================*/


#ifndef __XML_PARSER_H
#define __XML_PARSER_H

#include <fstream>
#include <string>
#include <vector>
using std::string;
using std::vector;

#include "OFELI_Config.h"
#include "io/xmlsp/xmlsp.h"
#include "io/Prescription.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/DMatrix.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

typedef std::pair<std::string,std::string> SString;

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

 enum {
    DOMAIN_,
    MESH_,
    GRID,
    VECTOR,
    MATRIX,
    PROJECT,
    INPUT,
    MATERIAL,
    PRESCRIBE,
    FUNCTION
 };
 enum { ASCII, BINARY };

   XMLParser();
   XMLParser(string file, int type=MESH_);
   XMLParser(string file, Mesh& ms, int type=MESH_, int format=ASCII);
   XMLParser(string file, Grid& g, int format=ASCII);
   XMLParser(string file, Vect<real_t>& v, int format=ASCII);
   XMLParser(const XMLParser &p);
   ~XMLParser();
   void setMaterialNumber(int m);
   int scan();
   int scan(vector<real_t>& t, int type);
   void setFile(string file);
   void set(Mesh& ms, int format=ASCII);
   void set(Grid& gr);
   int get(EqDataType type, vector<PrescriptionPar>& p);
   int get(Mesh& ms, int format=ASCII);
   int get(Mesh& ms, vector<vector<real_t> >& v, string& name);
   int get(Mesh& ms, Vect<real_t>& v, real_t time=-1, string name="ANYTHING", int format=ASCII);
   int get(Grid& gr, vector<vector<real_t> >& v, string& name);
   int get(Grid& gr, Vect<real_t>& v, real_t time=-1, string name="ANYTHING", int format=ASCII);
   int get(Vect<real_t>& v, real_t time=-1, string name="ANYTHING", int format=ASCII);
   int get(Vect<real_t>& v, const string& name);
   int get(IPF& ipf);
   int get(Domain& dm);
   int get(Tabulation& t);
   int get(DMatrix<real_t>& A);
   int get(Matrix<real_t>* A);
   int getArray(Vect<real_t>& A);
   int getMaterial();
   int get(int type, Vect<real_t>& v, real_t time=0);
   size_t getNbDOF() const { return _nb_dof; }

 protected:
   std::ifstream _is;
   bool _is_opened, _is_closed, _set_mesh, _set_grid, _set_vector, _set_file;
   bool _set_prescription, _set_matrix, _set_domain, _prescription_opened, _compact, _value;
   real_t _x, _y, _z, _time, _sought_time, _val;
   int _access, _type, _cm, _format, _var, _code, _rtype;
   EqDataType _prescription_type;
   string _file, _mesh_file, _vect_file, _el_shape, _sd_shape, _name, _sought_name, _tag_name, _xml, _mat;
   size_t _dof, _label, _nb_dof, _dim, _nb_nodes, _nb_elements, _nb_sides, _nb_edges, _tab_size, _vect_size;
   bool _scan;
   size_t _nb_el_nodes, _nb_sd_nodes, _nb_mat, _all_steps, _nb_funct , _nb_rows, _nb_cols;
   size_t _ik1, _ik2, _dk1, _dk2, _ck, _mk, _pk, _dk, _nb_var, _nx, _ny, _nz, _nt;
   DOFSupport _dof_support;
   mutable Mesh *_theMesh;
   Grid *_theGrid;
   Vect<real_t> *_v, *_theArray;
   XMLParser *_parser;
   IPF *_ipf;
   Domain *_theDomain;
   Tabulation *_theTabulation;
   Matrix<real_t> *_theMatrix;
   PrescriptionPar _par;
   vector<PrescriptionPar> *_vp;
   vector<real_t> *_ft;
   vector<vector<real_t> > *_V;
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
   void read_matrix(const SString& ai);
   void read_matrix_data(const vector<string> &tokens, vector<string>::iterator &it);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif
