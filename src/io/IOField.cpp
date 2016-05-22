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

                         Implementation of class 'IOField'

  ==============================================================================*/

#include "io/IOField.h"
#include "io/saveField.h"

namespace OFELI {

IOField::IOField()
{
   _is_closed = false;
   _field_opened = false;
   _no_mesh_file = true;
   _compact = true;
}

  /*
IOField::IOField(const string &file, int access, bool compact)
{
   _is_closed = false;
   _compact = compact;
   _field_opened = false;
   _no_mesh_file = true;
   open(file,access);
}*/



IOField::IOField(const string& file,
                 AccessType    access,
                 bool          compact)
        : _field_opened(false), _compact(compact), _no_mesh_file(true), _theMesh(NULL)
{
   _is_opened = false;
   _ipf = NULL;
   _set_mesh = false;
   _set_field = true;
   _set_file = true;
   _scan = false;
   _nb_dof = 1;
   _dof_support = 0;
   _nb_nodes = _nb_elements = _nb_sides = _nb_edges = 0;
   _v = NULL;
   _file = file;
   _parser = NULL;
   _verb = 0;
   _is_closed = false;
   open(file,access);
   XMLParser::open();
}


IOField::IOField(const string& file,
                 AccessType    access,
                 const string& name)
        : _field_name(name), _field_opened(false), _compact(true),
          _no_mesh_file(true), _theMesh(NULL)
{
   _is_opened = false;
   _ipf = NULL;
   _set_mesh = false;
   _set_field = true;
   _set_file = true;
   _scan = false;
   _nb_dof = 1;
   _dof_support = 0;
   _nb_nodes = _nb_elements = _nb_sides = _nb_edges = 0;
   _v = NULL;
   _file = file;
   _parser = NULL;
   _verb = 0;
   _is_closed = false;
   open(file,access);
   XMLParser::open();
}


IOField::IOField(const string& mesh_file,
                 const string& file,
                 Mesh&         ms,
                 AccessType    access,
                 bool          compact)
        : _field_opened(false), _compact(compact), _no_mesh_file(false), _theMesh(&ms)
{
   _is_opened = false;
   _ipf = NULL;
   _set_mesh = true;
   _set_field = true;
   _set_file = true;
   _scan = false;
   _nb_dof = 1;
   _dof_support = 1;
   _nb_nodes = _theMesh->getNbNodes();
   _nb_elements = _theMesh->getNbElements();
   _nb_sides = _theMesh->getNbSides();
   _nb_edges = _theMesh->getNbEdges();
   _v = NULL;
   _file = file;
   _parser = NULL;
   _verb = 0;
   _is_closed = false;
   open(file,access);
   setMeshFile(mesh_file);
   set(ms);
   XMLParser::open();
}


IOField::IOField(const string& file,
                 Mesh&         ms,
                 AccessType    access,
                bool           compact)
        : _field_opened(false), _compact(compact), _no_mesh_file(true), _theMesh(&ms)
{
   _is_opened = false;
   _ipf = NULL;
   _set_mesh = true;
   _set_field = true;
   _set_file = true;
   _scan = false;
   _nb_dof = 1;
   _dof_support = 1;
   _nb_nodes = _theMesh->getNbNodes();
   _nb_elements = _theMesh->getNbElements();
   _nb_sides = _theMesh->getNbSides();
   _nb_edges = _theMesh->getNbEdges();
   _v = NULL;
   _file = file;
   _parser = NULL;
   _verb = 0;
   _is_closed = false;
   open(file,access);
   set(ms);
   XMLParser::open();
}


IOField::~IOField()
{
   if (_access==OUT && _is_closed==false)
      close();
}


void IOField::saveGMSH(string output_file,
                       string mesh_file)
{
   try {
      if (_access!=OUT)
         THROW_RT("saveGMSH(...): File must have been opened with the option OUT");
   }
   CATCH("IOField");
   close();
   _is_closed = true;
   saveGmsh(_file,output_file,mesh_file);
}


void IOField::setMeshFile(const string& file)
{
   _mesh_file = file;
   if (_access==OUT)
      *_of << "<Mesh file=\"" << _mesh_file << "\" />" << endl;
}


void IOField::open()
{
   if (_access == OUT)
      _of = new ofstream(_file.c_str());
   else if (_access == BIN_OUT)
      _of = new ofstream(_file.c_str(),ios::binary);
   else
      ;
   _state = 0;
   _of->setf(ios::right);
   *_of << setprecision(16) << scientific;
   XMLParser::open();
}


void IOField::open(const string& file,
                   AccessType    access)
{
   _file = file;
   _access = access;
   if (access==IN || access==BIN_IN)
      _state = 0;
   else if (_access == OUT) {
      _of = new ofstream(_file.c_str(),ios::out);
      _of->setf(ios::right|ios::scientific);
      *_of << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>" << endl;
      *_of << "<OFELI_File>" << endl;
      *_of << "<info>\n   <title></title>" << endl;
      *_of << "   <date>" << __DATE__ << "</date>" << endl;
      *_of << "   <author></author>\n</info>" << endl;
      _field_name = "abcdefghijklmnopqrstuvwxyz";
   }
   else if (_access==BIN_OUT) {
      _of = new ofstream(_file.c_str(),ios::binary|ios::out);
      _field_opened = false;
      _of->setf(ios::right);
      *_of << setprecision(16) << scientific;
   }
   else
      _is_closed = false;
   _state = 0;
}


void IOField::close()
{
   if (!_is_closed && _access==OUT) {
      _is_closed = true;
      if (_field_opened==true)
         *_of << "</Field>" << endl;
      *_of << "</OFELI_File>" << endl;
   }
   _of->close();
}


void IOField::put(Mesh& ms)
{
   _theMesh = &ms;
   string sh[10] = {"line","triangle","quadrilateral","tetrahedron","hexahedron","pentahedron"};
   size_t n = _theMesh->getNbDOF()/_theMesh->getNbNodes();
   *_of << "<Mesh dim=\"" << ms.getDim() << "\" nb_dof=\"" << n << "\">" << endl;
   if (ms.getNbNodes()>0) {
      size_t k = 0;
      *_of << "   <Nodes>" << endl;
      mesh_nodes(*_theMesh) {
         for (size_t i=1; i<=_theMesh->getDim(); i++)
            *_of << "  " << The_node.getCoord(i);
         size_t m = 0;
         for (size_t j=1; j<=n; j++)
            m += The_node.getCode(j)*size_t(pow(10.,real_t(n-j)));
         *_of << setw(6) << m << "   ";
         k++;
         if (k%5==0)
            *_of << endl;
      }
      if (k%5!=0)
         *_of << endl;
      *_of << "   </Nodes>" << endl;
   }
   if (_theMesh->getNbElements()>0) {
      size_t k = 0;
      size_t nbn = _theMesh->getPtrElement(1)->getNbNodes();
      string shape = sh[_theMesh->getPtrElement(1)->getShape()];
      *_of << "   <Elements shape=\"" << shape << "\"  nodes=\"" << nbn << "\">" << endl;
      mesh_elements(*_theMesh) {
         for (size_t i=1; i<=nbn; i++)
            *_of << setw(10) << The_element.getNodeLabel(i);
         *_of << setw(10) << The_element.getCode() << "   ";
         k++;
         if (k%5==0)
            *_of << endl;
      }
      if (k%5!=0)
         *_of << endl;
      *_of << "   </Elements>" << endl;
   }
   if (_theMesh->getNbSides()>0) {
      size_t k = 0;
      size_t nbn = _theMesh->getPtrSide(1)->getNbNodes();
      string shape = sh[_theMesh->getPtrSide(1)->getShape()];
      *_of << "   <Sides shape=\"" << shape << "\"  nodes=\"" << nbn << "\">" << endl;
      size_t m = 0;
      mesh_sides(*_theMesh) {
         for (size_t i=1; i<=nbn; i++)
            *_of << setw(8) << the_side->getNodeLabel(i);
         m = 0;
         for (size_t j=1; j<=n; j++)
            m += The_side.getCode(j)*size_t(pow(10.,real_t(n-j)));
         *_of << setw(6) << m << "   ";
         k++;
         if (k%5==0)
            *_of << endl;
      }
      if (k%5!=0)
         *_of << endl;
      *_of << "   </Sides>" << endl;
   }
   /*   if (_nb_mat>1 || ms.getMaterial().getName(1)!="Generic") {
      *_of << "   <Material>" << endl;
      for (size_t i=0; i<_nb_mat; i++)
         *_of << setw(9) << ms.getMaterial().getCode(i+1) << "   "
              << ms.getMaterial().getName(i+1) << endl;
      *_of << "   </Material>" << endl;
      }*/
   *_of << "</Mesh>" << endl;
}


void IOField::put(const Vect<real_t>& v)
{
   try {
      if (_access==OUT || _access==BIN_OUT) {
         if (v.getName() != _field_name) {
            if (_field_opened==true)
               *_of << "</Field>" << endl;
            _field_name = v.getName();
            if (v.getDOFType()==NONE) {
               *_of << "<Field name=\"" << v.getName() << "\" type=\"None\"";
               *_of << " nb_dof=\"1\">" << endl;
            }
            else {
               if (v.getDOFType()==NODE_FIELD)
                  *_of << "<Field name=\"" << v.getName() << "\" type=\"Node\"";
               else if (v.getDOFType()==ELEMENT_FIELD)
                  *_of << "<Field name=\"" << v.getName() << "\" type=\"Element\"";
               else if (v.getDOFType()==SIDE_FIELD)
                  *_of << "<Field name=\"" << v.getName() << "\" type=\"Side\"";
               *_of << " nb_dof=\"" << v.getNbDOF() << "\">" << endl;
	    }
            _field_opened = true;
         }
         if (_compact)
            *_of << v.getTime() << endl;
         else
            *_of << "   <Step time=\"" << v.getTime() << "\">" << endl;
         if (_access==BIN_OUT) {
            string s; 
	 //         s = *(v.getBinary());
            *_of << s << endl;
         }
         else {
            size_t l = 0;
            for (size_t i=1; i<=v.getNb(); i++) {
               for (size_t j=1; j<=v.getNbDOF(); j++)
                  *_of << setprecision(8) << setw(18) << v(i,j);
               l++;
               if (l%(10/v.getNbDOF())==0)
                  *_of << endl;
            }
            if (l%(10/v.getNbDOF())!=0)
               *_of << endl;
         }
         if (!_compact)
            *_of << "   </Step>" << endl;
      }
      else if (_access==BIN_OUT)
         _of->write(reinterpret_cast<const char *>(&v),sizeof(v));
      else
         THROW_RT("put(Vect,real_t): This instance of IOField was constructed for input");
   }
   CATCH("IOField");
   _state++;
   _field_opened = true;
}


#ifdef USE_PETSC
void IOField::put(const PETScVect<real_t>& v)
{
   try {
      if (_access==OUT || _access==BIN_OUT) {
         if (v.getName() != _field_name) {
            if (_field_opened==true)
               *_of << "</Field>" << endl;
            _field_name = v.getName();
            if (v.getDOFType()==NONE) {
               *_of << "<Field name=\"" << v.getName() << "\" type=\"None\"";
               *_of << " nb_dof=\"1\">" << endl;
            }
            else {
               if (v.getDOFType()==NODE_FIELD)
                  *_of << "<Field name=\"" << v.getName() << "\" type=\"Node\"";
               else if (v.getDOFType()==ELEMENT_FIELD)
                  *_of << "<Field name=\"" << v.getName() << "\" type=\"Element\"";
               else if (v.getDOFType()==SIDE_FIELD)
                  *_of << "<Field name=\"" << v.getName() << "\" type=\"Side\"";
               *_of << " nb_dof=\"" << v.getNbDOF() << "\">" << endl;
            }
            _field_opened = true;
         }
         if (_compact)
            *_of << v.getTime() << endl;
         else
            *_of << "   <Step time=\"" << v.getTime() << "\">" << endl;
         if (_access==BIN_OUT) {
            string s; 
	 //         s = *(v.getBinary());
            *_of << s << endl;
         }
         else {
            size_t l = 0;
            for (size_t i=1; i<=v.getNb(); i++) {
               for (size_t j=1; j<=v.getNbDOF(); j++)
                  *_of << setprecision(8) << setw(18) << v(i,j);
               l++;
               if (l%(10/v.getNbDOF())==0)
                  *_of << endl;
            }
            if (l%(10/v.getNbDOF())!=0)
               *_of << endl;
         }
         if (!_compact)
            *_of << "   </Step>" << endl;
      }
      else if (_access==BIN_OUT)
         _of->write(reinterpret_cast<const char *>(&v),sizeof(v));
      else
         THROW_RT("put(Vect,real_t): This instance of IOField was constructed for input");
   }
   CATCH("IOField");
   _state++;
   _field_opened = true;
}
#endif


real_t IOField::get(Vect<real_t>& v)
{
   try {
      if (_access != IN)
         THROW_RT("get(Vect): This instance of IOField was constructed for output");
   }
   CATCH("IOField");
   if (_no_mesh_file) {
      return XMLParser::get(v);
   }
   else
      return XMLParser::get(*_theMesh,v);
}


int IOField::get(Vect<real_t>& v,
                 real_t        t)
{
  //   if (_access != IN)
  //     _e.Message(__FILE__,__LINE__,61);
   return XMLParser::get(*_theMesh,v,t);
}


int IOField::get(Vect<real_t>& v,
                 const string& name)
{
   return XMLParser::get(v,name);
}


int IOField::get(DMatrix<real_t>& A,
                 const string&    name)
{
   Vect<real_t> v(A.getNbRows(),A.getNbColumns());
   int ret = XMLParser::get(v,name);
   for (size_t i=1; i<=A.getNbRows(); i++)
      for (size_t j=1; j<=A.getNbColumns(); j++)
         A(i,j) = v(i,j);
   return ret;
}


int IOField::get(DSMatrix<real_t>& A,
                 const string&     name)
{
   Vect<real_t> v(A.getNbRows(),A.getNbRows());
   int ret = XMLParser::get(v,name);
   for (size_t i=1; i<=A.getNbRows(); i++)
      for (size_t j=1; j<=i; j++)
         A(i,j) = v(i,j);
   return ret;
}

} /* namespace OFELI */
