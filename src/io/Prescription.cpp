/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                         Implementation of class 'Prescription'

  ==============================================================================*/

#include "io/Prescription.h"
#include "mesh/Mesh.h"
#include "io/XMLParser.h"
#include "io/fparser/fparser.h"
#include "linear_algebra/Vect_impl.h"
#include "OFELIException.h"

extern FunctionParser theParser;

namespace OFELI {

Prescription::Prescription()
{
}


Prescription::Prescription(Mesh&              mesh,
                           const std::string& file)
{
   _theMesh = &mesh;
   _file = file;
}


Prescription::~Prescription() { }


int Prescription::get(EqDataType    type,
                      Vect<real_t>& v,
                      real_t        time,
                      size_t        dof)
{
   _time = time;
   PrescriptionPar par;
   par.dof = dof;
   int ret;
   _v = &v;
   _v->clear();
   XMLParser p(_file,*_theMesh,XMLParser::PRESCRIBE);
   p.get(type,_p);
   for (size_t k=0; k<_p.size(); k++) {
      if ((ret=PARSE(_p[k].fct,"x,y,z,t")) != -1)
         throw OFELIException("In Prescription::get(int,Vect<real_t>,real_t,size_t):"
                              " Error in expression parsing.");
      if (dof) {
         if (type==BOUNDARY_CONDITION)
            get_boundary_condition(k,dof);
         else if (type==BOUNDARY_FORCE || type==TRACTION || type==FLUX)
            get_boundary_force(k,dof);
         else if (type==INITIAL_FIELD || type==BODY_FORCE || type==SOLUTION)
            get_vector(k,dof);
         else if (type==POINT_FORCE)
            get_point_force(k,dof);
         else
            throw OFELIException("In Prescription::get(int,Vect<real_t>,real_t,size_t):"
                                 " Type "+itos(type)+" not found.");
      }
      else {
         if (type==BOUNDARY_CONDITION)
            get_boundary_condition(k);
         else if (type==BOUNDARY_FORCE || type==TRACTION || type==FLUX)
            get_boundary_force(k);
         else if (type==INITIAL_FIELD || type==BODY_FORCE || type==SOLUTION)
            get_vector(k);
         else if (type==POINT_FORCE)
            get_point_force(k);
         else
            throw OFELIException("In Prescription::get(int,Vect<real_t>,real_t,size_t):"
                                 " Type "+itos(type)+" not found.");
      }
   }
   return 0;
}


void Prescription::get_point_force(size_t k)
{
   int err;
   bool p = _p[k].bx || _p[k].by || _p[k].bz;
   MESH_ND {
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         real_t e = 0;
         if (_p[k].bx)
            e += (_p[k].x-The_node.getX())*(_p[k].x-The_node.getX());
         if (_p[k].by)
            e += (_p[k].y-The_node.getY())*(_p[k].y-The_node.getY());
         if (_p[k].bz)
            e += (_p[k].z-The_node.getZ())*(_p[k].z-The_node.getZ());
         if (i==_p[k].dof && sqrt(e)<OFELI_EPSMCH && p) {
            (*_v)(node_label,i) = EVAL_XT(The_node.getCoord(),_time);
            if ((err=EVAL_ERR))
               throw OFELIException("In Prescription::get_point_force(size_t): "
                                    "Error in expression evaluation.");
         }
      }
   }
}


void Prescription::get_point_force(size_t k,
                                   size_t dof)
{
   int err;
   MESH_ND {
      real_t e = 0;
      if (_p[k].bx)
         e += (_p[k].x-The_node.getX())*(_p[k].x-The_node.getX());
      if (_p[k].by)
         e += (_p[k].y-The_node.getY())*(_p[k].y-The_node.getY());
      if (_p[k].bz)
         e += (_p[k].z-The_node.getZ())*(_p[k].z-The_node.getZ());
      if (sqrt(e)<OFELI_EPSMCH && _p[k].dof==dof) {
         (*_v)(node_label,dof) = EVAL_XT(The_node.getCoord(),_time);
         if ((err=EVAL_ERR))
            throw OFELIException("In Prescription::get_point_force(size_t,size_t): "
                                 "Error in expression evaluation.");
      }
   }
}


void Prescription::get_boundary_condition(size_t k,
                                          size_t dof)
{
   int err;
   MESH_ND {
      if (the_node->getCode(dof)==_p[k].code) {
         (*_v)(node_label,dof) = EVAL_XT(The_node.getCoord(),_time);
         if ((err=EVAL_ERR))
            throw OFELIException("In Prescription::get_boundary_condition(size_t,size_t): "
                                 "Error in expression evaluation.");
      }
   }
}


void Prescription::get_boundary_condition(size_t k)
{
   size_t l=0;
   int err;
   MESH_ND {
      for (size_t i=1; i<=The_node.getNbDOF(); i++) {
         l++;
         if (The_node.getCode(_p[k].dof)==_p[k].code && i==_p[k].dof) {
            (*_v)(l) = EVAL_XT(The_node.getCoord(),_time);
            if ((err=EVAL_ERR))
               throw OFELIException("In Prescription::get_boundary_condition(size_t): "
                                    "Error in expression evaluation.");
         }
      }
   }
}


void Prescription::get_boundary_force(size_t k,
                                      size_t dof)
{
   int err;
   MESH_SD {
      if (The_side.getCode(_p[k].dof)==_p[k].code && (_p[k].dof==dof||dof==0)) {
         real_t z = EVAL_XT(The_side.getCenter(),_time);
         for (size_t i=1; i<=The_side.getNbNodes(); ++i)
            (*_v)(The_side(i)->n()) = z;
         if ((err=EVAL_ERR))
            throw OFELIException("In Prescription::get_boundary_force(size_t,size_t): "
                                 "Error in expression evaluation.");
      }
   }
}


void Prescription::get_boundary_force(size_t k)
{
  int err=0;
   MESH_SD {
      if (The_side.getGlobalCode()==9998) {
         (*_v)(The_side.getDOF(_p[k].dof)) = EVAL_XT(The_side.getCenter(),_time);
         if ((err=EVAL_ERR))
            throw OFELIException("In Prescription::get_boundary_force(size_t): "
                                 "Error in expression evaluation.");
      }
      else {
         for (size_t j=1; j<=The_side.getNbDOF(); j++) {
            if (The_side.getCode(_p[k].dof)==_p[k].code && j==_p[k].dof) {
               (*_v)(The_side.getDOF(j)) = EVAL_XT(The_side.getCenter(),_time);
               if ((err=EVAL_ERR))
                  throw OFELIException("In Prescription::get_boundary_force(size_t): "
                                       "Error in expression evaluation.");
            }
         }
      }
   }
}


void Prescription::get_vector(size_t k,
                              size_t dof)
{
   int err;
   MESH_ND {
      if (_p[k].dof==dof || dof==0) {
         (*_v)(node_label) = EVAL_XT(The_node.getCoord(),_time);
         if ((err=EVAL_ERR))
            throw OFELIException("In Prescription::get_vector(size_t,size_t): "
                                 "Error in expression evaluation.");
      }
   }
}


void Prescription::get_vector(size_t k)
{
   int err;
   MESH_ND {
      for (size_t i=1; i<=the_node->getNbDOF(); ++i) {
         if (i==_p[k].dof) {
            (*_v)(node_label,i) = EVAL_XT(The_node.getCoord(),_time);
            if ((err=EVAL_ERR))
               throw OFELIException("In Prescription::get_vector(size_t): "
                                    "Error in expression evaluation.");
         }
      }
   }
}

} /* namespace OFELI */
