/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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
#include "linear_algebra/Vect_impl.h"
#include "OFELIException.h"

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
   _v = &v;
   _v->clear();
   XMLParser p(_file,*_theMesh,XMLParser::PRESCRIBE);
   p.get(type,_p);
   _theFct.resize(_p.size());
   for (size_t k=0; k<_p.size(); k++) {
      _theFct[k].set(_p[k].fct);
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
                                 " Type "+std::to_string(type)+" not found.");
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
                                 " Type "+std::to_string(type)+" not found.");
      }
   }
   return 0;
}


Vect<real_t> &Prescription::get(EqDataType type,
                                real_t     time,
                                size_t     dof)
{
   _time = time;
   PrescriptionPar par;
   par.dof = dof;
   _v = new Vect<real_t>(*_theMesh);
   _v->clear();
   XMLParser p(_file,*_theMesh,XMLParser::PRESCRIBE);
   p.get(type,_p);
   _theFct.resize(_p.size());
   for (size_t k=0; k<_p.size(); k++) {
      _theFct[k].set(_p[k].fct);
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
            throw OFELIException("In Prescription::get(int,real_t,size_t): Type "+std::to_string(type)+" not found.");
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
            throw OFELIException("In Prescription::get(int,real_t,size_t): Type "+std::to_string(type)+" not found.");
      }
   }
   return *_v;
}


void Prescription::get_point_force(size_t k)
{
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
         if (i==_p[k].dof && sqrt(e)<OFELI_EPSMCH && p)
            (*_v)(node_label,i) = _theFct[k](The_node.getCoord());
      }
   }
}


void Prescription::get_point_force(size_t k,
                                   size_t dof)
{
   MESH_ND {
      real_t e = 0;
      if (_p[k].bx)
         e += (_p[k].x-The_node.getX())*(_p[k].x-The_node.getX());
      if (_p[k].by)
         e += (_p[k].y-The_node.getY())*(_p[k].y-The_node.getY());
      if (_p[k].bz)
         e += (_p[k].z-The_node.getZ())*(_p[k].z-The_node.getZ());
      if (sqrt(e)<OFELI_EPSMCH && _p[k].dof==dof)
         (*_v)(node_label,dof) = _theFct[k](The_node.getCoord());
   }
}


void Prescription::get_boundary_condition(size_t k,
                                          size_t dof)
{
   MESH_ND {
      if (the_node->getCode(dof)==_p[k].code)
         (*_v)(node_label,dof) = _theFct[k](The_node.getCoord());
   }
}


void Prescription::get_boundary_condition(size_t k)
{
   size_t l=1;
   MESH_ND {
      for (size_t i=1; i<=The_node.getNbDOF(); i++, l++) {
         if (The_node.getCode(_p[k].dof)==_p[k].code && i==_p[k].dof)
            (*_v)(l) = _theFct[k](The_node.getCoord());
      }
   }
}


void Prescription::get_boundary_force(size_t k,
                                      size_t dof)
{
   MESH_SD {
      if (The_side.getCode(_p[k].dof)==_p[k].code && (_p[k].dof==dof||dof==0)) {
         for (size_t i=1; i<=The_side.getNbNodes(); ++i)
            (*_v)(The_side(i)->n()) = _theFct[k](The_side.getCenter());
      }
   }
}


void Prescription::get_boundary_force(size_t k)
{
   MESH_SD {
      if (The_side.getGlobalCode()==9998)
         (*_v)(The_side.getDOF(_p[k].dof)) = _theFct[k](The_side.getCenter());
      else {
         for (size_t j=1; j<=The_side.getNbDOF(); j++) {
            if (The_side.getCode(_p[k].dof)==_p[k].code && j==_p[k].dof)
               (*_v)(The_side.getDOF(j)) = _theFct[k](The_side.getCenter());
         }
      }
   }
}


void Prescription::get_vector(size_t k,
                              size_t dof)
{
   MESH_ND {
      if (_p[k].dof==dof || dof==0)
         (*_v)(node_label) = _theFct[k](The_node.getCoord());
   }
}


void Prescription::get_vector(size_t k)
{
   MESH_ND {
      for (size_t i=1; i<=the_node->getNbDOF(); ++i) {
         if (i==_p[k].dof)
            (*_v)(node_label,i) = _theFct[k](The_node.getCoord());
      }
   }
}

} /* namespace OFELI */
