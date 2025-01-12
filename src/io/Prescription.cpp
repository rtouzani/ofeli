/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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
#include "io/XMLParser.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect_impl.h"
#include "OFELIException.h"

namespace OFELI {

Prescription::Prescription()
{
}


Prescription::Prescription(Mesh&         mesh,
                           const string& file)
{
   _theMesh = &mesh;
   _file = file;
}


Prescription::~Prescription() { }


int Prescription::get(EType         type,
                      Vect<real_t>& v,
                      real_t        time,
                      size_t        dof)
{
   _time = time;
   PPar par;
   par.dof = dof;
   _v = &v;
   _v->clear();
   XMLParser p(_file,*_theMesh,EType::PRESCRIBE);
   p.get(type,_p);
   _theFct.resize(_p.size());
   for (size_t k=0; k<_p.size(); k++) {
      _theFct[k].set(_p[k].fct);
      if (type==EType::BOUNDARY_CONDITION)
         get_boundary_condition(k,dof);
      else if (type==EType::BOUNDARY_FORCE)
         get_boundary_force(k,dof);
      else if (type==EType::INITIAL || type==EType::BODY_FORCE || type==EType::SOLUTION)
         get_vector(k,dof);
      else if (type==EType::POINT_FORCE)
         get_point_force(k,dof);
      else
         throw OFELIException("In Prescription::get(EType,Vect<real_t>,real_t,size_t):"
                              " Type not found.");
   }
   return 0;
}


Vect<real_t> &Prescription::get(EType  type,
                                real_t time,
                                size_t dof)
{
   _time = time;
   PPar par;
   par.dof = dof;
   _v = new Vect<real_t>(*_theMesh);
   _v->clear();
   XMLParser p(_file,*_theMesh,EType::PRESCRIBE);
   p.get(type,_p);
   _theFct.resize(_p.size());
   for (size_t k=0; k<_p.size(); k++) {
      _theFct[k].set(_p[k].fct);
      if (type==EType::BOUNDARY_CONDITION)
         get_boundary_condition(k,dof);
      else if (type==EType::BOUNDARY_FORCE)
         get_boundary_force(k,dof);
      else if (type==EType::INITIAL || type==EType::BODY_FORCE || type==EType::SOLUTION)
         get_vector(k,dof);
      else if (type==EType::POINT_FORCE)
         get_point_force(k,dof);
      else
         throw OFELIException("In Prescription::get(EType,real_t,size_t): Type not found.");
   }
   return *_v;
}


void Prescription::getBoundaryCondition(Vect<real_t>& v,
                                        real_t        time,
                                        size_t        dof)
{
   get(EType::BOUNDARY_CONDITION,v,time,dof);
}
 

void Prescription::getDirichlet(Vect<real_t>& v,
                                real_t        time,
                                size_t        dof)
{
   get(EType::DIRICHLET,v,time,dof);
}

 
void Prescription::getInitial(Vect<real_t>& v,
                              real_t        time,
                              size_t        dof)
{
   get(EType::INITIAL,v,time,dof);
}

 
void Prescription::getSolution(Vect<real_t>& v,
                               real_t        time,
                               size_t        dof)
{
   get(EType::SOLUTION,v,time,dof);
}

 
void Prescription::getBodyForce(Vect<real_t>& v,
                                real_t        time,
                                size_t        dof)
{
   get(EType::BODY_FORCE,v,time,dof);
}


void Prescription::getBoundaryForce(Vect<real_t>& v,
                                    real_t        time,
                                    size_t        dof)
{
   get(EType::BOUNDARY_FORCE,v,time,dof);
}


void Prescription::getTraction(Vect<real_t>& v,
                               real_t        time,
                               size_t        dof)             
{
   get(EType::BOUNDARY_FORCE,v,time,dof);
}


void Prescription::getNeumann(Vect<real_t>& v,
                              real_t        time,
                              size_t        dof)
{
   get(EType::BOUNDARY_FORCE,v,time,dof);
}


void Prescription::getFlux(Vect<real_t>& v,
                           real_t        time,
                           size_t        dof)
{
   get(EType::FLUX,v,time,dof);
}


void Prescription::getPointForce(Vect<real_t>& v,
                                 real_t        time,
                                 size_t        dof)
{
   get(EType::POINT_FORCE,v,time,dof);
}


Vect<real_t> &Prescription::getBoundaryCondition(real_t time,
                                                 size_t dof)
{
   return get(EType::BOUNDARY_CONDITION,time,dof);
}


Vect<real_t> &Prescription::getDirichlet(real_t time,
                                         size_t dof)
{
   return get(EType::BOUNDARY_CONDITION,time,dof);
}


Vect<real_t> &Prescription::getInitial(real_t time,
                                       size_t dof)
{
   return get(EType::INITIAL,time,dof);
}


Vect<real_t> &Prescription::getBodyForce(real_t time,
                                         size_t dof)
{
   return get(EType::BODY_FORCE,time,dof);
}


Vect<real_t> &Prescription::getBoundaryForce(real_t time,
                                             size_t dof)
{
   return get(EType::BOUNDARY_FORCE,time,dof);
}


Vect<real_t> &Prescription::getTraction(real_t time,
                                        size_t dof)
{
   return get(EType::BOUNDARY_FORCE,time,dof);
}


Vect<real_t> &Prescription::getFlux(real_t time,
                                    size_t dof)
{
   return get(EType::BOUNDARY_FORCE,time,dof);
}


Vect<real_t> &Prescription::getNeumann(real_t time,
                                       size_t dof)
{
   return get(EType::BOUNDARY_FORCE,time,dof);
}


Vect<real_t> &Prescription::getPointForce(real_t time,
                                          size_t dof)
{
   return get(EType::POINT_FORCE,time,dof);
}


void Prescription::get_point_force(size_t k,
                                   size_t dof)
{
   Fct &f = _theFct[k];
   if (dof) {
      MESH_ND {
         real_t e = 0;
         if (_p[k].bx)
            e += (_p[k].x-The_node.getX())*(_p[k].x-The_node.getX());
         if (_p[k].by)
            e += (_p[k].y-The_node.getY())*(_p[k].y-The_node.getY());
         if (_p[k].bz)
            e += (_p[k].z-The_node.getZ())*(_p[k].z-The_node.getZ());
         if (sqrt(e)<OFELI_EPSMCH && _p[k].dof==dof)
            (*_v)(node_label,dof) = f(The_node.getCoord());
      }
   }
   else {
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
               (*_v)(node_label,i) = f(The_node.getCoord());
         }
      }
   }
}


void Prescription::get_boundary_condition(size_t k,
                                          size_t dof)
{
   Fct &f = _theFct[k];
   int c = _p[k].code;
   if (dof) {
      MESH_ND {
         if (The_node.getCode(dof)==c)
            (*_v)(node_label,dof) = f(The_node.getCoord());
      }
   }
   else {
      size_t l=1;
      MESH_ND {
         for (size_t i=1; i<=The_node.getNbDOF(); i++, l++) {
            if (The_node.getCode(_p[k].dof)==c && i==_p[k].dof)
               (*_v)(l) = f(The_node.getCoord());
         }
      }
   }
}


void Prescription::get_boundary_force(size_t k,
                                      size_t dof)
{
   Fct &f = _theFct[k];
   int c = _p[k].code;
   if (dof) {
      MESH_SD {
         if (The_side.getCode(_p[k].dof)==c && _p[k].dof==dof) {
            for (size_t i=1; i<=The_side.getNbNodes(); ++i)
               (*_v)(The_side(i)->n()) = f(The_side.getCenter());
         }
      }
   }
   else {
      MESH_SD {
         if (The_side.getGlobalCode()==BCType::CONTACT_BC)
            (*_v)(The_side.getDOF(_p[k].dof)) = f(The_side.getCenter());
         else {
            for (size_t j=1; j<=The_side.getNbDOF(); j++) {
               if (The_side.getCode(_p[k].dof)==c && j==_p[k].dof)
                  (*_v)(The_side.getDOF(j)) = f(The_side.getCenter());
            }
         }
      }
   }
}


void Prescription::get_vector(size_t k,
                              size_t dof)
{
   Fct &f = _theFct[k];
   if (dof) {
      MESH_ND {
         if (_p[k].dof==dof)
            (*_v)(node_label) = f(The_node.getCoord());
      }
   }
   else {
      MESH_ND {
         for (size_t i=1; i<=The_node.getNbDOF(); ++i) {
            if (i==_p[k].dof)
               (*_v)(node_label,i) = f(The_node.getCoord());
         }
      }
   }
}

} /* namespace OFELI */
