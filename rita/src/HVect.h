/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2024 Rachid Touzani

    This file is part of rita.

    rita is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rita is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

  ==============================================================================

                  Definition and implementation of class 'HVect'

  ==============================================================================*/

#pragma once

#include "linear_algebra/Vect_impl.h"
#include "io/IOField.h"
using std::vector;
using namespace OFELI;

namespace RITA {

class HVect
{

 private:

    string name;
    vector<Vect<double> > vs;
    vector<double> ts;

 public:

    HVect() : nt(0), size(0) { }

    HVect(int n) : nt(n), size(0) { }

    void setSize(size_t n) { size = n; }

    void set(Vect<double>& v, double t)
    {
       nt++;
       ts.push_back(t);
       v.setTime(t);
       vs.push_back(v);
    }

    Vect<double> *get(int n)
    {
       if (n<=nt)
          return &vs[n-1];
       else
          return nullptr; 
    }

    Vect<double> *get(double t)
    {
       for (int i=0; i<nt; ++i) {
          if (ts[i]==t)
             return &vs[i];
       }
       return nullptr;
    }

    double getTime(int i) const { return ts[i-1]; } 

    ~HVect() {}

    int nt, size;

   int saveOFELI(const string &file, int e=1)
   {
      IOField ff(file,IOField::OUT);
      for (int n=0; n<nt; n+=e)
         ff.put(vs[n]);
      return 0;
   }

   int saveGnuplot(const string &file, int e=1)
   {
      ofstream ff(file);
      for (int n=0; n<nt; n+=e) {
         Vect<double> &v = vs[n];
         ff << ts[n];
         for (size_t i=0; i<v.size(); ++i)
            ff << "  " << v[i];
         ff << endl;
      }
      return 0;
   }

   int saveGmsh(const string &file, int e=1)
   {
      size_t nb_en=0;
      ofstream ff(file);

      Vect<double> &v = *get(1);
      if (v.WithMesh()==false)
         return 1;
      Mesh &ms = v.getMesh();
      size_t nb_dof = v.getNbDOF();

      ff << "View \"" << v.getName() << "\" {" << endl;
      switch (ms.getDim()) {

         case 1:
            for (auto const& el: ms.theElements) {
               ff << "SL(";
               ff << (*el)(1)->getX() <<  ", 0., 0., " << (*el)(2)->getX() <<  ", 0., 0. ) {" << endl;
               for (int n=0; n<nt; n+=e) {
                  Vect<double> &v = *get(n+1);
                  ff << v((*el)(1)->n(),1) << "," << v((*el)(2)->n(),1);
                  if (n<nt-1 && n+e<nt)
                     ff << ",";
                  ff << endl;
               }
               ff << "};" << endl;
            }
            ff << "};" << endl;
            break;

         case 2:
            for (auto const& el: ms.theElements) {
               if (nb_dof==1)
                  ff << 'S';
               else
                  ff << 'V';
               if ((nb_en=el->getNbNodes())==3)
                  ff << "T(";
               else if (nb_en==4)
                  ff << "Q(";
               for (size_t k=1; k<nb_en; ++k)
                  ff << (*el)(k)->getX() << "," << (*el)(k)->getY() << ",0.,";
               ff << (*el)(nb_en)->getX() << "," << (*el)(nb_en)->getY() << ",0.) {" << endl;
               for (int n=0; n<nt; n+=e) {
                  Vect<double> &v = *get(n+1);
                  for (size_t k=1; k<nb_en; ++k) {
                     ff << v((*el)(k)->n(),1);
                     if (nb_dof > 1)
                        ff << "," << v((*el)(k)->n(),2) << ",0.0";
                     ff << "," << endl;
                  }
                  ff << v((*el)(nb_en)->n(),1);
                  if (nb_dof > 1)
                     ff << "," << v((*el)(nb_en)->n(),2) << ",0.0";
                  if (n<nt-1 && n+e<nt)
                     ff << ",";
                  ff << endl;
               }
               ff << "};" << endl;
            }
            ff << "};" << endl;
            break;

         case 3:
            for (auto const& el: ms.theElements) {
               if (nb_dof==1)
                  ff << 'S';
               else
                  ff << 'V';
               if ((nb_en=el->getNbNodes())==4)
                  ff << "S(";
               else if (nb_en==8)
                  ff << "H(";
               else if (nb_en==6)
                  ff << "I(";
               else if (nb_en==5)
                  ff << "Y(";
               for (size_t k=1; k<nb_en; ++k)
                  ff << (*el)(k)->getX() << "," << (*el)(k)->getY() << "," << (*el)(k)->getZ() << ",";
               ff << (*el)(nb_en)->getX() << "," << (*el)(nb_en)->getY() << "," << (*el)(nb_en)->getZ() << ") {" << endl;
               for (int n=0; n<nt; n+=e) {
                  Vect<double> &v = *get(n+1);
                  for (size_t k=1; k<nb_en; ++k) {
                     ff << v((*el)(k)->n(),1);
                     if (nb_dof>1)
                        ff << "," << v((*el)(k)->n(),2) << "," << v((*el)(k)->n(),3);
                     ff << "," << endl;
                  }
                  ff << v((*el)(nb_en)->n(),1);
                  if (nb_dof>1)
                     ff << "," << v((*el)(nb_en)->n(),2) << "," << v((*el)(nb_en)->n(),3);
                  if (n<nt-1 && n+e<nt)
                     ff << ",";
                  ff << endl;
               }
               ff << "};" << endl;
            }
            ff << "};" << endl;
            break;
      }
      return 0;
   }

   int saveVTK(const string &file, int e=1)
   {
      string of="", proj=file.substr(0,file.rfind("."));
      map<int,int> nen = {{LINE,2},{TRIANGLE,3},{QUADRILATERAL,4},{TETRAHEDRON,4},
                          {HEXAHEDRON,8},{PENTAHEDRON,6}};
      map<int,int> ShCode = {{LINE,3},{TRIANGLE,5},{QUADRILATERAL,9},{TETRAHEDRON,10},
                             {HEXAHEDRON,12},{PENTAHEDRON,13}};

      Vect<double> &v = *get(1);
      if (v.WithMesh()==false)
         return 1;
      size_t nb_dof = v.getNbDOF();
      Mesh &ms = v.getMesh();

      int sz=0, nn=0;
      for (auto const& el: ms.theElements)
         sz += nen[el->getShape()] + 1;
      for (int n=0; n<nt; n+=e) {
         Vect<double> &v = *get(n+1);
         of = proj + ".vtk";
         if (nt>1)
            of = proj + "-" + zeros(nn++) + ".vtk";
         cout << "   Storing time step " << n << " in file " << of << endl;
         ofstream pf(of.c_str());
         pf << setprecision(16) << std::scientific;
         pf << "# vtk DataFile Version 2.0\n# Imported from OFELI files\nASCII" << endl;
         pf << "DATASET UNSTRUCTURED_GRID\nPOINTS " << ms.getNbNodes() << " double" << endl;
         for (auto const& nd: ms.theNodes)
            pf << nd->getX() << "  " << nd->getY() << "  " << nd->getZ() << endl;
         pf << "\nCELLS " << ms.getNbElements() << setw(10) << sz << endl;
         for (auto const& el: ms.theElements) {
            pf << setw(8) << nen[el->getShape()];
            for (int i=1; i<=nen[el->getShape()]; i++)
               pf << setw(10) << (*el)(i)->n()-1;
            pf << endl;
         }
         pf << "\nCELL_TYPES  " << ms.getNbElements() << endl;
         int k=0;
         for (auto const& el: ms.theElements) {
            pf << setw(4) << ShCode[el->getShape()];
            if (++k%30 == 0)
               pf << endl;
         }
         pf << "\nPOINT_DATA  " << ms.getNbNodes() << endl;
         if (nb_dof==1)
            pf << "SCALARS  " << v.getName() << "  double  1\nLOOKUP_TABLE  default" << endl;
         else
            pf << "VECTORS  " << v.getName() << "  double" << endl;

         for (auto const& nd: ms.theNodes) {
            pf << v(nd->n(),1) << " ";
            if (nb_dof>1) {
               pf << v(nd->n(),2) << " ";
               if (nb_dof>2)
                  pf << v(nd->n(),3) << " ";
               else
                  pf << 0. << " ";
            }
            pf << endl;
         }
      }
      return 0;
   }

   int saveTecplot(const string& file, int e=1)
   {
      map<int,string> shape = {{LINE,"LINESEG"},{QUADRILATERAL,"QUADRILATERAL"},{TRIANGLE,"TRIANGLE"},
                               {TETRAHEDRON,"TETRAHEDRON"},{HEXAHEDRON,"HEXAHEDRON"},{PENTAHEDRON,"HEXAHEDRON"}};
      Vect<double> &v = *get(1);
      if (v.WithMesh()==false)
         return 1;
      size_t nb_dof = v.getNbDOF();
      Mesh &ms = v.getMesh();

      ofstream ff(file.c_str());
      ff.setf(ios::right|ios::scientific);
      int count = 0;
      for (int n=0; n<nt; n+=e) {
         Vect<double> &v = *get(n+1);
         if (count==0) {
            ff << "TITLE = \" \"\n" << endl;
            ff << "VARIABLES = \"X\", \"Y\"";
            if (ms.getDim()==3)
               ff << ", \"Z\"";
            if (nb_dof==1)
               ff << ", \"T\"";
            else if (nb_dof==2)
               ff << ", \"UX\", \"UY\"";
            else if (nb_dof==3)
               ff << ", \"UX\", \"UY\", \"UZ\"";
            ff << endl;
         }
         ff << "\nZONE T=\"" << "step-" << n << "\", N=" << ms.getNbNodes() << ", E="
            << ms.getNbElements() << ", F=FEPOINT, ET=" << shape[ms.getShape()]
            << ", SOLUTIONTIME=" << v.getTime();
         if (count) {
            ff << ", D=(1,";
            if (ms.getDim()>1)
               ff << "2,";
            if (ms.getDim()==3)
               ff << "3,";
            ff << "FECONNECT)";
         }
         ff << endl;
         for (auto const& nd: ms.theNodes) {
            if (count==0)
               for (size_t i=1; i<=ms.getDim(); i++)
                  ff << "  " << nd->getCoord(i);
            for (size_t j=1; j<=nb_dof; j++)
               ff << "  " << v(nd->n(),j);
            ff << endl;
         }
         if (count==0) {
            for (auto const& el: ms.theElements) {
               for (size_t i=1; i<=el->getNbNodes(); ++i)
                  ff << setw(10) << (*el)(i)->n();
               ff << endl;
            }
         }
         count++;
      }
      ff.close();
      return 0;
   }
};


inline ostream &operator<<(ostream&     s,
                           const HVect& v)
{
   s << "Contents of history vector:" << endl;
   s << "Number of time steps: " << v.nt << endl;
   return s;
}

} /* namespace RITA */