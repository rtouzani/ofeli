
/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2025 Rachid Touzani

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

                        Definition of class 'mesh'

  ==============================================================================*/

#pragma once

#include <fstream>
#include <map>
#include "mesh/Mesh.h"
#include "mesh/Domain.h"
#include "configure.h"

using std::map;

namespace RITA {

class rita;
class data;

class mesh
{

 public:

    mesh(rita *r, cmd* command, configure* config);
    ~mesh();
    int run();
    void setVerbose(int verb) { _verb = verb; }
    void set(cmd* command) { _cmd = command; }
    void set(OFELI::Mesh* ms) { _theMesh = ms; }
    OFELI::Mesh* get() const { return _theMesh; }
    OFELI::Mesh* getMesh() const { return _theMesh; }
    void set(configure* config) { _configure = config; }
    bool MeshIsOK() { return _generated; }
    string mesh_name;

 private:

   struct Point { int n; double x, y, z, h; };
   struct Entity { int nb; int type; vector<int> l; };
   struct Subdomain { int ln, orientation, code; };

   rita *_rita;
   data *_data;
   string _pr, _PR=">mesh>";
   bool _saved, _generated, _geo;
   int _verb, _dim, _nb_dof, _key, _generator;
   OFELI::Mesh *_theMesh;
   OFELI::Domain *_theDomain;
   string _mesh_file;
   typedef void (mesh::* MeshData_Ptr)();
   static MeshData_Ptr MESH_DATA[20];
   std::vector<string> _kw;
   int _nb_Ccontour, _nb_Scontour, _nb_Vcontour;
   int _nb_sub_domain, _nb_point, _nb_curve, _nb_surface, _nb_volume;
   Point _point;
   Entity _curv;
   Subdomain _sd;
   map<int,Point> _points;
   map<int,Entity> _curve, _surface, _volume;
   map<int,Entity> _Ccontour, _Scontour, _Vcontour;
   map<int,Entity> _Pcode, _Ccode, _Scode, _Vcode;
   vector<Subdomain> _subdomains;
   configure *_configure;
   cmd *_cmd;

   void List();
   int set1D();
   int setRectangle();
   int setCube();
   int setPoint();
   int setCurve();
   int setSurface();
   int setVolume();
   int setLineContour();
   int setContour();
   int setSubDomain();
   int setCode();
   void saveDomain(const string& file);
   int Generate();
   int setNbDOF();
   int Plot();
   int Read();
   int Save();
   void saveGeo(const string& file);
   void setConfigure();
};

} /* namespace RITA */