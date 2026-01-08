/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2026 Rachid Touzani

    This file is part of rita.

    rita is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or

    rita is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

  ==============================================================================
 
                      Implementation of class 'mesh'

  ==============================================================================*/

#include <iostream>
#include <stdlib.h>
#include "mesh.h"
#include "mesh/saveMesh.h"
#include "rita.h"
#include "data.h"
#include "calc.h"
#include "defs.h"
#include "helps.h"

#ifdef USE_GMSH
#include <gmsh.h>
#endif

namespace RITA {

using std::cout;

mesh::mesh(rita*      r,
           cmd*       command,
           configure* config)
     : _rita(r), _saved(false), _generated(false), _geo(false), _generator(0),
       _theMesh(nullptr), _theDomain(nullptr), _nb_Ccontour(0), _nb_Scontour(0),
       _nb_Vcontour(0), _nb_sub_domain(0), _nb_point(0), _nb_curve(0), _nb_surface(0),
       _nb_volume(0), _configure(config), _cmd(command)
{
   mesh_name = "";
}


mesh::~mesh()
{
   if (_theMesh!=nullptr)
      delete _theMesh, _theMesh = nullptr;
   if (_theDomain!=nullptr)
      delete _theDomain, _theDomain = nullptr;
}


int mesh::run()
{
   int ret=0, key=0;
   string fn="";
   _pr = _PR;
   _nb_dof = 1;
   _data = _rita->_data;
   static const vector<string> kw {"1d","rect$angle","cube","point","curve","surface","volume",
                                   "contour","code","gen$erate","nbdof","list","plot","save","read"};
#ifndef USE_GMSH
   _theDomain = new OFELI::Domain;
#endif
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args>1,_pr,"Illegal number of arguments")

   if (nb_args==1) {
      int nb=0;
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available sub-commands: " << Mesh_help << endl;
            RET(1)

         case 101:
            cout << Mesh_Help << endl;
            RET(1)

         case 104:
         case 105:
            *_rita->ofh << " end" << endl;
            RET(0)

         case 106:
            _rita->setEcho(nb_args);
            RET(0)

         default:
            UNKNOWN_ARG(_pr)
      }
   }

   *_rita->ofh << "mesh" << endl;
   while (1) {
      READLINE
      switch (key) {

         case   0:
            ret = set1D();
            break;

         case   1:
            ret = setRectangle();
            break;

         case   2:
            ret = setCube();
            break;

         case   3:
            ret = setPoint();
            break;

         case   4:
            ret = setCurve();
            break;

         case   5:
            ret = setSurface();
            break;

         case   6:
            ret = setVolume();
            break;

         case   7:
            ret = setContour();
            break;

         case   8:
            ret = setCode();
            break;

         case   9:
            Generate();
            break;

         case  10:
            setNbDOF();
            break;

         case  11:
            List();
            break;

         case  12:
            Plot();
            break;

         case  13:
            Save();
            break;

         case  14:
            Read();
            break;

         case 100:
             cout << "Available commands: " << Mesh_help << endl;
             break;

         case 101:
            cout << Mesh_Help << endl;
            break;

         case 102:
            _rita->getLicense();
            break;

         case 103:
            ret = _configure->run();
            _verb = _configure->getVerbose();
            break;

         case 104:
         case 105:
            *_rita->ofh << "  end" << endl;
            RET(ret)

         case 106:
            _rita->setEcho(nb_args);
            break;

         default:
            CHK_MSGB(_cmd->checkEq()==string::npos,_pr,"Unrecognized command: "+_cmd->getToken())
            ret = _rita->_calc->run();
            break;
      }
   }
   RET(1)
}


void mesh::List()
{
   if (_generated==0) {
      cout << "Geometry data\n" << endl;
      cout << "Space dimension:    " << _theMesh->getDim() << endl;
      cout << "Number of nodes:    " << _theMesh->getNbNodes() << endl;
      cout << "Number of elements: " << _theMesh->getNbElements() << endl;
      cout << "Number of sides:    " << _theMesh->getNbSides() << endl;
   }
   else {
      cout << "Mesh data\n" << endl;
      cout << "Space dimension:    " << _theMesh->getDim() << endl;
      cout << "Number of nodes:    " << _theMesh->getNbNodes() << endl;
      cout << "Number of elements: " << _theMesh->getNbElements() << endl;
      cout << "Number of sides:    " << _theMesh->getNbSides() << endl;
   }
}


void mesh::setConfigure()
{
   if (_verb)
      cout << "Setting configuration parameters ..." << endl;
   _configure->setVerbose(_verb);
   _configure->run();
   _verb = _configure->getVerbose();
}


int mesh::set1D()
{
   string fn="";
   _pr = _PR + "1d>";
   _dim = 1;
   _nb_dof = 1;
   double xmin=0., xmax=1.;
   int nb=0, ret=0, ne=10, cmin=0, cmax=0;
   _mesh_file = "rita-1d.m";
   _saved = false;
   _generated = false;
   static const vector<string> kw {"domain","ne","codes","nbdof"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();

   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Mesh_1d_help << endl;
            RET(0)

         case 101:
            cout << Mesh_1d_Help << endl;
            RET(0)

         case   0:
            if (nb==1) {
               xmin = 0.;
               DEF_PAR_R(0,_pr,xmax)
            }
            else if (nb==2) {
               DEF_PAR_R(0,_pr,xmin)
               DEF_PAR_R(1,_pr,xmax)
            }
            else
               MSGR(_pr,"This argument requires 1 or 2 parameters.")
            break;

         case   1:
            DEF_PAR_R(0,_pr,ne)
            break;

         case   2:
            if (nb==1) {
               DEF_PAR_R(0,_pr,cmin)
               cmax = cmin;
            }
            else if (nb==2) {
               DEF_PAR_R(0,_pr,cmin)
               DEF_PAR_R(1,_pr,cmax)
            }
            else
               MSGR(_pr,"This argument requires 1 or 2 parameters.")
            break;

         case   3:
            DEF_PAR_R(0,_pr,_nb_dof)
            break;

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   if (nb_args>0) {
      CHK_MSGR(xmax<=xmin,_pr,"Error in values of xmin and xmax: "+to_string(xmin)+", "+to_string(xmax))
      CHK_MSGR(ne<2,_pr,"Number of elements must be > 1")
      _theMesh = new OFELI::Mesh(xmin,xmax,ne,cmin,cmax,1,size_t(_nb_dof));
      _theMesh->removeImposedDOF();
      _data->addMesh(_theMesh,mesh_name);
      _saved = false;
      _generator = 1, _generated = true;
      *_rita->ofh << "  1d domain=" << xmin << "," << xmax << " codes=" << cmin
                  << "," << cmax << " ne=" << ne << " nbdof=" << _nb_dof << endl;
      if (_verb)
         cout << "1-D mesh complete. Mesh Name: "+mesh_name << endl;
      _pr = _PR;
   }

   else {
      if (_verb) {
         cout << "Default interval: (0,1)\n";
         cout << "Default number of elements: 10\n";
         cout << "Default codes for end nodes: 0 0\n";
         cout << "Default nb of dof: 1" << endl;
      }
      *_rita->ofh << "  1d" << endl;
      int key=0;
      for (;;) {
         READLINE
         switch (key) {

            case   0:
               if (_verb>1)
                  cout << "Setting interval bounds ..." << endl;
               CHK_MSGB(_cmd->setNbArg(2,"Give xmin and xmax."),_pr+"domain>","Missing values of xmin and xmax.")
               DEF_PAR_R(-1,_pr+"domain>",xmin)
               DEF_PAR_R(-1,_pr+"domain>",xmax)
               CHK_MSGB(xmax<=xmin,_pr+"domain>","Value of xmax: "+to_string(xmax)+" must be > xmin: "+to_string(xmin))
               *_rita->ofh << "    domain " << xmin << " " << xmax << endl;
               break;

            case   1:
               if (_verb>1)
                  cout << "Setting interval subdivision ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(1,"Give number of elements."),_pr+"ne>","Missing number of elements.","",1)
               DEF_PAR_R(-1,_pr+"ne>",ne)
               *_rita->ofh << "  ne " << ne << endl;
               break;

            case   2:
               if (_verb>1)
                  cout << "Setting boundary codes ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(2,"Give cmin and cmax."),_pr+"codes>","Missing codes for end nodes.","",1)
               ret  = _cmd->get(cmin);
               ret += _cmd->get(cmax);
               if (!ret)
                  *_rita->ofh << "  codes " << cmin << " " << cmax << endl;
               break;

            case   3:
               if (_verb>1)
                  cout << "Setting number of degrees of freedom ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(1,"Give number of dof per node."),_pr+"nbdof>","Missing number of dof.","",1)
               ret = _cmd->get(_nb_dof);
               if (!ret)
                  *_rita->ofh << "  nbdof " << _nb_dof << endl;
               break;

            case 100:
               cout << "Available commands: " << Mesh_1d_help << endl;
               break;

            case 101:
               cout << Mesh_1d_Help << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               ret = _configure->run();
               _verb = _configure->getVerbose();
               break;

            case 104:
            case 105:
               if (_verb>1)
                  cout << "Getting back to higher level ..." << endl;
               if (!_saved) {
                  if (_theMesh!=nullptr)
                     delete _theMesh, _theMesh = nullptr;
                  _theMesh = new OFELI::Mesh(xmin,xmax,ne,cmin,cmax,1,size_t(_nb_dof));
                  _theMesh->removeImposedDOF();
                  _data->addMesh(_theMesh,mesh_name);
               }
               *_rita->ofh << "  end" << endl;
               if (_verb)
                  cout << "1-D mesh complete. Mesh Name: "+mesh_name << endl;
               _generator = 1;
               _generated = true;
               RET(0)

            case 106:
               _rita->setEcho(nb_args);
               break;

            case -2:
               break;

            default:
               CHK_MSGB(_cmd->checkEq()==0,_pr,"Unrecognized command: "+_cmd->getToken())
               ret = _rita->_calc->run();
               break;
         }
      }
   }
   RET(ret)
}


int mesh::setRectangle()
{
   string fn="";
   _pr = _PR + "rectangle>";
   double xmin=0., xmax=1., ymin=0., ymax=1.;
   int c[4] = {0,0,0,0}, cv[4] = {0,0,0,0};
   int nb=0, nx=10, ny=10, ret=0, key=0;
   _nb_dof = 1;
   _dim = 2;
   _saved = false;
   _mesh_file = "rita-rectangle.m";
   static const vector<string> kw {"min","max","ne","codes","nbdof"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   for (size_t i=0; i<nb_args; ++i) {
      size_t n = _cmd->getArgs(nb);
      switch (n) {

         case    0:
            if (nb==1) {
               DEF_PAR_R(0,_pr,xmin)
               ymin = xmin;
            }
            if (nb==2) {
               DEF_PAR_R(0,_pr,xmin)
               DEF_PAR_R(1,_pr,ymin)
            }
            else
               MSGR(_pr,"This argument requires 1 or 2 parameters.")
            break;

         case   1:
            if (nb==1) {
               DEF_PAR_R(0,_pr,xmax)
               ymax = xmax;
            }
            if (nb==2) {
               DEF_PAR_R(0,_pr,xmax)
               DEF_PAR_R(1,_pr,ymax)
            }
            else
               MSGR(_pr,"This argument requires 1 or 2 parameters.")
            break;

         case   2:
            if (nb==1) {
               DEF_PAR_R(0,_pr,nx)
               ny = nx;
            }
            else {
               DEF_PAR_R(0,_pr,nx)
               DEF_PAR_R(1,_pr,ny)
            }
            break;

         case   3:
            if (nb==1) {
               DEF_PAR_R(0,_pr,c[0])
               c[1] = c[2] = c[3] = c[0];
            }
            else if (nb==4) {
               DEF_PAR_R(0,_pr,c[0])
               DEF_PAR_R(1,_pr,c[1])
               DEF_PAR_R(2,_pr,c[2])
               DEF_PAR_R(3,_pr,c[3])
            }
            else
               MSGR(_pr,"This argument requires 1 or 4 parameters.")
            break;

         case   4:
            DEF_PAR_R(0,_pr,_nb_dof)
            break;

         case  100:
            cout << "Available arguments: " << Mesh_Rect_help << endl;
            RET(0)

         case  101:
            cout << Mesh_Rect_Help << endl;
            RET(0)

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   if (nb_args>0) {
      CHK_MSGR(xmax<=xmin,_pr,"xmax: "+to_string(xmax)+" must be > than xmin: "+ to_string(xmin))
      CHK_MSGR(ymax<=ymin,_pr,"ymax: "+to_string(ymax)+" must be > ymin: "+to_string(ymin))
      _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,nx,ny,c[3],c[1],c[0],c[2],TRIANGLE,size_t(_nb_dof));
      _data->addMesh(_theMesh,mesh_name);
      if (_verb)
         cout << "2-D mesh complete. Mesh Name: "+mesh_name << endl;
      _generated = true, _generator = 1;
      *_rita->ofh << "  rectangle min=" << xmin << "," << ymin << " max=" << xmax << "," << ymax;
      *_rita->ofh << " ne=" << nx << "," << ny << " codes=" << c[0] << "," << c[1] << "," << c[2];
      *_rita->ofh << "," << c[3] << " nbdof=" << _nb_dof << endl;
      RET(0)
   }

   else {
      if (_verb) {
         cout << "Default rectangle: (0,1)*(0,1)\n";
         cout << "Default node codes: 0\n";
         cout << "Default subdivision: 10*10\n";
         cout << "Default nb of dof: 1" << endl;
      }
      *_rita->ofh << "  rectangle" << endl;

      while (1) {
         READLINE
         switch (key) {

            case   0:
               if (_verb>1)
                  cout << "Setting xmin and ymin ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(2,"Give values of xmin and ymin."),_pr+"min>","Missing xmin and ymin values","",1)
               DEF_PAR_R(-1,_pr+"min>",xmin)
               DEF_PAR_R(-1,_pr+"min>",ymin)
               *_rita->ofh << "    min " << xmin << " " << ymin << endl;
               break;

            case   1:
               if (_verb>1)
                  cout << "Setting xmax and ymax ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(2,"Give values of xmax and ymax."),_pr+"max>","Missing xmax and ymax values","",1)
               DEF_PAR_R(-1,_pr+"max>",xmax)
               DEF_PAR_R(-1,_pr+"max>",ymax)
               *_rita->ofh << "    max " << xmax << " " << ymax << endl;
               break;

            case   2:
               if (_verb>1)
                  cout << "Setting mesh subdivisions ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(2,"Give subdivision in the x and y-directions."),_pr+"ne>",
                                        "Missing subdivisions in x and y directions","",1)
               DEF_PAR_R(-1,_pr+"ne>",nx)
               DEF_PAR_R(-1,_pr+"ne>",ny)
               *_rita->ofh << "    ne " << nx << " " << ny << endl;
               break;

            case   3:
               if (_verb>1)
                  cout << "Setting boundary codes ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(4,"Code to assign on the line y=ymin."),_pr+"codes>","Missing code to assign on the line y=ymin.","",1)
               DEF_PAR_R(-1,_pr+"codes>",c[0])
               DEF_PAR_R(-1,_pr+"codes>",c[1])
               DEF_PAR_R(-1,_pr+"codes>",c[2])
               DEF_PAR_R(-1,_pr+"codes>",c[3])
               *_rita->ofh << "    codes " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << endl;
               if (cv[0]==0)
                  cv[0] = c[0];
               if (cv[1]==0)
                  cv[1] = c[0];
               if (cv[1]==0)
                  cv[1] = c[1];
               if (cv[2]==0)
                  cv[2] = c[1];
               if (cv[2]==0)
                  cv[2] = c[2];
               if (cv[3]==0)
                  cv[3] = c[2];
               if (cv[3]==0)
                  cv[3] = c[3];
               if (cv[0]==0)
                  cv[0] = c[3];
               break;

            case   4:
               if (_verb>1)
                  cout << "Setting number of degrees of freedom ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(1,"Give number of dof per node."),_pr+"nbdof>","Missing number of dof.","",1)
               DEF_PAR_B(-1,_pr+"nbdof>",_nb_dof)
               *_rita->ofh << "    nbdof " << _nb_dof << endl;
               break;

            case   5:
               if (_verb>1)
                  cout << "Saving mesh in file ..." << endl;
               _saved = false;
               CHK_MSG1B(_saved,_pr+"save>","Trying to delete a created mesh.","retype command 'save' to confirm.",1)
               _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,nx,ny,c[3],c[1],c[0],c[2],TRIANGLE,size_t(_nb_dof));
               _saved = true;
               _generator = 2;
               if (_cmd->getNbArgs()>0)
                  ret = _cmd->get(_mesh_file);
               if (!ret) {
                  _theMesh->put(_mesh_file);
                  _data->addMesh(_theMesh,mesh_name);
                  _saved = true;
                  cout << "2-D mesh complete. Mesh Name: "+mesh_name << endl;
               }
               break;

            case 100:
               cout << "Available commands:\n" << Mesh_Rect_help << endl;
               break;

            case 101:
               cout << Mesh_Rect_Help << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               ret = _configure->run();
               _verb = _configure->getVerbose();
               break;

            case 104:
            case 105:
               if (_verb>1)
                  cout << "Getting back to higher level ..." << endl;
               *_rita->ofh << "    end" << endl;
               if (!_saved) {
                  _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,nx,ny,c[3],c[1],c[0],c[2],TRIANGLE,size_t(_nb_dof));
                  _data->addMesh(_theMesh,mesh_name);
                  _generator = 2, _generated = true;
                  _saved = true;
               }
               if (_verb && !_saved)
                  cout << "Mesh 'rectangle' complete. Mesh Name: "+mesh_name << endl;
               RET(90)

            case 106:
               _rita->setEcho(nb_args);
               break;

            case -2:
               break;

            default:
               CHK_MSGB(_cmd->checkEq()==0,_pr,"Unrecognized command: "+_cmd->getToken())
               ret = _rita->_calc->run();
               break;
         }
      }
   }
   RET(ret)
}


int mesh::setCube()
{
   string fn="";
   _pr = _PR + "cube>";
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, ret=0, key=0;
   int nx=10, ny=10, nz=10, cxmin=0, cxmax=0, cymin=0, cymax=0, czmin=0, czmax=0;
   _nb_dof = 1;
   _dim = 3;
   _saved = false;
   if (_verb) {
      cout << "Default cube: (0,1)*(0,1)*(0,1)\n";
      cout << "Default node codes: 0\n";
      cout << "Default subdivision: 10*10*10\n";
      cout << "Default nb of dof: 1" << endl;
   }
   static const vector<string> kw {"min","max","ne","codes","nbdof"};

   _cmd->set(kw);
   size_t nb_args = _cmd->getNbArgs();
   for (size_t i=0; i<nb_args; ++i) {
      size_t n = _cmd->getArgs(nb);
      switch (n) {

         case 100:
            cout << "Available arguments:\n" << Mesh_cube_help << endl;
            RET(0)

         case 101:
            cout << Mesh_cube_Help << endl;
            RET(0)

         case   0:
            if (nb==1) {
               DEF_PAR_R(0,_pr,xmin)
               ymin = zmin = xmin;
            }
            if (nb==3) {
               DEF_PAR_R(0,_pr,xmin)
               DEF_PAR_R(1,_pr,ymin)
               DEF_PAR_R(2,_pr,zmin)
            }
            else
               MSGR(_pr,"This argument requires 1 or 3 parameters.")
            break;

         case   1:
            if (nb==1) {
               DEF_PAR_R(0,_pr,xmax)
               ymax = zmax = xmax;
            }
            if (nb==3) {
               DEF_PAR_R(0,_pr,xmax)
               DEF_PAR_R(1,_pr,ymax)
               DEF_PAR_R(2,_pr,zmax)
            }
            else
               MSGR(_pr,"This argument requires 1 or 3 parameters.")
            break;

         case   2:
            if (nb==1) {
               DEF_PAR_R(0,_pr,nx)
               ny = nz = nx;
            } 
            else {
               DEF_PAR_R(0,_pr,nx)
               DEF_PAR_R(1,_pr,ny)
               DEF_PAR_R(2,_pr,nz)
            }
            break;

         case   3:
            if (nb==1) {
               DEF_PAR_R(0,_pr,cxmin)
               cxmax = cymin = cymax = czmin = czmax = cxmin;
            }
            else if (nb==6) {
               DEF_PAR_R(0,_pr,cxmin)
               DEF_PAR_R(1,_pr,cxmax)
               DEF_PAR_R(2,_pr,cymin)
               DEF_PAR_R(3,_pr,cymax)
               DEF_PAR_R(4,_pr,czmin)
               DEF_PAR_R(5,_pr,czmax)
            }
            else
               MSGR(_pr,"This argument requires 1 or 6 parameters.")
            break;

         case   4:
            DEF_PAR_R(0,_pr,_nb_dof)
            break;

         case   5:
            _mesh_file = _cmd->string_token(0);
            break;

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   if (nb_args>0) {
      CHK_MSGR(xmax<=xmin,_pr,"Value of xmax: "+to_string(xmax)+" must be > xmin: "+to_string(xmin))
      CHK_MSGR(ymax<=ymin,_pr,"Value of ymax: "+to_string(ymax)+" must be > ymin: "+to_string(ymin))
      CHK_MSGR(zmax<=zmin,_pr,"Value of zmax: "+to_string(zmax)+" must be > zmin: "+to_string(zmin))
      if (!_saved) {
         _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,cxmin,cxmax,cymin,
                                    cymax,czmin,czmax,HEXAHEDRON,size_t(_nb_dof));
         _data->addMesh(_theMesh,mesh_name);
         if (_verb)
            cout << "3-D mesh complete. Mesh Name: "+mesh_name << endl;
         _saved = true;
      }
      _generator = 1, _generated = true;
      *_rita->ofh << "  cube min=" << xmin << "," << ymin << "," << zmin << " max=" << xmax
                  << "," << ymax << "," << zmax << " ne=" << nx << "," << ny << "," << nz
                  << " codes=" << cxmin << "," << cxmax << "," << cymin << "," << cymax
                  << "," << czmin << "," << czmax << " nbdof=" << _nb_dof << endl;
   }

   else {
      while (1) {
         *_rita->ofh << "  cube" << endl;
         READLINE
         switch (key) {

            case   0:
               if (_verb>1)
                  cout << "Setting xmin, ymin and zmin ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(3,"Give values of xmin, ymin and zmin."),_pr+"min>","Missing values of xmin, ymin and zmin.","",1)
               DEF_PAR_B(-1,_pr+"min>",xmin)
               DEF_PAR_B(-1,_pr+"min>",ymin)
               DEF_PAR_B(-1,_pr+"min>",zmin)
               *_rita->ofh << "    min " << xmin << " " << ymin << " " << zmin << endl;
               break;

            case   1:
               if (_verb>1)
               cout << "Setting xmax, ymax and zmax ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(3,"Give values of xmax, ymax and zmax."),_pr+"max>","Missing values of xmax, ymax and zmax.","",1)
               DEF_PAR_B(-1,_pr+"max>",xmax)
               DEF_PAR_B(-1,_pr+"max>",ymax)
               DEF_PAR_B(-1,_pr+"max>",zmax)
               *_rita->ofh << "    max " << xmax << " " << ymax << " " << zmax << endl;
               break;

            case   2:
               if (_verb>1)
                  cout << "Setting mesh subdivisions ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(3,"Give subdivision in the x-, y- and z-directions."),_pr+"ne>","Missing subvdivisions in x, y and z directions","",1)
               DEF_PAR_B(-1,_pr+"ne>",nx)
               DEF_PAR_B(-1,_pr+"ne>",ny)
               DEF_PAR_B(-1,_pr+"ne>",nz)
               *_rita->ofh << "    ne " << nx << " " << ny << " " << nz << endl;
               break;
   
            case   3:
               if (_verb>1)
                  cout << "Setting boundary codes ..." << endl;
               CHK_MSG1B(_cmd->setNbArg(6,"Codes to assign to faces."),_pr+"codes>","Missing codes to assign to faces.","",1)
               DEF_PAR_B(-1,_pr+"codes>",cxmin)
               DEF_PAR_B(-1,_pr+"codes>",cxmax)
               DEF_PAR_B(-1,_pr+"codes>",cymin)
               DEF_PAR_B(-1,_pr+"codes>",cymax)
               DEF_PAR_B(-1,_pr+"codes>",czmin)
               DEF_PAR_B(-1,_pr+"codes>",czmax)
               *_rita->ofh << "    codes " << cxmin << " " << cxmax << " " << cymin << " " << cymax
                           << " " << czmin << " " << czmax << endl;
               break;

            case   4:
               if (_verb>1)
                  cout << "Setting number of degrees of freedom ..." << endl;
               CHK_MSGB(_cmd->setNbArg(1,"Give number of dof per node."),_pr+"nbdof>","Missing number of dof.")
               DEF_PAR_B(-1,_pr+"nbdof>",_nb_dof)
               *_rita->ofh << "    nbdof " << _nb_dof << endl;
               break;

            case   5:
               if (_verb>1)
                  cout << "Saving mesh in file ..." << endl;
               if (_saved) {
                  _rita->msg(_pr+"save>","Trying to delete a created mesh.",
                             "You are trying to delete an existing mesh.\n"
                             "retype command 'save' to confirm.");
                  _saved = false;
                  break;
               }
               _mesh_file = "rita-cube.m";
               if (_cmd->getNbArgs()>0)
                  _cmd->get(_mesh_file);
               _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,cxmin,cxmax,cymin,
                                          cymax,czmin,czmax,HEXAHEDRON,size_t(_nb_dof));
               _theMesh->put(_mesh_file);
               _data->addMesh(_theMesh,mesh_name);
               _generator = 3;
               _generated = true;
               cout << "Mesh of cube complete. Mesh Name: M-"+to_string(_data->nbAllMesh()) << endl;
               _saved = true;
               break;

            case 100:
               cout << "Available commands:\n" << Mesh_cube_help << endl;
               break;

            case 101:
               cout << Mesh_cube_Help << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               ret = _configure->run();
               _verb = _configure->getVerbose();
               break;

            case 104:
            case 105:
               if (_verb>1)
                  cout << "Getting back to higher level ..." << endl;
               if (!_saved) {
                  _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,cxmin,
                                             cxmax,cymin,cymax,czmin,czmax,HEXAHEDRON,size_t(_nb_dof));
                  _data->addMesh(_theMesh,mesh_name);
                  _saved = true;
               }
               *_rita->ofh << "    end" << endl;
               if (_verb && !_saved)
                  cout << "Mesh 'cube' complete. Mesh Name: "+mesh_name << endl;
               _generator = 1, _generated = true;
               RET(90)

            case 106:
               _rita->setEcho(nb_args);
               break;

            case -2:
               break;

            default:
               CHK_MSGB(_cmd->checkEq()==0,_pr,"Unrecognized command: "+_cmd->getToken())
               ret = _rita->_calc->run();
               break;
         }
      }
   }
   RET(ret)
}


int mesh::setCode()
{
   _pr = _PR + "code>";
   int nb=0, c=0, d=0, np=0, nc=0, ns=0, nv=0;
   int c_ok=0, points_ok=0, curves_ok=0, surfaces_ok=0, volumes_ok=0;
   vector<int> points, curves, surfaces, volumes;
   CHK_MSGR(_generator>0 && _generator<4,_pr,"Keyword not allowed for generated mesh")
   _generator = 10;
   if (_verb>1) {
      cout << "Default codes for generated nodes:    0" << endl;
      cout << "Default codes for generated elements: 1" << endl;
   }

   static const vector<string> kw {"value","points","curves","surfaces","volumes"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(_pr)
   for (size_t i=0; i<nb_args; ++i) {
      size_t n = _cmd->getArgs(nb);
      switch (n) {
  
         case   0:
            DEF_PAR_R(0,_pr,c)
            if (c<0)
               c = MY_RANDOM - c;
            c_ok++;
            break;

         case   1:
            np = nb;
            for (int j=0; j<np; ++j) {
               DEF_PAR_R(j,_pr,d)
               points.push_back(d);
            }
            points_ok++;
            break;

         case   2:
            nc = nb;
            for (int j=0; j<nc; ++j) {
               DEF_PAR_R(j,_pr,d)
               curves.push_back(d);
            }
            curves_ok++;
            break;

         case   3:
            ns = nb;
            for (int j=0; j<ns; ++j) {
               DEF_PAR_R(j,_pr,d)
               surfaces.push_back(d);
            }
            surfaces_ok++;
            break;

         case   4:
            nv = nb;
            for (int j=0; j<nv; ++j) {
               DEF_PAR_R(j,_pr,d)
               volumes.push_back(d);
            }
            volumes_ok++;
            break;

         case 100:
            cout << "Available arguments: " << Mesh_code_help << endl;
            RET(0)

         case 101:
            cout << Mesh_code_Help << endl;
            RET(0)

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   CHK_MSGR(!c_ok,_pr,"Missing code value.")
   CHK_MSGR(points_ok>1 || curves_ok>1 || surfaces_ok>1 || volumes_ok>1,_pr,"Each entity must be given only once.")
   CHK_MSGR(np+nc+ns+nv==0,_pr,"At least one entity must be given.")
   *_rita->ofh << "  code value=" << c;
   if (points_ok) {
      *_rita->ofh << " points=";
      for (int j=0; j<np; ++j)
         _Pcode[c].l.push_back(points[j]), _Pcode[c].nb++;
      for (int j=0; j<np-1; ++j)
         *_rita->ofh << points[j] << ",";
      *_rita->ofh << points[np-1] << endl;
   }
   if (curves_ok) {
      *_rita->ofh << " curves=";
      for (int j=0; j<nc; ++j)
         _Ccode[c].l.push_back(curves[j]), _Ccode[c].nb++;
      for (int j=0; j<nc-1; ++j)
         *_rita->ofh << curves[j] << ",";
      *_rita->ofh << curves[nc-1] << endl;
   }
   if (surfaces_ok) {
      *_rita->ofh << " surfaces=";
      for (int j=0; j<ns; ++j)
         _Scode[c].l.push_back(surfaces[j]), _Scode[c].nb++;
      for (int j=0; j<ns-1; ++j)
         *_rita->ofh << surfaces[j] << ",";
      *_rita->ofh << surfaces[ns-1] << endl;
   }
   if (volumes_ok) {
      *_rita->ofh << " volumes=";
      for (int j=0; j<nv; ++j)
         _Vcode[c].l.push_back(volumes[j]), _Vcode[c].nb++;
      for (int j=0; j<nv-1; j++)
         *_rita->ofh << volumes[j] << ",";
      *_rita->ofh << volumes[nv-1] << endl;
   }
   RET(0)
}


int mesh::setPoint()
{
   _pr = _PR + "point>";
   int nb=0;
   bool n_ok=false, x_ok=false, y_ok=false, z_ok=false, h_ok=false, del_ok=false;
   static int it = 0;
   if (it==0) {
      _point.n = 1;
      _point.x = _point.y = _point.z = 0.;
      _point.h = 0.1;
   }
   it++;
   static const vector<string> kw {"label","n","coord","size"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(_pr)
   for (size_t i=0; i<nb_args; ++i) {
      size_t n = _cmd->getArgs(nb);
      switch (n) {

         case   0:
         case   1:
            NO_VALUE_ARG(_pr)
            DEF_PAR_R(0,_pr,_point.n)
            n_ok = true;
            break;

         case   2:
            NO_VALUE_ARG(_pr)
            _point.x = _point.y = _point.z = 0.;
            DEF_PAR_R(0,_pr,_point.x)
            x_ok = true;
            if (nb>1) {
               DEF_PAR_R(1,_pr,_point.y)
               y_ok = true;
            }
            if (nb>2) {
               DEF_PAR_R(2,_pr,_point.z)
               z_ok = true;
            }
            break;

         case   3:
            DEF_PAR_R(0,_pr,_point.h)
            h_ok = true;
            break;

         case 100:
            cout << "Available arguments: " << Mesh_point_help << endl;
            RET(0)

         case 101:
            cout << Mesh_point_help << endl;
            RET(0)

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   if (del_ok) {
      CHK_MSGR(!n_ok,_pr,"No point label given.")
      _points.erase(_point.n);
      *_rita->ofh << "  point  label="<< _point.n << "  delete" << endl;
      RET(0)
   }
   if (!n_ok) {
      _point.n++;
      if (_verb>1)
         cout << "Assumed point label: " << _point.n << endl;
   }
   CHK_MSGR(_point.n<=0,_pr,"Label must be positive.")
   if (!x_ok && _verb>1)
      cout << "Warning: No x-coordinate given. Assumed x-coordinate: " << _point.x << endl;
   if (!y_ok && _verb>1)
      cout << "Warning: No y-coordinate given. Assumed y-coordinate: " << _point.y << endl;
   if (!z_ok && _verb>1)
      cout << "Warning: No z-coordinate given. Assumed z-coordinate: " << _point.z << endl;
   if (!h_ok && _verb)
      cout << "Warning: No mesh size given. Assumed mesh size: " << _point.h << endl;
   if (_generator>0 && _generator<=3) {
      _rita->msg(_pr,"Mesh already generated. Needs retyping command to confirm.");
      _generator = 10;
   }
   _points[_point.n] = _point;
   *_rita->ofh << "  point  label=" << _point.n << "  x=" << _point.x << "  y=" << _point.y
               << "  z=" << _point.z << "  size=" << _point.h << endl;
#ifndef USE_GMSH
   _theDomain->insertVertex(_point.x,_point.y,_point.z,_point.h,0);
#endif
   RET(0)
}


int mesh::setCurve()
{
   _pr = _PR + "curve>";
   static int nn = 0;
   int nb=0, n1=0, n2=0, n3=0;
   bool n_ok=false, line_ok=false, circle_ok=false, del_ok=false;
   if (_generator>0 && _generator<=3) {
      _generator = 0;
      MSGR(_pr,"Mesh already generated. Needs retyping command to confirm.")
   }
   _generator = 10;
   static const vector<string> kw {"label","n","line","circle","del$ete"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(_pr)
   for (size_t i=0; i<nb_args; ++i) {
      size_t n = _cmd->getArgs(nb);
      switch (n) {

         case   0:
         case   1:
            NO_VALUE_ARG(_pr)
            DEF_PAR_R(0,_pr,nn)
            n_ok = true;
            break;

         case   2:
            NO_VALUE_ARG(_pr)
            CHK_MSGR(nb!=2,_pr,"Illegal number of arguments for a line")
            DEF_PAR_R(0,_pr,n1)
            DEF_PAR_R(1,_pr,n2)
            line_ok = true;
            break;

         case   3:
            NO_VALUE_ARG(_pr)
            CHK_MSGR(nb!=3,_pr,"Illegal number of arguments for a circle")
            DEF_PAR_R(0,_pr,n1)
            DEF_PAR_R(1,_pr,n2)
            DEF_PAR_R(2,_pr,n3)
            circle_ok = true;
            break;

         case   4:
            del_ok = true;
            break;

         case 100:
            cout << "Available arguments: " << Mesh_curve_help << endl;
            RET(0)

         case 101:
            cout << Mesh_curve_Help << endl;
            RET(0)

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }

   if (del_ok) {
      CHK_MSGR(!n_ok,_pr,"No curve label given.")
      _curve.erase(nn);
      *_rita->ofh << "  curve label="<< nn << "  delete" << endl;
      RET(0)
   }
   CHK_MSGR(nn<=0,_pr,"Label must be positive.")
   if (!n_ok) {
      nn++;
      if (_verb>1)
         cout << "Assumed curve label: " << nn << endl;
   }
   CHK_MSGR(!line_ok && !circle_ok,_pr,"A line or a circle must be defined.")
   CHK_MSGR(line_ok && circle_ok,_pr,"Curve cannot be defined as a line and a circle simultaneously.")
   _curve[nn].nb = nn;
   *_rita->ofh << "  curve label=" << nn;
   if (line_ok) {
      CHK_MSGR(_points.find(n1)==_points.end(),_pr,"Undefined end point: "+to_string(n1))
      CHK_MSGR(_points.find(n2)==_points.end(),_pr,"Undefined end point: "+to_string(n2))
      _curve[nn].l.clear();
      _curve[nn].l.push_back(n1), _curve[nn].l.push_back(n2);
      _curve[nn].nb = 2;
      _curve[nn].type = 1;
      *_rita->ofh << "  line=" << n1 << "," << n2 << endl;
   }
   else if (circle_ok) {
      CHK_MSGR(_points.find(n1)==_points.end(),_pr,"Undefined point: "+to_string(n1))
      CHK_MSGR(_points.find(n2)==_points.end(),_pr,"Undefined point: "+to_string(n2))
      CHK_MSGR(_points.find(n3)==_points.end(),_pr,"Undefined point: "+to_string(n3))
      _curve[nn].l.clear();
      _curve[nn].l.push_back(n1);
      _curve[nn].l.push_back(n2);
      _curve[nn].l.push_back(n3);
      _curve[nn].nb = 3;
      _curve[nn].type = 2;
      *_rita->ofh << "  circle=" << n1 << "," << n2 << "," << n3 << endl;
   }
#ifndef USE_GMSH
   _theDomain->insertLine(_curve[nn].l[0],_curve[nn].l[1],0);
#endif
   RET(0)
}


int mesh::setContour()
{
   _pr = _PR + "contour>";
   int nb=0, nn=0, s=1, n1=0, n2=0, n=0;
   vector<int> curv, surf;
   int n_ok=0, curves_ok=0, surfaces_ok=0, del_ok=0;
   if (_generator>0 && _generator<=3) {
      _generator = 0;
      MSGR(_pr,"Mesh already generated")
   }
   _generator = 10;
   static const vector<string> kw {"label","n","curv$es","surf$aces"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(_pr)
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case   0:
         case   1:
            NO_VALUE_ARG(_pr)
            DEF_PAR_R(0,_pr,nn)
            n_ok++;
            break;

         case   2:
            NO_VALUE_ARG(_pr)
            for (int i=0; i<nb; i++) {
               DEF_PAR_R(i,_pr,n)
               curv.push_back(n);
            }
            curves_ok++;
            break;

         case   3:
            NO_VALUE_ARG(_pr)
            for (int j=0; j<nb; j++) {
               DEF_PAR_R(j,_pr,nn)
               surf.push_back(nn);
            }
            surfaces_ok++;
            break;

         case   4:
            del_ok++;
            break;

         case 100:
            cout << "Available arguments: " << Mesh_contour_help << endl;
            RET(0)

         case 101:
            cout << Mesh_contour_Help << endl;
            RET(0)

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   CHK_MSGR(!n_ok,_pr,"Missing contour label.")
   CHK_MSGR(n_ok>1,_pr,"Contour label cannot be given more than once.")
   if (del_ok) {
      CHK_MSGR(del_ok>1,_pr,"Keyword delete cannot be given more than once.")
      _Ccontour.erase(nn);
      _Scontour.erase(nn);
      *_rita->ofh << "  contour  label="<< nn << "  delete" << endl;
      RET(0)
   }
   if (nn==0) {
      cout << "Contour label cannot be zero." << endl;
      MSGR(_pr,"Contour label cannot be zero.")
   }
   CHK_MSGR(!curves_ok && !surfaces_ok,_pr,"Contour must be defined either by curves or surfaces.")
   CHK_MSGR(curves_ok && surfaces_ok,_pr,"Contour cannot be defined by curves or surfaces simultaneously.")
   s = 1;
   if (nn<0)
      s = -1, nn = -nn;
   if (curves_ok) {
      CHK_MSGR(curves_ok>1,_pr,"A curve contour has already been defined.")
      _Ccontour[nn].nb = 0;
      for (int i=0; i<nb; ++i) {
         int l = curv[i];
         CHK_MSGR(_curve.find(l)==_curve.end(),_pr,"Undefined curve: "+to_string(l))
         _Ccontour[nn].l.push_back(l), _Ccontour[nn].nb++;
      }
      n1 = _curve[_Ccontour[nn].l[0]].l[0];
      if (s==-1)
         n1 = _curve[_Ccontour[nn].l[0]].l[1];
      n2 = _curve[_Ccontour[nn].l[_Ccontour[nn].nb-1]].l[1];
      if (s==-1)
         n2 = _curve[_Ccontour[nn].l[_Ccontour[nn].nb-1]].l[0];
      CHK_MSGR(n1!=n2,_pr,"Non closed contour.")
      *_rita->ofh << "  contour  label=" << nn << "  curves=";
      for (int i=0; i<nb-1; ++i)
         *_rita->ofh << curv[i] << ",";
      *_rita->ofh << curv[nb-1] << endl;
      _nb_Ccontour = _Ccontour.size();
   }
   if (surfaces_ok) {
      CHK_MSGR(surfaces_ok>1,_pr,"A surface contour already defined.")
      _Scontour[nn].nb = 0;
      for (int i=0; i<nb; ++i) {
         int l = surf[i];
         CHK_MSGR(_surface.find(l)==_surface.end(),_pr,"Undefined surface: "+to_string(l))
         _Scontour[nn].l.push_back(l); _Scontour[nn].nb++;
      }
      n1 = _curve[_Scontour[nn].l[0]].l[0];
      if (s==-1)
         n1 = _surface[_Scontour[nn].l[0]].l[1];
      n2 = _surface[_Scontour[nn].l[_Scontour[nn].nb-1]].l[1];
      if (s==-1)
         n2 = _surface[_Scontour[nn].l[_Scontour[nn].nb-1]].l[0];
      CHK_MSGR(n1!=n2,_pr,"Non closed contour.")
      *_rita->ofh << "  contour  label=" << nn << "  surfaces=";
      for (int i=0; i<nb-1; ++i)
         *_rita->ofh << surf[i] << ",";
      *_rita->ofh << curv[nb-1] << endl;
      _nb_Scontour = _Scontour.size();
   }
   RET(0)
}


int mesh::setSurface()
{
   _pr = _PR + "surface>";
   static int nn = 1;
   int m=0, nb=0, n=_surface.size()+1;
   bool n_ok=false, cont_ok=false, del_ok=false;
   vector<int> cont;
   if (_generator>0 && _generator<=3) {
      _generator = 0;
      MSGR(_pr,"Mesh already generated. Needs retyping command to confirm.")
   }
   _generator = 10;
   static const vector<string> kw {"label","n","contours"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(_pr)
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case   0:
         case   1:
            DEF_PAR_R(0,_pr,m)
            n_ok = true;
            break;

         case   2:
            NO_VALUE_ARG(_pr)
            for (int j=0; j<nb; ++j) {
               DEF_PAR_R(j,_pr,m)
               cont.push_back(m);
            }
            cont_ok = true;
            break;

         case   3:
            del_ok = true;
            break;

         case 100:
            cout << "Available arguments: " << Mesh_surf_help << endl;
            RET(0)

         case 101:
            cout << Mesh_surf_Help << endl;
            RET(0)

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   if (del_ok) {
      CHK_MSGR(!n_ok,_pr,"No surface label given.")
      _surface.erase(nn);
      *_rita->ofh << "  surface  label="<< nn << "  delete" << endl;
      RET(0)
   }
   if (!n_ok) {
      nn++;
      if (_verb)
         cout << "Assumed surface label: " << nn << endl;
   }
   CHK_MSGR(nn==0,_pr,"Label 0 is forbidden")
   CHK_MSGR(!cont_ok,_pr,"Surface contours must be given.")
   _surface[nn].nb = 0;
   for (auto const& c: cont) {
      CHK_MSGR(_Ccontour.find(c)==_Ccontour.end(),_pr,"Undefined curve contour: "+to_string(n))
      _surface[nn].l.push_back(c); _surface[nn].nb++;
   }
   *_rita->ofh << "  surface  label=" << nn << "  contours=";
   for (int i=0; i<nb-1; ++i)
      *_rita->ofh << cont[i] << ",";
   *_rita->ofh << cont[nb-1] << endl;
   _nb_surface = _surface.size();
   RET(0)
}


int mesh::setVolume()
{
   _pr = _PR + "volume>";
   _generator = 10;
   _nb_volume = 0;
   MSGR(_pr,"Volume generation not implemented yet.")
}


int mesh::setSubDomain()
{
   _pr = _PR + "subdomain>";
   string fn="";
   size_t nb_args = 0;
   if (_generator>0 && _generator<=3) {
      _rita->msg(_pr,"Mesh already generated. Needs retyping command to confirm.");
      _generator = 0;
      RET(0)
   }
   _generator = 10;
   int ret=0;
   *_rita->ofh << "   subdomain " << endl;   
   _sd.ln = 1, _sd.orientation = 1, _sd.code = 1;
   if (_verb) {
      cout << "Default line in subdomain: 1\n";
      cout << "Default orientation: 1\n";
      cout << "Default code: 10" << endl;
   }

   static const vector<string> kw {"line","orient$ation"};
   int key=0;
   while (1) {
      READLINE
      switch (key) {

         case   0:
            if (_verb>1)
               cout << "Setting line ..." << endl;
            CHK_MSG1B(_cmd->setNbArg(1,"Give label of a line in subdomain."),_pr+"line>","Missing label of vertex in subdomain.","",1)
            ret = _cmd->get(_sd.ln);
            if (!ret)
               *_rita->ofh << "    line " << _sd.ln << endl;
            break;

         case   1:
            if (_verb>1)
               cout << "Setting orientation ..." << endl;
            CHK_MSG1B(_cmd->setNbArg(1,"Give orientation of subdomain (1/-1)."),_pr+"orientation>","Missing orientation.","",1)
            ret = _cmd->get(_sd.orientation);
            if (!ret)
               *_rita->ofh << "    orientation " << _sd.orientation << endl;   
            break;

         case   2:
            if (_verb>1)
               cout << "Setting code ..." << endl;
            CHK_MSG1B(_cmd->setNbArg(1,"Give code to associate to subdomain"),_pr+"code>","Missing code to associate to subdomain.","",1)
            ret = _cmd->get(_sd.code);
            if (!ret)
               *_rita->ofh << "    code " << _sd.code << endl;   
            break;

         case   3:
            if (_verb>1)
               cout << "Saving subdomain ..." << endl;
            _cmd->setNbArg(0);
    //        _theDomain->insertSubDomain(_sd.ln,_sd.orientation,_sd.code);
            _nb_sub_domain++;
            _subdomains.push_back(_sd);
            cout << "Subdomain Added." << endl;
            cout << "Line " << _sd.ln << ", Orientation " << _sd.orientation
                 << ", Code " << _sd.code << endl;
            cout << "Total number of subdomains: " << _nb_sub_domain << endl;
            *_rita->ofh << "    save" << endl;
            _saved = true;
            RET(0)

         case   4:
         case   5:
            if (_verb>1)
               cout << "Getting back to mesh menu ..." << endl;
            _cmd->setNbArg(0);
            CHK_MSG1B(!_saved,_pr,"Subdomain not saved.","Type again 'subdomain' and save generated one",1)
            *_rita->ofh << "    end" << endl;
            RET(0)

         case 100:
            cout << "Available arguments: " << Mesh_subdomain_help << endl;
            break;

         case 101:
            cout << Mesh_subdomain_Help << endl;
            RET(0)

         case 102:
            _rita->getLicense();
            RET(0)

         case 103:
            ret = _configure->run();
            _verb = _configure->getVerbose();
            break;

         case 104:
         case 105:
            RET(100)

         case 106:
            _rita->setEcho(nb_args);
            break;

         case -2:
         case -3:
         case -4:
            break;

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
       }
   }
   RET(0)
}


void mesh::saveGeo(const string& file)
{
   _generator = 10;
   ofstream s(file.c_str());
   for (auto const& v: _points)
      s << "Point(" << v.first << ") = {" << v.second.x << ", " << v.second.y
        << ", " << v.second.z << ", " << v.second.h << "};" << endl;

   for (auto const& v: _Pcode) {
      s << "Physical Point(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }
                   
   for (auto const& v: _curve) {
      Entity curve = v.second;
      if (curve.type==1)
         s << "Line(" << v.first << ") = {" << curve.l[0] << ", " << curve.l[1] << "};" << endl;
      else if (curve.type==2)
         s << "Circle(" << v.first << ") = {" << curve.l[0] << ", " << curve.l[2] << ", " 
           << curve.l[1] << "};" << endl;
   }

   for (auto const& v: _Ccode) {
      Entity code = v.second;
      s << "Physical Curve(" << v.first << ") = {";
      for (int i=0; i<code.nb-1; ++i)
         s << code.l[i] << ", ";
      s << code.l[code.nb-1] << "};" << endl;
   }

   for (auto const& v: _Ccontour) {
      s << "Curve Loop(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }

   for (auto const& v: _surface) {
      s << "Plane Surface(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }

   for (auto const& v: _Scode) {
      s << "Physical Surface(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }

   for (auto const& v: _Vcode) {
      s << "Physical Volume(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }
   _geo = true;
}


void mesh::saveDomain(const string& file)
{
   ofstream id(file.c_str());
   id << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n<OFELI_File>\n<info>" << endl;
   id << "<title>Domain file created by rita</title>" << endl;
   id << "<date></date>\n<author></author>\n</info>" << endl;
   id << "<Project name=\"\">" << endl;
   id << "  <Domain dim=\"2\">" << endl;
   for (auto const& p: _points)
      id << "    <vertex>" << p.second.x << " " << p.second.y << " " << p.second.z << " " << p.second.h << "</vertex>" << endl;
   for (int i=1; i<=_nb_curve; ++i)
      id << "    <line>" << _curve[i].l[0] << " " << _curve[i].l[1] << " "
         << 0 << "</line>" << endl;
/*   for (int i=0; i<_nb_circle; ++i)
      id << "    <circle>" << _circles[i].pt1 << " " << _circles[i].c << " " << _circles[i].pt2
         << " " << 0 << "</circle>" << endl;*/
   for (auto const& s: _subdomains)
      id << "    <Subdomain>" << s.ln << " " << s.orientation << " " << s.code << "</Subdomain>" << endl;
   for (int i=0; i<_nb_sub_domain; ++i)
      id << "    <Subdomain>" << _subdomains[i].ln << " " << _subdomains[i].orientation << " "
         << _subdomains[i].code << "</Subdomain>" << endl;
   id << "  </Domain>\n</Project>\n</OFELI_File>" << endl;
}


int mesh::Generate()
{
   _pr = _PR + "generate>";
   if (_generator>0 && _generator<=3) {
      CHK_MSGR(_generator==1,_pr,"A 1-D mesh has already been generated.")
      CHK_MSGR(_generator==2,_pr,"A 2-D mesh has already been generated.")
      CHK_MSGR(_generator==3,_pr,"A 3-D mesh has already been generated.")
   }
   _generator = 10;
   if (_theMesh!=nullptr)
      delete _theMesh, _theMesh = nullptr;
   _mesh_file = "rita.m";
   
// Save geo gmsh file and generate gmsh mesh
   if (!_geo)
      saveGeo("rita.geo");
#ifdef USE_GMSH
   if (_verb)
      cout << "Starting mesh generation using Gmsh ..." << endl;
   gmsh::initialize();
   gmsh::option::setNumber("General.Terminal",1);
   gmsh::model::add("rita");
   for (auto const& v: _points)
      gmsh::model::geo::addPoint(v.second.x,v.second.y,v.second.z,v.second.h,v.first);
   for (auto const& v: _curve) {
      if (v.second.type==1)
         gmsh::model::geo::addLine(v.second.l[0],v.second.l[1],v.first);
      else if (v.second.type==2)
         gmsh::model::geo::addCircleArc(v.second.l[0],v.second.l[2],v.second.l[1],v.first);
   }
   for (auto const& v: _Ccontour)
      gmsh::model::geo::addCurveLoop(v.second.l,v.first);
   for (auto const& v: _surface)
      gmsh::model::geo::addPlaneSurface(v.second.l,v.first);
   for (auto const& v: _Pcode)
      gmsh::model::addPhysicalGroup(0,v.second.l,v.first);
   for (auto const& v: _Ccode)
      gmsh::model::addPhysicalGroup(1,v.second.l,v.first);
   for (auto const& v: _Scode)
      gmsh::model::addPhysicalGroup(2,v.second.l,v.first);
   for (auto const& v: _Vcode)
      gmsh::model::addPhysicalGroup(3,v.second.l,v.first);
   gmsh::model::setPhysicalName(2,4,"Domain");
   gmsh::model::geo::synchronize();
   gmsh::model::mesh::generate(2);
   gmsh::write("rita.msh");
   gmsh::finalize();
   if (_verb)
      cout << "Gmsh mesh generation complete." << endl;
   _theMesh = new OFELI::Mesh("rita.msh",false,NODE_DOF,_nb_dof);
   _generated = true;
#else
   _theDomain->setDim(_dim);
   _theDomain->setNbDOF(size_t(_nb_dof));
   saveDomain("rita.dom");
   _theDomain->genMesh(_mesh_file);
   _theMesh = new OFELI::Mesh(_mesh_file,false,NODE_DOF,2);
   _generated = true;
   _generator = 4;
#endif
   *_rita->ofh << "  generate" << endl;
      if (mesh_name=="")
         mesh_name = "M-" + to_string(_data->nbAllMesh()+1);
   _data->addMesh(_theMesh,mesh_name);
   if (_verb)
      cout << "Mesh complete. Mesh Name: "+mesh_name << endl;
   RET(0)
}


int mesh::setNbDOF()
{
   int n = 0;
   _pr = _PR;
   CHK_MSG1R(_cmd->setNbArg(1,"Give number of degrees of freedom per node."),_pr+"nbdof>",
             "Missing number of degrees of freedom.","",1)
   int ret = _cmd->get(n);
   if (!ret) {
      *_rita->ofh << "    nbdof " << n << endl;
      _nb_dof = n;
   }
   RET(ret)
}


int mesh::Plot()
{
   _pr = _PR + "plot>";
   string fn="";
   CHK_MSGR(_theMesh==nullptr,_pr,"No mesh to plot.")
   _cmd->setNbArg(0);
   string file = "rita.msh";
   _theMesh->save(file);
   string com = "gmsh " + file;
   CHK_MSGR(system(com.c_str()),_pr,"Unrecognizable system command.")
   _data->temp_file.push_back(file);
   RET(0)
}


int mesh::Read()
{
   _pr = _PR + "read>";
   int ret=0, key=0;
   ifstream ip;
   string dom_file, file, bamg_file, geo_file, out_file, Cmd, msh_file, fn="";
   int mesh_ok=0, geo_ok=0, gmsh_ok=0;
   static const vector<string> kw {"mesh","geo","gmsh"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   for (size_t i=0; i<nb_args; ++i) {
      size_t n = _cmd->getArg();
      switch (n) {

         case 100:
            cout << "Available arguments: " << Mesh_read_help << endl;
            RET(0)

         case 110:
            cout << Mesh_read_Help << endl;
            RET(0)

         case   0:
            file = _cmd->string_token();
            mesh_ok++;
            break;

         case   1:
            file = _cmd->string_token();
            geo_ok++;
            break;

         case   2:
            file = _cmd->string_token();
            gmsh_ok++;
            break;

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   if (nb_args>0) {
      CHK_MSGR(mesh_ok+geo_ok+gmsh_ok==0,_pr,"No input file provided.")
      CHK_MSGR(mesh_ok+geo_ok+gmsh_ok>1,_pr,"Only one input file must be provided.")
      ip.open(file);
      if (ip.is_open())
         ip.close();
      else
         MSGR(_pr,"Unable to open file: "+file)
      *_rita->ofh << "  read";
      if (mesh_ok) {
         CHK_MSGR(file.substr(file.find_last_of(".")+1)!="m",_pr,"File extension must be \".m\"")
         _theMesh = new OFELI::Mesh;
         _theMesh->get(file);
         CHK_MSGR(_theMesh->getNbNodes()==0,_pr,"Empty mesh")
         _mesh_file = file;
         _nb_dof = _theMesh->getNbDOF() / _theMesh->getNbNodes();
         *_rita->ofh << "  mesh=" << file;
         if (mesh_name=="")
            mesh_name = "M-1" + to_string(_data->nbAllMesh()+1);
         _data->addMesh(_theMesh,mesh_name);
         _generator = 1, _generated = true;
      }
      else if (geo_ok) {
         CHK_MSGR(file.substr(file.find_last_of(".")+1)!="geo",_pr,"File extension must be \".geo\"")
         msh_file = file.substr(0,file.rfind(".")) + "msh";
         Cmd = "gmsh -2 " + file + " -o " + msh_file;
         CHK_MSGR(system(Cmd.c_str()),_pr,"Unrecognizable system command.")
         _theMesh = new OFELI::Mesh;
         _theMesh->get(file,GMSH);
         _data->MeshName.push_back("M"+to_string(_data->theMesh.size()));
         _data->theMesh.push_back(_theMesh);
         _generated = true;
         *_rita->ofh << "  geo=" << file;
      }
      else if (gmsh_ok) {
         CHK_MSGR(file.substr(file.find_last_of(".")+1)!="msh",_pr,"File extension must be \".msh\"")
         _theMesh = new OFELI::Mesh;
         _theMesh->get(file,GMSH);
         if (mesh_name=="")
            mesh_name = "M-" + to_string(_data->nbAllMesh()+1);
         _data->addMesh(_theMesh,mesh_name);
         _generator = 1;
         _generated = true;
         *_rita->ofh << "  gmsh=" << file;
      }
      *_rita->ofh << endl;
   }
   else {
      _cmd->setNbArg(0);
      while (1) {
         READLINE
         switch (key) {

            case   0:
               CHK_MSG1B(_cmd->setNbArg(1,"Give OFELI mesh file name."),_pr+"mesh>","Missing Mesh file name.","",1)
               ret = _cmd->get(file);
               ip.open(file);
               if (ip.is_open()) {
                  ip.close();
                  _theMesh = new OFELI::Mesh(file);
                  *_rita->ofh << "  read mesh " << file << endl;
                  CHK_MSGR(_theMesh->getNbNodes()==0,_pr+"mesh>","Empty mesh")
                  _nb_dof = _theMesh->getNbDOF() / _theMesh->getNbNodes();
                  if (mesh_name=="")
                     mesh_name = "M-" + to_string(_data->nbAllMesh()+1);
                  _data->addMesh(_theMesh,mesh_name);
                  _generator = 1;
                  _generated = true;
                  ret = 90;
               }
               else {
                  _rita->msg(_pr+"mesh>","Unable to open file: "+file);
                  _generated = false;
                  ret = 0;
               }
               RET(ret)

            case   1:
               CHK_MSG1B(_cmd->setNbArg(1,"Give geo gmsh file name."),_pr+"geo>","Missing geo file name.","",1)
               ret = _cmd->get(file);
               if (!ret) {
                  CHK_MSGB(file.substr(file.find_last_of(".")+1)!=".geo",_pr+"geo>",
                           "File extension must be \".geo\"")
                  msh_file = file.substr(0,file.rfind(".")) + "msh";
                  Cmd = "gmsh -2 " + file + " -o " + msh_file;
                  CHK_MSGB(system(Cmd.c_str()),_pr+"geo>","Unrecognizable system command.")
                  _theMesh = new OFELI::Mesh;
                  _theMesh->get(msh_file,GMSH);
                  *_rita->ofh << "  read geo " << file << endl;
                  _data->MeshName.push_back("M"+to_string(_data->theMesh.size()));
                  _data->theMesh.push_back(_theMesh);
                  _geo = true;
                  _generated = false;
               }
               RET(90)

            case 2:
               CHK_MSG1B(_cmd->setNbArg(1,"Give gmsh file name where to save mesh."),_pr+"gmsh>","Missing gmsh file name.","",1)
               ret = _cmd->get(file);
               if (!ret) {
                  CHK_MSGB(file.substr(file.find_last_of(".")+1)!=".msh",_pr+"gmsh>","File extension must be \".msh\"")
                  _theMesh = new OFELI::Mesh;
                  _theMesh->get(file,GMSH);
                  if (mesh_name=="")
                     mesh_name = "M-" + to_string(_data->nbAllMesh()+1);
                  _data->addMesh(_theMesh,mesh_name);
                  *_rita->ofh << "  read gmsh " << file << endl;
               }
               _generated = true, _generator = 1;
               RET(90)

            case 100:
               cout << "Available commands: " << Mesh_read_help << endl;
               break;
            
            case 101:
               cout << Mesh_read_Help << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               ret = _configure->run();
               _verb = _configure->getVerbose();
               break;

            case 104:
               CHK_MSG1B(_cmd->setNbArg(1,"Data name to be given.",1),_pr+"end>","Missing data name.","",1)
               if (!_cmd->get(fn))
                  _data->print(fn);
               break;

            case 105:
               _data->Summary();
               break;

            case 106:
               _rita->setEcho(nb_args);
               break;

            case -2:
            case -3:
            case -4:
               break;

            default:
               MSGR(_pr,"Unknown argument: "+_cmd->getToken())
         }
      }
   }
   RET(0)
}


int mesh::Save()
{
   _pr = _PR + "save>";
   string domain_f="rita.dom", geo_f="rita.geo", mesh_f="rita.m", gmsh_f="rita.msh";
   string t="", vtk_f="rita.vtk", gnuplot_f="rita-gpl.dat", matlab_f="rita-matlab.m";
   string tecplot_f="rita-tecplot.dat";
   int domain_ok=0, geo_ok=0, mesh_ok=0, gmsh_ok=0, vtk_ok=0, gnuplot_ok=0, matlab_ok=0, tecplot_ok=0;
   int ret=0;
   static const vector<string> kw {"domain","geo","mesh","gmsh","vtk","gnuplot","matlab","tecplot"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   for (size_t i=0; i<nb_args; ++i) {
      size_t n = _cmd->getArg();
      switch (n) {

         case 0:
#ifndef USE_GMSH
            if ((t=_cmd->string_token())!="")
               domain_f = t;
            domain_ok++;
#endif
            break;

         case 1:
            if ((t=_cmd->string_token())!="")
               geo_f = t;
            geo_ok++;
            break;

         case 2:
            if ((t=_cmd->string_token())!="")
               mesh_f = t;
            mesh_ok++;
            break;

         case 3:
            if ((t=_cmd->string_token())!="")
               gmsh_f = t;
            gmsh_ok++;
            break;

         case 4:
            if ((t=_cmd->string_token())!="")
               vtk_f = t;
            vtk_ok++;
            break;

         case 5:
            if ((t=_cmd->string_token())!="")
               gnuplot_f = t;
            gnuplot_ok++;
            break;

         case 6:
            if ((t=_cmd->string_token())!="")
               matlab_f = t;
            matlab_ok++;
            break;

         case 7:
            if ((t=_cmd->string_token())!="")
               tecplot_f = t;
            tecplot_ok++;
            break;

         case 100:
            cout << "Available arguments: " << Mesh_save_help << endl;
            RET(0)

         case 101:
            cout << Mesh_save_Help << endl;
            RET(0)

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   CHK_MSGR(domain_ok+geo_ok+mesh_ok+gmsh_ok+vtk_ok+gnuplot_ok+matlab_ok+tecplot_ok==0,_pr,"Nothing to save.")
   CHK_MSGR(mesh_ok+gmsh_ok+vtk_ok+gnuplot_ok+matlab_ok+tecplot_ok>0 && !_generated,_pr,"No generated mesh to be saved")
   *_rita->ofh << "  save";
   if (domain_ok) {
      saveDomain(domain_f);
      *_rita->ofh << "  domain=" << domain_f;
   }
   if (geo_ok) {
      saveGeo(geo_f);
      *_rita->ofh << "  geo=" << geo_f;
   }
   if (mesh_ok) {
      _mesh_file = mesh_f;
      _theMesh->put(mesh_f);
      *_rita->ofh << "  mesh=" << mesh_f;
   }
   if (gmsh_ok) {
      saveMesh(gmsh_f,*_theMesh,GMSH);
      *_rita->ofh << "  gmsh=" << gmsh_f;
   }
   if (vtk_ok) {
      saveMesh(vtk_f,*_theMesh,VTK);
      *_rita->ofh << "  vtk=" << vtk_f;
   }
   if (gnuplot_ok) {
      saveMesh(gnuplot_f,*_theMesh,GNUPLOT);
      *_rita->ofh << "  gnuplot=" << gnuplot_f;
   }
   if (matlab_ok) {
      saveMesh(matlab_f,*_theMesh,MATLAB);
      *_rita->ofh << "  matlab=" << matlab_f;
   }
   if (tecplot_ok) {
      saveMesh(tecplot_f,*_theMesh,TECPLOT);
      *_rita->ofh << "  tecplot=" << tecplot_f;
   }
   *_rita->ofh << endl;
   RET(ret)
}

} /* namespace RITA */
