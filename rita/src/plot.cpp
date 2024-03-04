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

                    Implementation of function 'plot' in 'data'

  ==============================================================================*/

#include "data.h"
#include "configure.h"
#include "rita.h"
#include "linear_algebra/Matrix.h"
#include "io/IOField.h"
#include "calc.h"
#include "defs.h"
#include "helps.h"

using std::cout;
using std::endl;
using namespace OFELI;

namespace RITA {

int data::plot()
{
   int nb=0, pl=0, h=0, ret=0, k=0;
   _pl.lx = _pl.ly = _pl.component = 0;
   _pl.mx = _pl.Mx = _pl.my = _pl.My = 0.;
   _pl.title = _pl.mark = _pl.xaxist = _pl.yaxist = "";
   _pl.contour = 2;
   string f="", file="", tab="", ms="", vec="", his="", soft="gnuplot";
   const static vector<string> kw {"func$tion","file","tab$ulation","vect$or","hist$ory","mesh","iso$lines",
                                   "cont$our","x","y","log","lines","mark","title","label","use","compo$nent"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG("plot>")

   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case   0:
            NO_VALUE_ARG("plot>")
            f = _cmd->string_token(0);
            pl++;
            break;

         case   1:
            NO_VALUE_ARG("plot>")
            file = _cmd->string_token(0);
            pl++;
            break;

         case   2:
            NO_VALUE_ARG("plot>")
            tab = _cmd->string_token(0);
            pl++;
            break;

         case   3:
            NO_VALUE_ARG("plot>")
            vec = _cmd->string_token(0);
            pl++;
            break;

         case   4:
            NO_VALUE_ARG("plot>")
            his = _cmd->string_token(0);
            pl++;
            break;

         case   5:
            NO_VALUE_ARG("plot>")
            ms = _cmd->string_token(0);
            pl++;
            break;

         case   6:
            break;

         case   7:
            break;

         case   8:
            DDEF_PAR_R(0,"plot>",_pl.mx)
            DDEF_PAR_R(1,"plot>",_pl.Mx)
            break;

         case   9:
            DDEF_PAR_R(0,"plot>",_pl.my)
            DDEF_PAR_R(1,"plot>",_pl.My)
            break;

         case  10:
            DDEF_PAR_R(0,"plot>",_pl.lx)
            DDEF_PAR_R(1,"plot>",_pl.ly)
            break;

         case  11:
            break;

         case  12:
            _pl.mark = _cmd->string_token(0);
            break;

         case  13:
            _pl.title = _cmd->string_token(0);
            break;

         case  14:
            _pl.label = _cmd->string_token(0);
            break;

         case  15:
            soft = _cmd->string_token(0);
            break;

         case  16:
            DDEF_PAR_R(0,"plot>",_pl.component)
            break;

         case 100:
            cout << "Available arguments: " << Plot_help << endl;
            return 0;

         case 101:
            cout << Plot_Help << endl;
            return 0;

         default:
            UNKNOWN_ARG("plot>")
      }
   }
   if (nb_args>0) {
      if (pl==0)
         MSGR("plot>","No data have been defined to plot")
      else if (pl>1)
         MSGR("plot>","More than one data type has been defined")
      *_rita->ofh << "plot";
      if (_pl.mx!=0. && _pl.Mx!=0.)
         *_rita->ofh << " x=" << _pl.mx << "," << _pl.Mx;
      if (_pl.my!=0. && _pl.My!=0.)
         *_rita->ofh << " x=" << _pl.my << "," << _pl.My;
      if (_pl.lx!=0 && _pl.ly!=0)
         *_rita->ofh << " log=" << _pl.lx << "," << _pl.ly;
      if (_pl.mark!="")
         *_rita->ofh << " mark=" << _pl.mark;
      if (_pl.title!="")
         *_rita->ofh << " title=" << _pl.title;
      _pl.soft = soft;
      if (_pl.component>0)
         *_rita->ofh << " component=" << _pl.component;
      *_rita->ofh << endl;
      if (_pl.component==0)
         _pl.component = 1;
#ifdef USE_GMSH
      CHK_MSGR(soft!="gnuplot" && soft!="gmsh","plot>","Only gnuplot and gmsh are available for plotting")
#else
      CHK_MSGR(soft!="gnuplot","plot>","Only gnuplot is available for plotting")
#endif
      GnuplotPipe gp;
      if (f!="") {
         if (FctLabel.count(f)==0) {
            _rita->msg("plot>","Function "+f+" not found.");
            return 1;
         }
         _theFct = theFct[FctLabel[f]];
         if (_theFct->nb_var>2) {
            _rita->msg("plot>","No plotting available for functions with more than 2 variables.");
            return 1;
         }
         plot_fct();
         if (_pl.soft=="gnuplot") {
            for (auto const& p: _pl_cmd)
               gp.sendLine(p);
         }
         return 0;
      }
      if (file!="") {
         _rita->msg("plot>","This option has not been implemented yet.");
         return 1;
      }
      if (vec!="") {
         if (VectorLabel.count(vec)==0) {
            _rita->msg("plot>","Vector "+vec+" not found.");
            return 1;
         }
         _theVector = theVector[VectorLabel[vec]];
         plot_vect();
         if (_pl.soft=="gnuplot") {
            for (auto const& p: _pl_cmd)
               gp.sendLine(p);
         }
      }
      if (tab!="") {
         if (TabLabel.count(tab)==0) {
            _rita->msg("plot>","Tabulation "+tab+" not found.");
            return 1;
         }
         _theTab = theTab[TabLabel[tab]];
         *_rita->ofh << " tabulation=" << tab;
         plot_tab();
         if (_pl.soft=="gnuplot") {
            for (auto const& p: _pl_cmd)
               gp.sendLine(p);
         }
         return 0;
      }
      if (his!="") {
         if (HVectorLabel.count(his)==0) {
            _rita->msg("plot>","HVect "+his+" not found.");
            return 1;
         }
         _theHVector = theHVector[HVectorLabel[his]];
         *_rita->ofh << " history=" << his;
         plot_hist();
         if (_pl.soft=="gnuplot") {
            for (auto const& p: _pl_cmd)
               gp.sendLine(p);
         }
         return 0;
      }
      if (ms!="") {
         if (MeshLabel.count(ms)==0) {
            _rita->msg("plot>","Mesh "+ms+" not found.");
            return 1;
         }
         _theMesh = theMesh[MeshLabel[ms]];
         *_rita->ofh << " mesh=" << ms;
         plot_mesh();
         if (_pl.soft=="gnuplot") {
            for (auto const& p: _pl_cmd)
               gp.sendLine(p);
         }
         return 0;
      }
   }
   else {
      *_rita->ofh << "plot" << endl;
      for (;;) {
         if (_cmd->readline(sPrompt+"plot> ")<0)
            continue;
         int key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw);
         if (key>=200) {
            setDataExt(key);
            continue;
         }
         switch (key) {

            case   0:
               if (!_cmd->get(f)) {
                  *_rita->ofh << "  function " << f << endl;
                  pl++;
               }
               break;

            case   1:
               if (!_cmd->get(file)) {
                  *_rita->ofh << "  file " << file << endl;
                  pl++;
               }
               break;

            case   2:
               if (!_cmd->get(tab)) {
                  *_rita->ofh << "  tabulation " << tab << endl;
                  pl++;
               }
               break;

            case   3:
               if (!_cmd->get(vec)) {
                  *_rita->ofh << "  vector " << vec << endl;
                  pl++;
               }
               break;

            case   4:
               if (!_cmd->get(his)) {
                  *_rita->ofh << "  history " << his << endl;
                  pl++;
               }
               break;

            case   5:
               if (!_cmd->get(ms)) {
                  *_rita->ofh << "  mesh " << ms << endl;
                  pl++;
               }
               break;

            case   6:
               break;

            case   7:
               break;

            case   8:
               ret =  _cmd->get(_pl.mx);
               ret += _cmd->get(_pl.Mx);
               if (!ret)
                  *_rita->ofh << "x  " << _pl.mx << " " << _pl.Mx << endl;
               break;

            case   9:
               ret =  _cmd->get(_pl.my);
               ret += _cmd->get(_pl.My);
               if (!ret)
                  *_rita->ofh << "y  " << _pl.my << " " << _pl.My << endl;
               break;

            case  10:
               ret =  _cmd->get(_pl.lx);
               ret += _cmd->get(_pl.ly);
               if (!ret)
                  *_rita->ofh << "log " << _pl.lx << " " << _pl.ly << endl;
               break;

            case  11:
               ret = _cmd->get(_pl.mark);
               if (!ret)
                  *_rita->ofh << "mark " << _pl.mark << endl;
               break;

            case  12:
               ret = _cmd->get(_pl.mark);
               if (!ret)
                  *_rita->ofh << "mark " << _pl.mark << endl;
               break;

            case  13:
               ret = _cmd->get(_pl.title);
               if (!ret)
                  *_rita->ofh << "title  " << _pl.title << endl;
               break;

            case  14:
               ret = _cmd->get(_pl.label);
               if (!ret)
                  *_rita->ofh << "label  " << _pl.label << endl;
               break;

            case  15:
               ret = _cmd->get(soft);
#ifdef USE_GMSH
               CHK_MSGB(soft!="gnuplot" && soft!="gmsh","plot>","Only gnuplot and gmsh are available for plotting")
#else
               if (soft!="gnuplot") {
                  _rita->msg("plot>","Only gnuplot is available for plotting");
                  break;
               }
#endif
               if (!ret)
                  *_rita->ofh << "use  " << soft << endl;
               break;

            case  16:
               ret = _cmd->get(_pl.component);
               if (!ret)
                  *_rita->ofh << "component " << _pl.component << endl;
               break;

            case 100:
               h = 1;
               cout << "Available commands: " << Plot_help << endl;
               break;

            case 101:
               h = 1;
               cout << Plot_HHelp << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               _ret = _rita->_configure->run();
               break;

            case 104:
            case 105:
               if (h)
                  return 0;
               if (pl==0)
                  MSGR("plot>","No data have been defined to plot")
               else if (pl>1)
                  MSGR("plot>","More than one data type has been defined")
               {
                  if (_pl.component==0)
                     _pl.component = 1;
                  _pl.soft = soft;
                  GnuplotPipe gp;
                  if (f!="") {
                     k = FctLabel[f];
                     if (k<=0) {
                        _rita->msg("plot>","Function "+f+" not found.");
                        return 1;
                     }
                     _theFct = theFct[k];
                     plot_fct();
                  }
                  if (file!="") {

                  }
                  if (tab!="") {
                     k = TabLabel[tab];
                     if (k<=0) {
                        _rita->msg("plot>","Tabulation "+tab+" not found.");
                        return 1;
                     }
                     _theTab = theTab[k];
                     plot_tab();
                  }
                  if (vec!="") {
                     k = VectorLabel[vec];
                     if (k<=0) {
                        _rita->msg("plot>","Vector "+vec+" not found.");
                        return 1;
                     }
                     plot_vect();
                  }
                  if (his!="") {
                     k = HVectorLabel[vec];
                     if (k<=0) {
                        _rita->msg("plot>","History Vector "+his+" not found.");
                        return 1;
                     }
                     plot_hist();
                  }
                  if (ms!="") {
                     k = MeshLabel[ms];
                     if (k<=0) {
                        _rita->msg("plot>","Mesh "+ms+" not found.");
                        return 1;
                     }
                     _theMesh = theMesh[k];
                     plot_mesh();
                  }
                  if (_pl.soft=="gnuplot") {
                     for (auto const& p: _pl_cmd)
                        gp.sendLine(p);
                  }
               }
               return 0;

            case 106:
               _rita->setEcho(nb_args);
               break;

            default:
               DEFAULT_KW
         }
      }
   }
   return 0;
}


int data::plot_tab()
{
#ifdef USE_GMSH
   CHK_MSGR(_pl.soft=="gmsh","plot>","Tabulation plotting is not available with gmsh.");
#endif
   CHK_MSGR(_theTab->getNbVar(1)>1,"plot>","Plotting is available for one-variable cases only.")
   int np = _theTab->getSize(1,1);
   string file = "rita-gnuplot-tab.dat";
   temp_file.push_back(file);
   ofstream gm(file.c_str());
   double h = (_theTab->getMaxVar(1,1)-_theTab->getMinVar(1,1))/(np-1);
   double x = _theTab->getMinVar(1,1);
   gm << x << "  " << _theTab->Funct[0].Val(1) << endl;
   for (int i=1; i<np; ++i) {
      x += h;
      gm << x << "  " << _theTab->Funct[0].Val(i+1) << endl;
   }
   gm.close();
   _pl_cmd.push_back("plot '"+file+"' title \""+_pl.title+"\"");
   return 0;
}


int data::plot_mesh()
{
#ifdef USE_GMSH
   if (_pl.soft=="gmsh") {
      string file = "rita-temp.msh";
      _theMesh->save(file);
      temp_file.push_back(file);
      CHK_MSGR(system(("gmsh "+file).c_str()),"plot>","Unrecognizable system command.")
   }
#endif
   if (_pl.soft=="gnuplot") {
      CHK_MSGR(_theMesh->getDim()!=2,"plot>","Plotting is available for 2-D meshes only.")
      static map<int,size_t> sh = {{LINE,2},{QUADRILATERAL,4},{TRIANGLE,3}};
      string file = "rita-gnuplot-mesh.dat";
      temp_file.push_back(file);
      ofstream gm(file.c_str());
      for (auto const& e: _theMesh->theElements) {
         size_t m = sh[e->getShape()];
         if (m==0) {
            _rita->msg("plot>","Illegal element geometry.");
            return 1;
         }
         for (size_t n=1; n<=m; ++n)
            gm << (*e)(n)->getX() << "  " << (*e)(n)->getY() << endl;
         gm << (*e)(1)->getX() << "  " << (*e)(1)->getY() << endl << endl;
      }
      gm.close();
      _pl_cmd.push_back("plot '"+file+"' with lines");
   }
   return 0;
}


int data::plot_fct()
{
#ifdef USE_GMSH
   CHK_MSGR(_pl.soft=="gmsh","plot>","Function plotting is not available with gmsh.")
#endif
   string plt1 = _theFct->getExpression(), plt2="";
   if (_theFct->nb_var==1) {
      string v1 = _theFct->var[0];
      if (_pl.mx!=0. || _pl.Mx!=0.)
         _pl_cmd.push_back("set xrange ["+to_string(_pl.mx)+":"+to_string(_pl.Mx)+"]");
      if (_pl.my!=0. || _pl.My!=0.)
         _pl_cmd.push_back("set yrange ["+to_string(_pl.my)+":"+to_string(_pl.My)+"]");
      if (_pl.lx==1 && _pl.ly==0)
         _pl_cmd.push_back("set logscale x");
      else if (_pl.lx==0 && _pl.ly==1)
         _pl_cmd.push_back("set logscale y");
      else if (_pl.lx==1 && _pl.ly==1)
         _pl_cmd.push_back("set logscale xy");
      modify_exp(plt1,plt2);
      _pl_cmd.push_back("plot ["+v1+"=*:*] "+plt2+" title \""+_pl.title+"\"");
   }
   else if (_theFct->nb_var==2) {
      string v1 = _theFct->var[0], v2 = _theFct->var[1];
      _pl_cmd.push_back("set contour");
      _pl_cmd.push_back("set hidden3d");
//      _pl_cmd.push_back("set isosamples 10,10; set samples 10,10;");
      modify_exp(plt1,plt2);
      _pl_cmd.push_back("splot ["+v1+"=*:*] ["+v2+"=*:*]"+plt2);
//      _pl_cmd.push_back("set samples 4,10; replot");
//      _pl_cmd.push_back("set samples 10,4; replot");
   }
   else
      MSGR("plot>","No plotting available for functions with more than 2 variables.")
   return 0;
}


int data::plot_hist()
{
#ifdef USE_GMSH
   if (_pl.soft=="gmsh") {
      _rita->msg("plot>","History plotting is not available with gmsh.");
      return 1;
   }
#endif
   string file = "rita-gnuplot-hist.dat";
   temp_file.push_back(file);
   ofstream gm(file.c_str());
   for (int n=1; n<=_theHVector->nt; ++n) {
      Vect<double> *v = _theHVector->get(n);
      v = _theHVector->get(n);
         gm << _theHVector->getTime(n) << "  " << (*v)[_pl.component-1] << endl;
   }
   gm.close();
   if (_pl.mx!=0. || _pl.Mx!=0.)
      _pl_cmd.push_back("set xrange ["+to_string(_pl.mx)+":"+to_string(_pl.Mx)+"]");
   if (_pl.my!=0. || _pl.My!=0.)
      _pl_cmd.push_back("set yrange ["+to_string(_pl.my)+":"+to_string(_pl.My)+"]");
   if (_pl.lx==1 && _pl.ly==0)
      _pl_cmd.push_back("set logscale x");
   else if (_pl.lx==0 && _pl.ly==1)
      _pl_cmd.push_back("set logscale y");
   else if (_pl.lx==1 && _pl.ly==1)
      _pl_cmd.push_back("set logscale xy");
   _pl_cmd.push_back("plot 'rita-gnuplot-hist.dat' with lines");
   return 0;
}


int data::plot_vect()
{
   Vect<double> &v = *_theVector;
#ifdef USE_GMSH
   if (_pl.soft=="gmsh") {
      if (v.WithMesh()==false) {
         _rita->msg("plot>","Vector plotting with gmsh is available for mesh defined vectors only.");
         return 1;
      }
      string file = "rita-temp.pos";
      temp_file.push_back(file);
      saveGmsh(file,*_theVector);
      CHK_MSGR(system(("gmsh "+file).c_str()),"plot>","Unrecognizable system command.")
      return 0;
   }
#endif
   if (v.WithMesh()) {
      _theMesh = &_theVector->getMesh();
      string file = "rita-gnuplot.dat";
      temp_file.push_back(file);
      if (_theMesh->getDim()==1) {
         ofstream gm(file.c_str());
         for (size_t n=1; n<=_theMesh->getNbNodes(); ++n)
            gm << (*_theMesh)[n]->getX() << "  " << v(n) << endl;
         gm.close();
         if (_pl.title!="")
            _pl_cmd.push_back("set title "+_pl.title);
         if (_pl.mx!=0. || _pl.Mx!=0.)
            _pl_cmd.push_back("set xrange ["+to_string(_pl.mx)+":"+to_string(_pl.Mx)+"]");
         if (_pl.my!=0. || _pl.My!=0.)
            _pl_cmd.push_back("set yrange ["+to_string(_pl.my)+":"+to_string(_pl.My)+"]");
         if (_pl.lx==1 && _pl.ly==0)
            _pl_cmd.push_back("set logscale x");
         else if (_pl.lx==0 && _pl.ly==1)
            _pl_cmd.push_back("set logscale y");
         else if (_pl.lx==1 && _pl.ly==1)
            _pl_cmd.push_back("set logscale xy");
         _pl_cmd.push_back("plot '"+file+"' with lines");
      }
      else if (_theMesh->getDim()==2)
         MSGR("plot>","Plotting contours of mesh associated vectors is not available using gnuplot.")
/*
         ofstream gm(file.c_str());
         for (auto const& e: _theMesh->theElements) {
            for (size_t i=1; i<=e->getNbNodes(); ++i)
               gm << (*e)(i)->getX() << "  " << (*e)(i)->getY() << "  " << v((*e)(i)->n()) << endl;
            gm << endl;
         }
         gm.close();
         if (_pl.title!="")
            _pl_cmd.push_back("set title "+_pl.title);
         if (_pl.contour==1) {
            _pl_cmd.push_back("set cntrparam levels 20");
            _pl_cmd.push_back("set contour base");
            _pl_cmd.push_back("unset sur");
            _pl_cmd.push_back("set view map");
//            _pl_cmd.push_back("set xrange [0:4]");
//            _pl_cmd.push_back("set yrange [0:4]");
            _pl_cmd.push_back("set isosamples 20");
//            _pl_cmd.push_back("set samp 100");
//            _pl_cmd.push_back("set key rmargin");
            _pl_cmd.push_back("splot '"+file+"' with lines");
         }
         else if (_pl.contour==2) {
            _pl_cmd.push_back("set pm3d map");
//            _pl_cmd.push_back("set xrange [0:4]");
//            _pl_cmd.push_back("set yrange [0:4]");
//            _pl_cmd.push_back("set zrange [-1:1]");
//            _pl_cmd.push_back("set cbrange [-1:1]");
            _pl_cmd.push_back("set palette rgbformulae 22,13,10");
            _pl_cmd.push_back("set cntrparam cubicspline");
            _pl_cmd.push_back("splot '"+file+"'");
         }
      }*/
   }
   else if (v.WithGrid()) {
      string file = "rita-gnuplot.dat";
      temp_file.push_back(file);
      _theGrid = &_theVector->getGrid();
      ofstream gm(file.c_str());
      if (_theGrid->getDim()==1) {
         for (size_t n=1; n<=_theGrid->getNx()+1; ++n)
            gm << _theGrid->getX(n) << "  " << v(n) << endl;
         gm.close();
         if (_pl.title!="")
            _pl_cmd.push_back("set title "+_pl.title);
         if (_pl.mx!=0. || _pl.Mx!=0.)
            _pl_cmd.push_back("set xrange ["+to_string(_pl.mx)+":"+to_string(_pl.Mx)+"]");
         if (_pl.my!=0. || _pl.My!=0.)
            _pl_cmd.push_back("set yrange ["+to_string(_pl.my)+":"+to_string(_pl.My)+"]");
         if (_pl.lx==1 && _pl.ly==0)
            _pl_cmd.push_back("set logscale x");
         else if (_pl.lx==0 && _pl.ly==1)
            _pl_cmd.push_back("set logscale y");
         else if (_pl.lx==1 && _pl.ly==1)
            _pl_cmd.push_back("set logscale xy");
         _pl_cmd.push_back("plot '"+file+"' with lines");
      }
      else if (_theGrid->getDim()==2) {
         for (size_t i=1; i<=_theGrid->getNx()+1; ++i) {
            for (size_t j=1; j<=_theGrid->getNy()+1; ++j)
               gm << _theGrid->getX(i) << "  " << _theGrid->getY(j) << "  " << v(i,j) << endl;
            gm << endl;
         }
         gm.close();
         if (_pl.title!="")
            _pl_cmd.push_back("set title "+_pl.title);
         if (_pl.contour==1) {
            _pl_cmd.push_back("set cntrparam levels 10");
            _pl_cmd.push_back("set contour base");
            _pl_cmd.push_back("unset sur");
            _pl_cmd.push_back("set view map");
//            _pl_cmd.push_back("set xrange [0:4]");
//            _pl_cmd.push_back("set yrange [0:4]");
            _pl_cmd.push_back("set iso 100");
            _pl_cmd.push_back("set samp 100");
            _pl_cmd.push_back("set key rmargin");
         }
         else if (_pl.contour==2) {
            _pl_cmd.push_back("set pm3d map");
//            _pl_cmd.push_back("set xrange [0:4]");
//            _pl_cmd.push_back("set yrange [0:4]");
//            _pl_cmd.push_back("set zrange [-1:1]");
//            _pl_cmd.push_back("set cbrange [-1:1]");
            _pl_cmd.push_back("set palette rgbformulae 22,13,10");
         }
         _pl_cmd.push_back("splot '"+file+"'");
      }
      else
         MSGR("plot>","Plotting is available for vectors associated to meshes or grids only.")
   }
   return 0;
}


void data::modify_exp(const string& s1, string& s2)
{
   for (size_t i=0; i<s1.size(); ++i) {
      if (s1[i]=='^')
         s2.push_back('*'), s2.push_back('*');
      else
         s2.push_back(s1[i]);
   }
}

} /* namespace RITA */