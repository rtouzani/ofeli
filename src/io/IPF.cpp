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

            Implementation of class 'IPF' to read parameter input files

  ==============================================================================*/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "io/XMLParser.h"
#include "util/banner.h"
#include "io/IPF.h"
#include "linear_algebra/Vect_impl.h"
#include "OFELIException.h"

namespace OFELI {

int Verbosity = 0;

IPF::IPF()
{
}


IPF::IPF(const string& prog,
         const string& file)
{
   _file = file;
   std::ifstream inf(file.c_str());
   string cc;
   inf >> cc;
   inf.close();
   init();
   XMLParser p(file,EType::PROJECT);
   p.get(*this);
   if (_restart_file.size()==0)
      _restart_file = _project + ".res";
   if (_save_file.size()==0)
      _save_file = _project + ".sav";
//   if (_plot_file.size()==0)
//      _plot_file.push_back(_project+"-1.pl");
   Verbosity = getVerbose();
}


IPF::IPF(const string& file)
{
   _file = file;
   std::ifstream inf(file.c_str());
   string cc;
   inf >> cc;
   inf.close();
   init();
   XMLParser p(file,EType::PROJECT);
   p.get(*this);
   if (_restart_file.size()==0)
      _restart_file = _project + ".res";
   if (_save_file.size()==0)
      _save_file = _project + ".sav";
//   if (_plot_file.size()==0)
//      _plot_file.push_back(_project+"-1.pl");
   Verbosity = getVerbose();
}


IPF::~IPF() { }


void IPF::init()
{
   _verbose = 0;
   _save = _plot = 0;
   _output = 1;
   _nb_iter = 100;
   _tolerance = 1.e-6;
   _nb_steps = 10;
   _time_step = 0.1; 
   _max_time = 1;
   _nb_mat = _nb_int_par = _nb_real_par = _nb_complex_par = 0;
   _nb_aux_files = _nb_data_files = _nb_mesh_files = _nb_plot_files = _nb_point_double_par = 0;
   _bc = _bf = _sf = _ini = _data = 0;
   _init_file = _file;
   _restart_file = _file;
   _domain_file = _file;
   if (_mesh_file.size()==0)
      _mesh_file.push_back(_file);
   if (_data_file.size()==0)
      _data_file.push_back(_file);
}


string IPF::getString(const string& label)
{
   if (_param.find(label)!=_param.end())
      return _param[label];
   throw OFELIException("In IPF::getString(string): Parameter " + label + " unfound in project file.");
   return "";
}


string IPF::getString(const string& label,
                      string        def)
{
   if (_param.find(label)!=_param.end())
      return _param[label];
   return def;
}


int IPF::getInteger(const string& label)
{
   if (_param.find(label)!=_param.end())
      return stoi(_param[label]);
   throw OFELIException("In IPF::getInteger(string): Parameter " + label + " unfound in project file.");
   return 0;
}


int IPF::getInteger(const string& label,
                    int           def)
{
   if (_param.find(label)!=_param.end())
      return stoi(_param[label]);
   return def;
}


real_t IPF::getDouble(const string& label)
{
   if (_param.find(label)!=_param.end())
      return stod(_param[label]);
   throw OFELIException("In IPF::getDouble(string): Parameter " + label + " unfound in project file.");
   return 0;
}


real_t IPF::getDouble(const string& label,
                      real_t        def)
{
   if (_param.find(label)!=_param.end())
      return stod(_param[label]);
   return def;
}


complex_t IPF::getComplex(const string& label)
{
   if (_cparam.find(label)!=_cparam.end())
      return complex_t(stod(_cparam[label].first),stod(_cparam[label].second));
   throw OFELIException("In IPF::getComplex(string): Parameter " + label + " unfound in project file.");
   return 0;
}
   
   
complex_t IPF::getComplex(const string& label,
                          complex_t     def)
{
   if (_cparam.find(label)!=_cparam.end())
      return complex_t(stod(_cparam[label].first),stod(_cparam[label].second));
   return def;
}


void IPF::get(const string& label,
              Vect<real_t>& a)
{
   if (_vparam.find(label)!=_vparam.end()) {
      a.setSize(_vparam[label].size());
      for (size_t i=0; i<a.size(); ++i)
         a[i] = stod(_vparam[label][i]);
      return;
   }
   throw OFELIException("In IPF::get(string,Vect<real_t>): Parameter " + label + " unfound in project file.");
}


void IPF::get(const string& label,
              int&          a)
{
   a = getInteger(label);
}


void IPF::get(const string& label,
              real_t&       a)
{
   a = getDouble(label);
}


void IPF::get(const string& label,
              complex_t&    a)
{
   a = getComplex(label);
}


void IPF::get(const string& label,
              string&       a)
{
   a = getString(label);
}


void IPF::set_data_file(const string& s)
{
   if (_data_file.size()==1)
      _data_file[0] = s;
   else
      _data_file.push_back(s);
}


void IPF::set_mesh_file(const string& s)
{
   if (_mesh_file.size()==1)
      _mesh_file[0] = s;
   else
      _mesh_file.push_back(s);
}

} /* namespace OFELI */
