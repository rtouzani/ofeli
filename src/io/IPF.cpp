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

            Implementation of class 'IPF' to read parameter input files

  ==============================================================================*/

#include "io/XMLParser.h"
#include "util/banner.h"
#include "io/IPF.h"

namespace OFELI {

IPF::IPF()
{
}


IPF::IPF(const string& prog,
         const string& file)
{
   _file = file;
   ifstream inf(file.c_str());
   string cc;
   inf >> cc;
   inf.close();
   init();
   XMLParser p(file,XMLParser::PROJECT);
   p.get(*this);
   if (_restart_file.size()==0)
      _restart_file = _project + ".res";
   if (_save_file.size()==0)
      _save_file = _project + ".sav";
   for (size_t i=0; i<MAX_NB_PAR; i++) {
      if (_plot_file[i].size()==0) {
         _plot_file[i] = _project;
         _plot_file[i] += itos(i+1) + ".pl";
      }
   }
}


IPF::IPF(const string& file)
{
   _file = file;
   ifstream inf(file.c_str());
   string cc;
   inf >> cc;
   inf.close();

   init();
   XMLParser p(file,XMLParser::PROJECT);
   p.get(*this);
   if (_restart_file.size()==0)
      _restart_file = _project + ".res";
   if (_save_file.size()==0)
      _save_file = _project + ".sav";
   for (size_t i=0; i<MAX_NB_PAR; i++) {
      if (_plot_file[i].size()==0) {
         _plot_file[i] = _project + "-";
         _plot_file[i] += itos(i+1) + ".pl";
      }
   }
}


IPF::~IPF() { }


void IPF::init()
{
   _verbose = 1; 
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
   for (size_t i=0; i<MAX_NB_PAR; i++) {
      _real_par[i] = 0;
      _int_par[i] = 0;
      _complex_par[i] = 0;
   }
   _init_file = _file;
   _restart_file = _file;
   _domain_file = _file;
   for (size_t i=0; i<MAX_NB_PAR; i++) {
      _mesh_file[i] = _file;
      _data_file[i] = _file;
   }
   for (size_t j=0; j<MAX_NB_MATERIALS; j++)
      _nmat[j] = 0;
}


int IPF::contains(const string& label) const
{
   for (size_t i=0; i<_param_label.size(); i++) {
      if (_param_label[i]==label)
         return int(i+1);
   }
   return 0;
}


string IPF::getString(const string& label) const
{
   int i=contains(label);
   try {
      if (i>0)
         return _param_value[i-1];
      THROW_RT("getString(string): Parameter " + label + " unfound in project file.");
   }
   CATCH("IPF");
   return " ";
}


string IPF::getString(const string& label,
                            string  def) const
{
   int i=contains(label);
   if (i>0)
      return _param_value[i-1];
   return def;
}
   
    
int IPF::getInteger(const string& label) const
{
   int i=contains(label);
   try {
      if (i>0)
         return atoi(_param_value[i-1].c_str());
      THROW_RT("getInteger(string): Parameter " + label + " unfound in project file.");
   }
   CATCH("IPF");
   return 0;
}
   
   
int IPF::getInteger(const string& label,
                          int     def) const
{
   int i=contains(label);
   if (i>0)
      return atoi(_param_value[i-1].c_str());
   return def;
}


real_t IPF::getDouble(const string& label) const
{
   int i=contains(label);
   try {
      if (i>0)
         return atof(_param_value[i-1].c_str());
      THROW_RT("getDouble(string): Parameter " + label + " unfound in project file.");
    }
    CATCH("IPF");
    return 0;
}
   
   
real_t IPF::getDouble(const string& label,
                            real_t  def) const
{
   int i=contains(label);
   if (i>0)
      return atof(_param_value[i-1].c_str());
   return def;
}


complex_t IPF::getComplex(const string& label) const
{
   real_t ar, ai;
   complex_t a;
   int i=contains(label);
   try {
      if (i>0) {
         ar = atof(_param_value[i-1].c_str());
         ai = atof(_param_value[i-1].c_str());
         a = complex_t(ar,ai);
         return a;
      }
      THROW_RT("getComplex(string): Parameter " + label + " unfound in project file.");
   }
   CATCH("IPF");
   return 0;
}
   
   
complex_t IPF::getComplex(const string&   label,
                                complex_t def) const
{
   real_t ar, ai;
   complex_t a;
   int i=contains(label);
   if (i>0) {
      ar = atof(_param_value[i-1].c_str());
      ai = atof(_param_value[i-1].c_str());
      a = complex_t(ar,ai);
      return a;
   }
   return def;
}


void IPF::get(const string&       label,
                    Vect<real_t>& a) const
{
   for (size_t i=0; i<_array_label.size(); i++) {
      if (_array_label[i]==label) {
         a.setSize(_array_size[i]);
         for (size_t j=0; j<_array_size[i]; j++)
            a[j] = (_array_value[j])[i];
         return;
      }
   }
   try {
      THROW_RT("get(string,Vect<real_t>): Parameter " + label + " unfound in project file.");
   }
   CATCH("IPF");
}


void IPF::get(const string& label,
                    int&    a) const
{
   a = getInteger(label);
}


void IPF::get(const string& label,
                    real_t &a) const
{
   a = getDouble(label);
}


void IPF::get(const string&    label,
                    complex_t& a) const
{
   a = getComplex(label);
}


void IPF::get(const string& label,
                    string& a) const
{
   a = getString(label);
}

} /* namespace OFELI */
