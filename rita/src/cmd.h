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

                        Definition of class Cmd

  ==============================================================================*/

#pragma once

#define USE_CTRL_D

#include <stdlib.h>
#include <ctype.h>
#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::endl;
using std::ostream;

#include <sstream>
#include <fstream>
using std::istringstream;
using std::ifstream;

#ifdef USE_CTRL_D
#include <termios.h> 
#include <csignal> 
#include <cstdlib>
#endif

#include <string>
using std::string;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace RITA {

/*! \file cmd.h
 *  \brief Definition file for class cmd.
 */

/*! \class cmd
 *  \brief 
 *
 * \author Rachid Touzani
 * \copyright GNU Public License
 */

class rita;

class cmd {

 private:
   rita *_rita;
   ifstream *_is;
   string _prompt, _buffer, _arg, _tok, _current_tok;
   vector<string> _toks;
   vector<string> _cmd_hist;
   bool _command, _comment, _with_kw;
   int _ch_script;
   vector<string> _word;
   size_t _ind;
   size_t _nb_args, _script_line_nb;
   int split();
   int get();
#ifdef USE_CTRL_D
   struct termios _old_termios, _new_termios;
   static void handler(int sig);
#endif

 public:

   cmd(rita *r);
   cmd(ifstream& is, rita *r);
   ~cmd();
   bool isInputFile() const { return (_is!=nullptr); }
   void setIFStream(ifstream *is) { _is = is; _script_line_nb = 0; }
   ifstream *getIFStream() const { return _is; }
   int readline(string p="");
   void setPrompt(string p) { _prompt = p; }
   void setErrorMsg(string e);
   int setNbArg(size_t n, const string& s="", int opt=0);
   int setNbArg(size_t n, size_t m, const string& s="", int opt=0);
   void setNbArg() { _nb_args = -1; }
   size_t getNbArgs() const { return _nb_args; }
   int getChScript() const { return _ch_script; }
   int get(string& s);
   int get(int &i);
   int get(double& d);
   int get(const vector<string>& kw, string& s);
   int getKW(const vector<string> &kw);
   int getKW(const vector<string> &kw1, const vector<string> &kw2);
   int getKW(const vector<string> &kw1, const vector<string> &kw2, const vector<string>& kw3);
   int getScriptLineNb() const { return _script_line_nb; }
   string getToken() const { return _arg; }
   void flush() { _arg = ""; _buffer.clear(); }
   size_t checkEq() const { return _buffer.find('='); }
   string token() const { return _word[_ind-1]; }
   string buffer() const { return _buffer; }
   string string_token() const { return _tok; }
   double double_token() const { return stod(_tok); }
   int int_token() const;
   string string_token(int i) const { return _toks[i]; }
   double double_token(int i) const { return stod(_toks[i]); }
   int int_token(int i) const { return stoi(_toks[i]); }
   int getArg(string delimiter="=");
   int getLRArgs(string& arg1, string& arg2, string delimiter="=");
   int getArgs(int &nb, string del1="=", string del2=",");
   string Arg() const { return _arg; }
   void set(const vector<string>& arg);
   void set(const vector<string>& arg1, const vector<string>& arg2);
   void set(const vector<string>& arg1, const vector<string>& arg2, const vector<string>& arg3);
   const vector<string> *_kw, *_gkw, *_dkw;
   string& trim(string& s, const string& chars = "\t\n\v\f\r ");
   int find_kw(const string &arg);
   bool isValidNumber(const string& s);
   bool isNumeric(const string& s);
   bool isNumeric(const string& s, double& d);
   bool isNumeric(const string& s, int& d);
   string getVar();

   friend class rita;
};

} /* namespace RITA */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
