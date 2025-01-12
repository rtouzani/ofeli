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

                            A class to manage dialogs

  ==============================================================================*/

#include "cmd.h"
#include "rita.h"
#include <algorithm>
#include <termios.h>
#include "defs.h"

namespace RITA {

cmd::cmd(rita *r)
    : _rita(r), _is(nullptr), _ch_script(0), _gkw(nullptr)
{
#ifdef USE_CTRL_D
   setvbuf(stdout,nullptr,_IONBF,0);
   tcgetattr(0,&_old_termios);
   signal(SIGINT,handler);
   _new_termios = _old_termios;
   _new_termios.c_cc[VINTR] = 4;
   tcsetattr(0,TCSANOW,&_new_termios);
#endif
}


cmd::cmd(ifstream& is, rita *r)
    : _rita(r), _is(&is), _ch_script(0), _gkw(nullptr)
{
}


cmd::~cmd()
{
}


#ifdef USE_CTRL_D
void cmd::handler(int sig)
{
   (void)sig;
   cout << "Do you really want to quit (y/n) ? ";
   string ans;
   std::cin >> ans;
   if (ans=="y" || ans=="Y")
      exit(0);
}
#endif


int cmd::readline(string p)
{
   _comment = _command = false;
   if (_is!=nullptr) {
      _ch_script = 1;
      getline(*_is,_buffer);
      _cmd_hist.push_back(_buffer);
      _script_line_nb++;
      if (_is->eof()) {
         _is = nullptr;
         _ch_script = -1;
         return -5;
      }
   }
   else {
      _ch_script = 0;
      cout << p;
      get();
   }
   trim(_buffer);
   if (_buffer.length()==0)
      readline(p);
   if (_rita->getEcho())
      cout << sPrompt+"> " << _buffer << endl;

// Execute a system command
   if (_buffer[0]=='!') {
      _command = true;
      if (_rita->getEcho())
         cout << _buffer << endl;
      if (system(_buffer.substr(1,_buffer.size()).c_str()))
         cout << "Error in system command." << endl;
      return -3;
   }

// Comment line
   else if (_buffer[0]=='#') {
      _comment = true;
      return -4;
   }

// Remove a space from line starting with '='
   else if (_buffer[0]=='=' && _buffer[1]!=' ')
      _buffer.insert(1," ");

// End
   if (_buffer=="exit" || _buffer=="quit") {
      if (_rita->getEcho())
         cout << sPrompt+"> " << _buffer << endl;
      _rita->finish();
   }

// Tokenize
   bool q = false;
   if (_buffer.find('\"')<_buffer.size())
      q = true;
   istringstream is(_buffer);
   _comment = _command = false;
   _word.clear();
   string token;
   while (!is.fail()) {
      is >> token;
      _word.push_back(token);
   }
   _ind = 0;
   _word.pop_back();
   _nb_args = _word.size() - 1;
   if (q)
      return split();
   return 0;
}


string cmd::getVar()
{
   int i=0;
   string s = "";
   while (_word[0][i]!=' ' && _word[0][i]!='=')
      s.push_back(_word[0][i++]);
   return s;
}


int cmd::get()
{
   _buffer = "";
   static struct termios t1, t2;
   tcgetattr(0,&t1);
   t2 = t1;
   t2.c_lflag &= ~(ICANON|ECHO);
   tcsetattr(0,TCSANOW,&t2);

   int k=_cmd_hist.size()-1;
   while (1) {
      char c = getchar();

      if (c==27) {
         if (getchar()==91 && getchar()==65 && k>=0) {
            for (size_t i=0; i<_buffer.length(); ++i)
               cout << "\b \b";
            _buffer = _cmd_hist[k--];
            cout << _buffer;
         }
         continue;
      }
      else if (c!=0x7f)
         _buffer += c;

//    if backspace
      if (c==0x7f && _buffer.length()!=0) {
         cout << "\b \b";
         _buffer.pop_back();
         continue;
      }
      if (_buffer.length()==0)
         continue;

//    End of line
      if (c=='\n') {
         _buffer.pop_back();
         break;
      }
      cout << c;
   }
   cout << endl;
   tcsetattr(0,TCSANOW,&t1);
   _cmd_hist.push_back(_buffer);
   return 0;
}


int cmd::split()
{
   int k=0;
   vector<string> w;
   bool q = false;
   w.push_back("");
   for (size_t i=0; i<_word.size(); ++i) {
      for (size_t j=0; j<_word[i].length(); ++j) {
         char z = _word[i][j];
         if (z=='\"')
            q = !q;
         else
            w[k] += z;
      }
      if (!q)
         k++, w.push_back("");
      else if (i<_word.size()-1)
         w[k] += " ";
   }
   w.pop_back();
   _word = w;
   _nb_args = _word.size() - 1;
   return 0;
}


int cmd::setNbArg(size_t        n,
                  size_t        m,
                  const string& s,
                  int           opt)
{
   if (n<_nb_args && !opt)
      IGNORED_ARGS
   else if (m>_nb_args)
      MISSING_ARGS(s)
   return 0;
}


int cmd::setNbArg(size_t        n,
                  const string& s,
                  int           opt)
{
   if (n<_nb_args && !opt)
      IGNORED_ARGS
   else if (n>_nb_args)
      MISSING_ARGS(s)
   return 0;
}


int cmd::get(const vector<string>& kw,
             string&               s)
{
   _kw = &kw;
   get(s);
   auto it = find(_kw->begin(),_kw->end(),s);
   if (it != _kw->end())
      return distance(_kw->begin(),it);
   else
      return -1;
}


int cmd::getKW(const vector<string> &kw1,
               const vector<string> &kw2)
{
   _kw = &kw1;
   _gkw = &kw2;
   if (++_ind>_word.size())
      return -1;
   return find_kw(_word[_ind-1]);
}


int cmd::getKW(const vector<string> &kw1,
               const vector<string> &kw2,
               const vector<string> &kw3)
{
   _kw = &kw1;
   _gkw = &kw2;
   _dkw = &kw3;
   if (++_ind>_word.size())
      return -1;
   return find_kw(_word[_ind-1]);
}


int cmd::getKW(const vector<string> &kw)
{
   _kw = &kw;
   if (++_ind>_word.size())
      return -1;
   return find_kw(_word[_ind-1]);
}


void cmd::set(const vector<string>& arg)
{
   _kw = &arg;
   _with_kw = true;
}


void cmd::set(const vector<string>& arg1,
              const vector<string>& arg2)
{
   _kw = &arg1;
   _gkw = &arg2;
   _with_kw = true;
}


void cmd::set(const vector<string>& arg1,
              const vector<string>& arg2,
              const vector<string>& arg3)
{
   _kw = &arg1;
   _gkw = &arg2;
   _dkw = &arg3;
   _with_kw = true;
}


int cmd::find_kw(const string &s)
{
   for (auto it=_kw->begin(); it!=_kw->end(); it++) {
      int n = int(it->find(SHORT_HAND));
      if (it->substr(0,n)==s.substr(0,n))
         return distance(_kw->begin(),it);
   }
   if (_gkw==nullptr)
      return -1;
   for (auto it=_gkw->begin(); it!=_gkw->end(); it++) {
      int n = int(it->find(SHORT_HAND));
      if (it->substr(0,n)==s.substr(0,n))
         return distance(_gkw->begin(),it)+100;
   }
   if (_dkw==nullptr)
      return -1;
   for (auto it=_dkw->begin(); it!=_dkw->end(); it++) {
      int n = int(it->find(SHORT_HAND));
      if (it->substr(0,n)==s.substr(0,n))
         return distance(_dkw->begin(),it)+200;
   }
   _arg = s;
   return -1;
}


int cmd::getArg(string delimiter)
{
   size_t pos=0, i=0;
   int ret=0;
   get(_tok);
   if ((pos=_tok.find(delimiter)) == string::npos) {
      ret = find_kw(_tok);
      _arg = _tok;
      _tok = "";
      return ret;
   }
   while ((pos=_tok.find(delimiter)) != string::npos) {
      _arg = _tok.substr(0,pos);
      _tok.erase(0,pos+delimiter.length());
      i++;
   }
   if (i>1)
      return 2;
   else if (i==0) {
      _tok = "";
      return 0;
   }
   return find_kw(_arg);
}


int cmd::int_token() const
{
   return stoi(_tok);
}


int cmd::getLRArgs(string& arg1,
                   string& arg2,
                   string  delimiter)
{
   size_t pos=0, i=0;
   int ret=0;
   get(_tok);
   if ((pos=_tok.find(delimiter)) == string::npos) {
      arg1 = _tok;
      return ret;
   }
   while ((pos=_tok.find(delimiter)) != string::npos) {
      arg2 = _tok.substr(0,pos);
      _tok.erase(0,pos+delimiter.length());
      i++;
   }
   if (i>1)
      return 2;
   else if (i==0) {
      _tok = "";
      return -1;
   }
   arg2 = _arg;
   return 0;
}


int cmd::getArgs(int&   nb,
                 string del1,
                 string del2)
{
   size_t pos=0, i=0;
   string tok;
   _toks.clear();
   get(tok);
   while (tok.find(del1) != string::npos) {
      pos = tok.find(del1);
      _arg = tok.substr(0,pos);
      tok.erase(0,pos+del1.length());
      i++;
   }
   if (pos==0) {
      nb = 0;
      return find_kw(tok);    
   }
   if (i>1)
      return 2;
   else if (i==0) {
      _tok = "";
      nb = 0;
      return 0;
   }
   while (tok.find(del2) != string::npos) {
      pos = tok.find(del2);
      _toks.push_back(tok.substr(0,pos));
      tok.erase(0,pos+1);
   }
   _toks.push_back(tok);
   nb = _toks.size();
   return find_kw(_arg);
}


int cmd::get(string& s)
{
   cout << _prompt;
   if (++_ind>_word.size()) 
      return -1;
   if (_word[_ind-1][0]=='\"') {
      s = "";
      size_t nb = 0;
      while (_word[_ind-1][_word[_ind-1].size()-1]!='\"') {
         nb++;
         s += _word[_ind-1] + " ";
         if (_ind==nb && _word[_ind][_word[_ind].size()-1]!='\"')
            return -1;
         _ind++;
      }
      s += _word[_ind-1];
      s.erase(0,1);
      s.erase(s.size()-1,1);
   }
   else
      s = _word[_ind-1];
   return 0;
}


int cmd::get(int& i)
{
   if (++_ind>_word.size()) 
      return -1;
   cout << _prompt;
   if (!isValidNumber(_word[_ind-1])) {
      cout << "Error in input: Expecting an integer. " << endl;
      return 1;
   }
   i = std::stoi(_word[_ind-1]);
   return 0;
}


int cmd::get(double& d)
{
   if (++_ind>_word.size())
      return -1;
   cout << _prompt;
   if (!isNumeric(_word[_ind-1])) {
      cout << "Error in input: Expecting a number. " << endl;
      return 1;
   }
   d = std::stof(_word[_ind-1]);
   return 0;
}


bool cmd::isValidNumber(const string& s)
{
   size_t i=0;
   if (s[0]=='-' || s[0]=='+')
      i++;
   for (; i<s.length(); ++i)
      if (s[i]<48 || s[i]>57)
         return false;
   return true;
}


bool cmd::isNumeric(const string& s)
{
   istringstream iss(s.c_str());
   double d;
   iss >> d;
   if (!iss)
      return false;
   return (iss.rdbuf()->in_avail()==0);
}


bool cmd::isNumeric(const string& s,
                    double&       d)
{
   istringstream iss(s.c_str());
   iss >> d;
   if (!iss)
      return false;
   return (iss.rdbuf()->in_avail()==0);
}


bool cmd::isNumeric(const string& s,
                    int&          d)
{
   istringstream iss(s.c_str());
   iss >> d;
   if (!iss)
      return false;
   return (iss.rdbuf()->in_avail() == 0);
}


string& cmd::trim(string&       s,
                  const string& chars)
{
   s.erase(0,s.find_first_not_of(chars));
   s.erase(s.find_last_not_of(chars)+1);
   return s;
}

} /* namespace RITA */
