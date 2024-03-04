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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

  ==============================================================================

                     Implementation of class 'configure'

  ==============================================================================*/

#include "configure.h"
#include "rita.h"
#include "helps.h"

namespace RITA {

configure::configure(rita *r, cmd *command)
          : _rita(r), _verb(1), _save_results(1), _his_file(".rita.his"),
            _log_file(".rita.log"), _cmd(command)
{
   init();
}


configure::~configure()
{
   _ocf.open((_HOME+"/.rita").c_str());
   save();
}


void configure::save()
{
   _ocf << "# rita configuration file" << endl;
   _ocf << "# " << currentDateTime() << "\n#\n";
   _ocf << "verbosity " << _verb << endl;
   _ocf << "save-results " << _save_results << endl;
   _ocf << "history-file " << _his_file << endl;
   _ocf << "log-file " << _log_file << endl;
   _ocf << "end" << endl;
   _ocf.close();
}


void configure::init()
{
   _HOME = getenv("HOME");
   _icf.open((_HOME+"/.rita").c_str());
   if (_icf.fail()) {
      _icf.close();
      _ocf.open((_HOME+"/.rita").c_str());
   }
   else {
      read();
      _ocf.open((_HOME+"/.rita.backup").c_str());
   }
   save();
   _ofl.open(_log_file);
   _ofl << "# rita log file" << endl;
   _ofl << "# " << currentDateTime() << "\n#\n";
   _ofh.open(_his_file);
   _ofh << "# rita history file" << endl;
   _ofh << "# " << currentDateTime() << "\n#\n";
}


int configure::read()
{
   cmd com(_icf,_rita);
   while (1) {
      if (com.readline()<0)
         continue;
      int key = com.getKW(_kw);
      switch (key) {

         case 0:
            com.get(_verb);
            break;

         case 1:
            com.get(_save_results);
            break;

         case 2:
            com.get(_his_file);
            break;

         case 3:
            com.get(_log_file);
            break;

         case 4:
            _icf.close();
            return 0;

         default:
            _rita->msg("set>:","Unknown setting: "+com.token(),
                       "Available settings: verbosity, save-results, history, log, end");
            return 1;
      }
   }
   return 0;
}


int configure::run()
{
   bool verb_ok=false, hist_ok=false, log_ok=false, save_ok=false;
   string hfile, lfile, buffer;
   ifstream is;
   _cmd->set(_kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(sPrompt)
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArg()) {

         case 0:
            _verb = _cmd->int_token();
            verb_ok = true;
            break;

         case 1:
            hist_ok = true;
            hfile = _his_file;
            _his_file = _cmd->string_token();
            break;

         case 2:
            log_ok = true;
            lfile = _log_file;
            _log_file = _cmd->string_token();
            break;

         case 100:
            cout << "Available arguments: " << Config_help << endl;
            return 0;
         
         case 101:
            cout << Config_Help << endl;
            return 0;

         case 106:
            _rita->setEcho(nb_args);
            break;

         case 108:
         case 109:
            return 500;
            break;

         case 110:
            return 600;
            break;

         default:
            UNKNOWN_ARG(sPrompt+"set>")
            return 1;
       }
   }
   if (nb_args>0) {
      _ofh << "set";
      if (verb_ok) {
         if (_verb<0 || _verb>10) {
            _rita->msg("set>","Illegal value of verbosity: "+to_string(_verb));
            return 1;
         }
         _ofh << " verbosity=" << _verb;
      }
      if (save_ok) {
         if (_save_results<0) {
            _rita->msg("set>","Illegal value of save: "+to_string(_save_results));
            return 1;
         }
         _ofh << " save-results=" << _save_results;
      }
      if (hist_ok) {
         _ofh.close();
         is.open(hfile);
         _ofh.open(_his_file.c_str());
         while (!is.eof()) {
            getline(is,buffer);
            _ofh << buffer << endl;
         }
      }
      if (log_ok) {
         _ofl.close();
         is.open(lfile);
         _ofl.open(_log_file.c_str());
         while (!is.eof()) {
            getline(is,buffer);
            _ofl << buffer << endl;
         }
      }
   }
   return 0;
}

} /* namespace RITA */