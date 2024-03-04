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

                         Definition of class 'configure'

  ==============================================================================*/

#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
using std::string;
using std::ifstream;
using std::ofstream;

#include "cmd.h"

namespace RITA {

#define RITA_PROMPT "rita>"
class rita;

class configure
{

 public:

    configure(rita *r, cmd *command);
    ~configure();
    int run();
    void setVerbose(int verb) { _verb = verb; }
    int getVerbose() const { return _verb; }
    std::ofstream* getOStreamLog() { return &_ofl; }
    std::ofstream* getOStreamHistory() { return &_ofh; }
    int getSaveResults() const { return _save_results; }
    void set(cmd* command) { _cmd = command; }
    void set(string cf);
    int read();
    void save();
    void init();
    
 private:

    const string currentDateTime() {
       time_t     now = time(0);
       struct tm  tstruct;
       char       buf[80];
       tstruct = *localtime(&now);
       strftime(buf,sizeof(buf),"%Y-%m-%d  %X",&tstruct);
       return buf;
    }

    rita *_rita;
    int _verb, _key, _save_results;
    string _HOME, _his_file, _log_file;
    ofstream _ofh, _ofl, _ocf;
    ifstream _icf;
    const vector<string> _kw {"verb$osity","save$-results","history$-file","log$-file","end"};
    cmd *_cmd;
};

} /* namespace RITA */