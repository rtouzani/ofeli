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

                           Definition of some macros

  ==============================================================================*/

#pragma once

#define SHORT_HAND "$"

#define CATCH catch(ritaException &e) {                                    \
                 std::cout << "RITA error: " << e.what() << endl;          \
                 return 1;                                                 \
              }                                                            \
              catch(runtime_error &e) {                                    \
                 std::cout << "RITA Runtime error: " << e.what() << endl;  \
                 return 1;                                                 \
              }                                                            \
              catch( ... ) {                                               \
                 std::cout << "RITA Unexpected error: " << endl;           \
                 return 1;                                                 \
              }

#define RET(i) { _pr = _PR; return i; }
#define MSG(p,m)       { _rita->msg(p,m); }
#define MSGR(p,m)      { _rita->msg(p,m); return 1; }
#define MSG1R(p,m,t,n) { _rita->msg(p,m,t,n); return 1; }
#define MSGB(p,m)      { ret = 1; _rita->msg(p,m); break; }
#define MSG1B(p,m,t,n) { ret = 1; _rita->msg(p,m,t,n); break; }
#define CHK_MSG(c,p,m) { if (c) { _rita->msg(p,m); } }
#define CHK_MSGR0(c,p,m)     { if (c) { _rita->msg(p,m); return 0; } }
#define CHK_MSGR(c,p,m)      { if (c) { _rita->msg(p,m); return 1; } }
#define CHK_MSG1R(c,p,m,t,n) { if (c) { _rita->msg(p,m,t,n); return 1; } }
#define CHK_MSGB(c,p,m)      { if (c) { ret = 1; _rita->msg(p,m); break; } }
#define CHK_MSG1B(c,p,m,t,n) { if (c) { ret = 1; _rita->msg(p,m,t,n); break; } }
#define MSG_PDE(m,s)         { _rita->msg(_pr,m); log.s = true; break; }
#define CHK_MSG_PDE(c,m,s)   { if (c) { _rita->msg(_pr,m); log.s = true; } break; }
#define CHK_MSGR_PDE(c,p,m)  { if (c) { _rita->msg(p,m); return 1; } }
#define IGNORED_ARGS { cout << "Command argument(s) ignored." << endl; return 0; }
#define MISSING_ARGS(s) { cout << "Missing argument(s): " << s << endl; return 1; }
#define MISSING_FUNCT_NAME(s) { _rita->msg(s,"Missing function name.","",1); break; }
#define REDEFINED(s) { _rita->msg("",s+name+" redefined"); }
#define ALREADY_USED(c,d) if (c.count(s) && dn[s].active) { \
                             _rita->msg("","Name "+s+" already used for "+d); return 1; }

#define READLINE if (_cmd->readline(_pr+" ")<0) continue; \
                 key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw); \
                 if (key>=200) { _data->setDataExt(key); continue; }

#define CHK_PB(s) _data->dn[s].dt!=DataType::LS && _data->dn[s].dt!=DataType::AE && _data->dn[s].dt!=DataType::ODE && \
                  _data->dn[s].dt!=DataType::PDE && _data->dn[s].dt!=DataType::OPTIM && _data->dn[s].dt!=DataType::EIGEN

#define DEF_PAR_R(n,p,x) CHK_MSGR(_data->getPar(n,p,x),p,"Illegal parameter value")
#define DEF_PAR_B(n,p,x) CHK_MSGB(_data->getPar(n,p,x),p,"Illegal parameter value")
#define DDEF_PAR_R(n,p,x) CHK_MSGR(getPar(n,p,x),p,"Illegal parameter value")
#define DDEF_PAR_B(n,p,x) CHK_MSGB(getPar(n,p,x),p,"Illegal parameter value")
#define UNKNOWN_ARG(s) MSGR(s,"Unknown argument: "+_cmd->getToken())
#define NO_ARG(s) CHK_MSGR(nb_args==0,s,"No argument given for command.")
#define NO_VALUE_ARG(s) CHK_MSGR(nb==0,s,"No value given for argument "+_cmd->getToken())
#define HELP_ARG2 if (_nb_args>1) _rita->msg(_pr,"The second argument is unnecessary.");
#define HELP_ARG3(s) if (_nb_args>2) _rita->msg(_pr,"Help topic "+string(s)+" needs only one argument.");
#define FCT_NOT_DEFINED(p,f) CHK_MSGR(k==-1,_pr+p,"Function "+f+" not defined.")
#define FCT_ALREADY_DEFINED(p,f) CHK_MSGR(k==0,_pr+p,"Function "+f+" already defined.")
#define DEFAULT_KW  { if (_cmd->getChScript()==-1) break;                \
                      if (_cmd->checkEq()==string::npos) { _rita->msg("","Unrecognized command: "+_cmd->getToken()); break;  } \
                      ret = _rita->_calc->run(); break; }

