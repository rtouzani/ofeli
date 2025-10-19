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
#define ALREADY_USED(c) if (c.count(name) && dn[name].active) { \
                           _rita->msg("","Name "+name+" already used for "+Type_st[d]); return 1; }
#define REPLACED(c) if (c.count(name) && dn[name].active) { \
                       _rita->msg("","\'"+name+"\' as "+Type_st[d]+" now used as "+Type_st[dt]); }

#define READLINE if (_cmd->readline(_pr+" ")<0) continue; \
                 key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw); \
                 if (key>=200) { _data->setDataExt(key); continue; }

#define CHK_PB(s) _data->dn[s].dt!=DataType::LS && _data->dn[s].dt!=DataType::AE && _data->dn[s].dt!=DataType::ODE && \
                  _data->dn[s].dt!=DataType::PDE && _data->dn[s].dt!=DataType::OPTIM && _data->dn[s].dt!=DataType::EIGEN

#define DEF_LS_NAME     "sys"
#define DEF_AE_NAME     "alg"
#define DEF_ODE_NAME    "ode"
#define DEF_PDE_NAME    "pde"
#define DEF_OPT_NAME    "opt"
#define DEF_EIG_NAME    "eig"
#define DEF_FCT_NAME    "fun"
#define DEF_GRID_NAME   "grd"
#define DEF_MESH_NAME   "msh"
#define DEF_VECTOR_NAME "v"
#define DEF_MATRIX_NAME "M"
#define DEF_TAB_NAME    "tab"

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
#define DEFAULT_KW(r)  { if (_cmd->getChScript()==-1) break; ret = r->_data->print(); if (ret) ret = r->_calc->run(); break; }

#define CHECK_VECT(s)          if (VectorLabel.count(s)==0) { _rita->msg("","Vector "+s+" undefined."); return 1; }
#define CHECK_ACTIVE(s)        if (dn[s].active==0) { _rita->msg("","Reference to non-active data "+s); return 2; }
#define CHECK_HVECT(s)         if (HVectorLabel.count(s)==0) { _rita->msg("","History vector "+s+" undefined."); return 1; }
#define CHECK_MATRIX(s)        if (MatrixLabel.count(s)==0) { _rita->msg("","Matrix "+s+" undefined."); return 1; }
#define CHECK_TAB(s)           if (TabLabel.count(s)==0) { _rita->msg("","Tabulation "+s+" undefined."); return 1; }
#define ILLEGAL_PREFIX         { _rita->msg("","Illegal name prefix for entity "+name); return 1; }

#define DEFAULT_LS_NAME(s)     if (s=="") s = DEF_LS_NAME+to_string(iLS)
#define DEFAULT_AE_NAME(s)     if (s=="") s = DEF_AE_NAME+to_string(iAE)
#define DEFAULT_ODE_NAME(s)    if (s=="") s = DEF_ODE_NAME+to_string(iODE)
#define DEFAULT_PDE_NAME(s)    if (s=="") s = DEF_PDE_NAME+to_string(iPDE)
#define DEFAULT_OPT_NAME(s)    if (s=="") s = DEF_OPT_NAME+to_string(iOpt)
#define DEFAULT_EIG_NAME(s)    if (s=="") s = DEF_EIG_NAME+to_string(iEig)
#define DEFAULT_FCT_NAME(s)    if (s=="") s = DEF_FCT_NAME+to_string(iFct)
#define DEFAULT_GRID_NAME(s)   if (s=="") s = DEF_GRID_NAME+to_string(iGrid)
#define DEFAULT_MESH_NAME(s)   if (s=="") s = DEF_MESH_NAME+to_string(iMesh)
#define DEFAULT_VECTOR_NAME(s) if (s=="") s = DEF_VECTOR_NAME+to_string(iVector)
#define DEFAULT_MATRIX_NAME(s) if (s=="") s = DEF_MATRIX_NAME+to_string(iMatrix)
#define DEFAULT_TAB_NAME(s)    if (s=="") s = DEF_TAB_NAME+to_string(iTab)
