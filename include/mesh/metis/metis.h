/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.1 1998/11/27 17:59:21 karypis Exp $
 */

#ifndef __METIS_H
#define __METIS_H

#include <stdio.h>
/*#ifdef __STDC__*/
#include <stdlib.h>
/*#else
#include <malloc.h>
#endif*/
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include "mesh/metis/defs.h"
#include "mesh/metis/struct.h"
#include "mesh/metis/metis_macros.h"
#include "mesh/metis/rename.h"
#include "mesh/metis/proto.h"

#endif
