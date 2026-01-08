# ==============================================================================
#
#                                    O  F  E  L  I
#
#                           Object  Finite  Element  Library
#
# ==============================================================================
#
#   Copyright (C) 1998 - 2024 Rachid Touzani
#
#   This file is part of OFELI.
#
#   OFELI is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   OFELI is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with OFELI. If not, see <http://www.gnu.org/licenses/>.
#
# ==============================================================================

target_sources (${PROJECT_NAME} PRIVATE
                balance.c
                bucketsort.c
                ccgraph.c
                coarsen.c
                compress.c
                debug.c
                estmem.c
                fm.c
                fortran.c
                frename.c
                graph.c
                initpart.c
                kmetis.c
                kvmetis.c
                kwayfm.c
                kwayrefine.c
                kwayvolfm.c
                kwayvolrefine.c
                match.c
                mbalance2.c
                mbalance.c
                mcoarsen.c
                memory.c
                mesh_metis.c
                meshpart.c
                mfm2.c
                mfm.c
                mincover.c
                minitpart2.c
                minitpart.c
                mkmetis.c
                mkwayfmh.c
                mkwayrefine.c
                mmatch.c
                mmd.c
                mpmetis.c
                mrefine2.c
                mrefine.c
                mutil.c
                myqsort.c
                ometis.c
                parmetis.c
                pmetis.c
                pqueue.c
                refine.c
                separator.c
                sfm.c
                srefine.c
                stat.c
                subdomains.c
                timing.c
                util.c
               )
