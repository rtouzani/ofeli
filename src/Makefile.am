## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libofeli.a
AUTOMAKE_OPTIONS = no-dependencies
SUBDIRS = equations io linear_algebra mesh post shape_functions solvers util

libofeli_a_SOURCES =

EQUATIONS_OBJS = \
			equations/DG.o

CL_OBJS = \
			equations/cl/Muscl.o \
			equations/cl/Muscl1D.o \
			equations/cl/Muscl2DT.o \
			equations/cl/Muscl3DT.o \
			equations/cl/LCL1D.o \
			equations/cl/LCL2DT.o \
			equations/cl/LCL3DT.o \
			equations/cl/ICPG1D.o \
			equations/cl/ICPG2DT.o \
			equations/cl/ICPG3DT.o

ELECTROMAGNETICS_OBJS =	\
			equations/electromagnetics/BiotSavart.o \
			equations/electromagnetics/EC2D1T3.o \
			equations/electromagnetics/EC2D2T3.o \
			equations/electromagnetics/EC2DB2T3.o \
			equations/electromagnetics/HelmholtzBT3.o

FLUID_OBJS = \
			equations/fluid/NSP2DQ41.o \
			equations/fluid/TINS2DT3B.o \
			equations/fluid/NS2DT3BT3.o

INTERFACE_OBJS = \
			equations/interface/FastMarching2D.o \
			equations/interface/FMM.o \
			equations/interface/FMM2D.o \
			equations/interface/FMM3D.o \
			equations/interface/FMMSolver.o \
			equations/interface/Heap.o \
			equations/interface/IPoint.o

LAPLACE_OBJS = \
			equations/laplace/Laplace1DL2.o \
			equations/laplace/Laplace1DL3.o \
			equations/laplace/Laplace2DT3.o \
			equations/laplace/Laplace3DT4.o \
			equations/laplace/Laplace2DMHRT0.o \
			equations/laplace/Laplace2DFVT.o \
			equations/laplace/LaplaceDG2DP1.o \
			equations/laplace/SteklovPoincare2DBE.o

POROUS_OBJS = \
			equations/porous/WaterPorous2D.o

SOLID_OBJS = \
			equations/solid/Bar2DL2.o \
			equations/solid/Beam3DL2.o \
			equations/solid/Elas2DT3.o \
			equations/solid/Elas2DQ4.o \
			equations/solid/Elas3DH8.o \
			equations/solid/Elas3DT4.o

THERM_OBJS = \
			equations/therm/DC1DL2.o \
			equations/therm/DC2DT3.o \
			equations/therm/DC2DT6.o \
			equations/therm/DC3DAT3.o \
			equations/therm/DC3DT4.o \
			equations/therm/PhaseChange.o

IO_OBJS = \
			io/FFI.o \
			io/IOField.o \
			io/IPF.o \
			io/Prescription.o \
			io/Tabulation.o \
			io/XMLParser.o \
			io/saveField.o

FPARSER_OBJS = \
			io/fparser/fparser.o \
			io/fparser/fpoptimizer.o

XMLSP_OBJS = \
			io/xmlsp/xmlsp.o

LINEAR_ALGEBRA_OBJS = \
			linear_algebra/GraphOfMatrix.o

MESH_OBJS = \
			mesh/Domain.o \
			mesh/Edge.o \
			mesh/Element.o \
			mesh/getMesh.o \
			mesh/Grid.o \
			mesh/Material.o \
			mesh/Mesh.o \
			mesh/MeshAdapt.o \
			mesh/MeshExtract.o \
			mesh/MeshUtil.o \
			mesh/Node.o \
			mesh/Partition.o \
			mesh/saveMesh.o \
			mesh/Side.o

BAMG_OBJS =	\
			mesh/bamg/bamg.o \
			mesh/bamg/Mesh2.o \
			mesh/bamg/MeshGeom.o \
			mesh/bamg/MeshQuad.o \
			mesh/bamg/MeshRead.o \
			mesh/bamg/MeshWrite.o \
			mesh/bamg/Meshio.o \
			mesh/bamg/Metric.o \
			mesh/bamg/QuadTree.o \
			mesh/bamg/SetOfE4.o

METIS_OBJS = \
			mesh/metis/coarsen.o \
			mesh/metis/fm.o \
			mesh/metis/initpart.o \
			mesh/metis/match.o \
			mesh/metis/ccgraph.o \
			mesh/metis/memory.o \
			mesh/metis/pmetis.o \
			mesh/metis/pqueue.o \
			mesh/metis/refine.o \
			mesh/metis/util.o \
			mesh/metis/timing.o \
			mesh/metis/debug.o \
			mesh/metis/bucketsort.o \
			mesh/metis/graph.o \
			mesh/metis/stat.o \
			mesh/metis/kmetis.o \
			mesh/metis/kwayrefine.o \
			mesh/metis/kwayfm.o \
			mesh/metis/balance.o \
			mesh/metis/ometis.o \
			mesh/metis/srefine.o \
			mesh/metis/sfm.o \
			mesh/metis/separator.o \
			mesh/metis/mincover.o \
			mesh/metis/mmd.o \
			mesh/metis/mesh_metis.o \
			mesh/metis/meshpart.o \
			mesh/metis/frename.o \
			mesh/metis/fortran.o \
			mesh/metis/myqsort.o \
			mesh/metis/compress.o \
			mesh/metis/parmetis.o \
			mesh/metis/estmem.o \
			mesh/metis/mpmetis.o \
			mesh/metis/mcoarsen.o \
			mesh/metis/mmatch.o \
			mesh/metis/minitpart.o \
			mesh/metis/mbalance.o \
			mesh/metis/mrefine.o \
			mesh/metis/mutil.o \
			mesh/metis/mfm.o \
			mesh/metis/mkmetis.o \
			mesh/metis/mkwayrefine.o \
			mesh/metis/mkwayfmh.o \
			mesh/metis/mrefine2.o \
			mesh/metis/minitpart2.o \
			mesh/metis/mbalance2.o \
			mesh/metis/mfm2.o \
			mesh/metis/kvmetis.o \
			mesh/metis/kwayvolrefine.o \
			mesh/metis/kwayvolfm.o \
			mesh/metis/subdomains.o

POST_OBJS = \
			post/Estimator.o \
			post/Reconstruction.o

SHAPE_FUNCTIONS_OBJS = \
			shape_functions/Hexa8.o \
			shape_functions/Line2.o \
			shape_functions/Line2H.o \
			shape_functions/Line3.o \
			shape_functions/Penta6.o \
			shape_functions/Quad4.o \
			shape_functions/Tetra4.o \
			shape_functions/Triang3.o \
			shape_functions/Triang6S.o

SOLVERS_OBJS =	\
			solvers/BSpline.o \
			solvers/EigenProblemSolver.o \
			solvers/ODESolver.o \
			solvers/OptimNM.o \
			solvers/OptimSA.o \
			solvers/OptimTN.o \
			solvers/OptimPG.o \
			solvers/OptSolver.o \
			solvers/Solver.o \
			solvers/TimeStepping.o

UTIL_OBJS = \
			util/Gauss.o

libofeli_a_LIBADD =	\
			$(EQUATIONS_OBJS) \
			$(CL_OBJS) \
			$(ELECTROMAGNETICS_OBJS) \
			$(FLUID_OBJS) \
			$(INTERFACE_OBJS) \
			$(LAPLACE_OBJS) \
			$(POROUS_OBJS) \
			$(SOLID_OBJS) \
			$(THERM_OBJS) \
			$(IO_OBJS) \
			$(FPARSER_OBJS) \
			$(XMLSP_OBJS) \
			$(LINEAR_ALGEBRA_OBJS) \
			$(MESH_OBJS) \
			$(BAMG_OBJS) \
			$(METIS_OBJS) \
			$(POST_OBJS) \
			$(SHAPE_FUNCTIONS_OBJS) \
			$(SOLVERS_OBJS) \
			$(UTIL_OBJS)
