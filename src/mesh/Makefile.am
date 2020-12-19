
noinst_LIBRARIES = libmesh.a
AUTOMAKE_OPTIONS = no-dependencies
AM_CPPFLAGS = -std=c++1y

SUBDIRS = bamg metis

libmesh_a_SOURCES =	Domain.cpp \
			Edge.cpp \
			Element.cpp \
			getMesh.cpp \
			Grid.cpp \
			Material.cpp \
			Mesh.cpp \
			MeshAdapt.cpp \
			MeshExtract.cpp \
			MeshUtil.cpp \
			Node.cpp \
			Partition.cpp \
			saveMesh.cpp \
			Side.cpp

libmesh_a_LIBADD =	bamg/libbamg.a \
			metis/libmetis.a

