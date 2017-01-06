
bin_PROGRAMS = ec2d1
AUTOMAKE_OPTIONS = no-dependencies

ec2d1_SOURCES = main.cpp ScaledMF.cpp UpdateMF.cpp User.cpp
ec2d1_HEADERS = UpdateMF.h User.h

ec2d1_LDADD  = $(top_builddir)/src/libofeli.a

ec2d1dir   =	$(datadir)/ofeli/demos/electromagnetics/EddyCurrent2d-1
ec2d1_DATA =	main.cpp \
					ScaledMF.cpp \
					UpdateMF.cpp \
					User.cpp \
					UpdateMF.h \
					User.h \
					proj.dat \
					proj-0.m \
					proj-1.m \
					Makefile.standalone \
					ec2d1.cbp \
					README
dist_ec2d1_DATA =	main.cpp \
						ScaledMF.cpp \
						UpdateMF.cpp \
						User.cpp \
						UpdateMF.h \
						User.h \
						proj.dat \
						proj-0.m \
						proj-1.m \
						Makefile.standalone \
						ec2d1.cbp \
						README
