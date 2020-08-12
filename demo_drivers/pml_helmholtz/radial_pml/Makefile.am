#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=radial_pml

#--------------------------------------

# Sources for executable
radial_pml_SOURCES =  radial_pml.cc


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
radial_pml_LDADD = -L@libdir@ -lpml_helmholtz -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

#EXTRA_DIST += annular_meshes.h
