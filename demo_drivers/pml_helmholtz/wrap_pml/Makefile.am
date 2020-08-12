#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=wrap_pml_test

#--------------------------------------

# Sources for executable
wrap_pml_test_SOURCES =  wrap_pml_test.cc


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
wrap_pml_test_LDADD = -L@libdir@ -lpml_helmholtz -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

#EXTRA_DIST += annular_meshes.h
