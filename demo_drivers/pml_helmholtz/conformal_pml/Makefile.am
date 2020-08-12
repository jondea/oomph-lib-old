#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=conformal_pml

#--------------------------------------

# Sources for executable
conformal_pml_SOURCES =  conformal_pml.cc


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
conformal_pml_LDADD = -L@libdir@ -lhelmholtz -lpml_helmholtz -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------
