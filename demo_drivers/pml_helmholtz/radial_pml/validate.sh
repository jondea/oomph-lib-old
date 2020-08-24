#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=1


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for non-adaptive scattering
#-----------------------------------
cd Validation

echo "Running scattering using radial PML"
mkdir RESLT
echo 0.0 0.5 > exact_soln_coefficients.dat
echo 1.0 0.0 >> exact_soln_coefficients.dat
echo 0.2 0.7 >> exact_soln_coefficients.dat
../radial_pml --n_angular_elements 32 --n_node 4 --n_pml 10 > OUTPUT
echo "done"
echo " " >> validation.log
echo "Running scattering using radial PML" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > pml_scattering_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/pml_scattering_results.dat.gz   \
    pml_scattering_results.dat >> validation.log
fi

# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..

#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10

