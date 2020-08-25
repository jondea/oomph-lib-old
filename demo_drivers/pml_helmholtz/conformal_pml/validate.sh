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

max_adapt=1
n_pml=10
problem_params="--n_node 3 --max_adapt $max_adapt --wavenumber 4.0"

echo "Running conformal PML tests"
mkdir RESLT
gzip -dkf ../validata/exact_soln_restart_file.dat.gz
mv ../validata/exact_soln_restart_file.dat ./
../conformal_pml $problem_params --calculate_error_from_file --n_pml $n_pml > OUTPUT
echo "done"
echo " " >> validation.log
echo "Running conformal PML tests" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/relative_error1.dat > conformal_pml_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/conformal_pml_results.dat.gz   \
    conformal_pml_results.dat >> validation.log
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

