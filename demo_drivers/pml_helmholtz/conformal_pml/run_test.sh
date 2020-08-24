#!/bin/bash

max_adapt=1
n_pml=10
problem_params="--n_node 3 --max_adapt $max_adapt --wavenumber 4.0"

echo "Running conformal PML fish scattering"
mkdir -p RESLT
./conformal_pml $problem_params --use_dtn --dump_problem > RESLT/output.log
cp -a RESLT/soln_restart_file$max_adapt.dat exact_soln_restart_file.dat
rm -rf RESLT_DTN
mv RESLT RESLT_DTN

mkdir -p RESLT
./conformal_pml $problem_params --calculate_error_from_file --n_pml $n_pml > RESLT/output.log
rm -rf RESLT_PML
mv RESLT RESLT_PML

mkdir -p RESLT
./conformal_pml $problem_params --calculate_error_from_file --n_pml $n_pml --use_sfb > RESLT/output.log
rm -rf RESLT_PML_SFB
mv RESLT RESLT_PML_SFB
