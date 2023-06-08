#!/bin/bash

# run parameters
a=10
gc=1
gl=`echo $gc*10 | bc`
gl=${gl%.*}
alpha=0
al=`echo $alpha*100 | bc`
al=${al%.*}

filepath="/cbica/home/angelonc/comp_space/contrast_glm_v2.0/_simulation"

# to test run:
# fn="$(printf "sim_a%02d_g%02d_alpha%02d-run%s.mat" $a $gl $al test)"
# $filepath/run_contrastGLM2_sim.sh $MATLAB_DIR $fn $a $gc $alpha
# /usr/bin/time -pv $filepath/run_contrastGLM2_sim.sh $MATLAB_DIR $fn $a $gc $alpha

for i in {1..100} ; do

    fn="$(printf "%s/_res/sim_a%02d_g%02d_alpha%02d-run%03d.mat" $filepath $a $gl $al $i)"
    # echo $filepath/run_contrastGLM2_sim.sh $MATLAB_DIR $fn $a $gc $alpha
    echo "submit $fn"
    qsub -l h_vmem=8G $filepath/run_contrastGLM2_sim.sh $MATLAB_DIR $fn $a $gc $alpha
    
done
