#!/bin/bash

# simulate neurons with different gain control, and different numbers of unique trials:
# 5 gain control levels, 3 different data conditions, 100 sims each = 1500 simulations

# gain values
gc=( -1.0 -0.5 0 0.5 1.0)

# exemplars
ex=( 5 50 100 )

# number of repeats
nt=( 100 10 5 )

filepath="/cbica/home/angelonc/comp_space/contrast_glm/_simulation"
respath="$filepath/_res"
parampath="$filepath/_params"

for j in "${gc[@]}" ; do

    # gain label format
    gl=`echo $j*10 | bc`
    gl=${gl%.*}

    for k in "${!ex[@]}" ; do

	#printf "%s\t%s\n" "$k" "${ex[$k]}"
	
	for i in {1..100} ; do

	    fn="$(printf "sim_noregweighted_gc%03d_ex%02d_nt%03d_run%03d.mat" $gl ${ex[$k]} ${nt[$k]} $i)"
	    fnp="$(printf "params_gc%03d_ex%03d_nt%03d.txt" $gl ${ex[$k]} ${nt[$k]})"

	    qsub -l h_vmem=12G -l short $filepath/run_contrastGLM_sim.sh $MATLAB_DIR $respath/$fn $parampath/$fnp
	    echo "submit $fn $fnp $gl ${ex[$k]} ${nt[$k]} $i with 16G ram"
  
	done
    done
done


# to submit a test run:
qsub -l h_vmem=12G $filepath/run_contrastGLM_sim.sh $MATLAB_DIR $respath/test.mat $parampath/params_gc010_ex050_nt010.txt

watch 'qstat | wc -l ; qstat'

rsync -avhP \
      angelonc@cubic-login.uphs.upenn.edu:/cbica/home/angelonc/comp_space/contrast_glm/_simulation/_res/*ex*.mat \
      /Users/chris/chris-lab/projects/contrast_glm/_simulation/_res/
