#!/bin/bash

filepath="/cbica/home/angelonc/comp_space/contrast_glm"

doneCnt=0
runCnt=0
win=0.3        # STRF window length (s)
gwin=1         # gain window length (s)
alpha=-1       # alpha parameter (regularization option)
dc=( 0 1 2 )   # dummy coding (0 = none, 1 = dummy beta2, 2 = dummy beta2&3)
cv=10          # cross validation folds
spline="[3,3]" # spline basis: [n_knots,degree]
spll="3-3"

#f_test="$filepath/../_data/K184_2021-04-21_13-48-12_009_025_mu.mat"
#qsub -l h_vmem=16G -N acute_2tau_spline $filepath/run_contrastGLM_acute.sh $MATLAB_DIR \
#     $f_test "qsub_test.mat" $win $gwin $alpha 2 10 $spline  

for f in $filepath/../_data/K184*.mat
do
    #echo "Input file: $f"
    fn="${f##*/}"

    for i in "${dc[@]}"; do

	fout="$(printf "%s/_res/acute_spline-%s_dc%d_%s" $filepath $spll $i $fn)"
	#fout="$filepath/_res/acute_$fn_dc$i"
	#echo "Output file: $fout"

	if test -f "$fout"; then
	    ((doneCnt=doneCnt+1))
	    echo "$fout exists ... $doneCnt"
	else
	    ((runCnt=runCnt+1))
	    echo "$fout doesnt exist, submitting job ... $runCnt"
	    qsub -l h_vmem=16G -N acute_2tau_spline $filepath/run_contrastGLM_acute.sh $MATLAB_DIR \
	    	   $f $fout $win $gwin $alpha $i $cv $spline
	fi

    done

done

echo "$doneCnt done. $runCnt to run."

watch 'qstat | wc -l ; qstat'


rsync -avhP \
      angelonc@cubic-login.uphs.upenn.edu:/cbica/home/angelonc/comp_space/contrast_glm/_res/acute_*.mat \
      /Users/chris/data/gain_opto/glm_res/


## submit failed run
dc=1
f_fail="$filepath/../_data/K184_2021-04-21_16-53-20_005_018_su.mat"
f_res="$(printf "%s/_res/acute_dc%d_%s" $filepath $dc "K184_2021-04-21_16-53-20_005_018_su.mat")"
qsub -l h_vmem=16G -N glm_acute_2tau $filepath/run_contrastGLM_acute.sh $MATLAB_DIR \
     $f_fail $f_res $win $gwin $alpha $dc 10 $spline
