#!/bin/bash

# simulate neurons with different gain control, and different numbers of unique trials:
# 5 gain control levels, 2 different data conditions, 100 sims each = 1500 simulations

# parameter values
mu=30                      # mean value of stimulus
sd="[1,3]"                 # contrast
ex=( 5 )                   # exemplars
nt=( 100 )                 # number of repeats
w=0.3                      # STRF window (s)
pad=0                      # trial pad
lags=40                    # time lags to estimate gain
alpha=-1                   # regularization option
dummycode=( 0 1 2 )        # dummy code each contrast separately
cv=10                      # cv runs (not used)
spline_basis="[7,3]"       # spline knots and degree
baseline=0.1               # baseline FR
gc=( -1.0 -0.5 0 0.5 1.0)  # amount of gain control (1 is perfect)
tau=( "[.5,.05]" "[.5,.5]" "[.05,.5]" )        # amount of adaptation (10 is effectively instant, .1 ~= 250ms)
tl=( "50-5" "50-50" "5-50" )

filepath="/cbica/home/angelonc/comp_space/contrast_glm/_simulation"
respath="$filepath/_res"

# for each gain level
for jj in "${gc[@]}" ; do

    # gain label format
    gl=`echo $jj*10 | bc`
    gl=${gl%.*}

    # for each exemplar
    for kk in "${!ex[@]}" ; do

	# for each tau
	for ll in "${!tau[@]}" ; do

	    for mm in "${dummycode[@]}" ; do

		# for each neuron
		for ii in {1..100} ; do

		    fn="$(printf 'sim_2tau-spline-5ex_gc%03d_ex%02d_nt%03d_tau%se-2_dc%d_run%03d.mat' $gl ${ex[$kk]} ${nt[$kk]} ${tl[$ll]} $mm $ii)"
		    qsub -l h_vmem=8G -l short -N spline-2tau $filepath/run_contrastGLM_sim.sh $MATLAB_DIR \
		    	 $respath/$fn $mu $sd ${ex[$kk]} ${nt[$kk]} $w $pad $lags $alpha $mm $cv $spline_basis $baseline $jj ${tau[$ll]}
		    echo "submit $fn $gl ${ex[$kk]} ${nt[$kk]} ${tau[$ll]} $mm $spline_basis $ii with 8G ram"
		    
		done
	    done
	done
    done
done

watch 'qstat | wc -l ; qstat'
	

# to submit a test run:
qsub -l h_vmem=8G -l short -N spline-2tau $filepath/run_contrastGLM_sim.sh $MATLAB_DIR \
     "test.mat" $mu $sd ${ex[0]} ${nt[0]} $w $pad $lags $alpha ${dummycode[0]} $cv $spline_basis $baseline 1 ${tau[0]}


rsync -avhP \
      angelonc@cubic-login.uphs.upenn.edu:/cbica/home/angelonc/comp_space/contrast_glm/_simulation/_res/*2tau*.mat \
      /Users/chris/data/gain_behavior/_simulations/
