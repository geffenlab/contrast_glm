#!/bin/bash

# simulate neurons with different gain control, and different numbers of unique trials:
# 5 gain control levels, 2 different data conditions, 100 sims each = 1500 simulations

# parameter values
mu=30                      # mean value of stimulus
sd="[1,3]"                 # contrast
ex=( 100 )                 # exemplars
nt=( 5 )                   # number of repeats
w=0.3                      # STRF window (s)
pad=0                      # trial pad
lags=40                    # time lags to estimate gain
alpha=-1                   # regularization option
dummycode="false"          # dummy code each contrast separately
baseline=0.1               # baseline FR
gc=( -1.0 -0.5 0 0.5 1.0)  # amount of gain control (1 is perfect)
tau=( "[.5,.05]" "[.5,.5]" "[.05,.5]" )        # amount of adaptation (10 is effectively instant, .1 ~= 250ms)

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
	for ll in "${tau[@]}" ; do

	    tl=`echo $ll*10 | bc`
	    tl=${tl%.*}

	    # for each neuron
	    for ii in {1..100} ; do

		fn="$(printf "sim_lohi-2tau_gc%03d_ex%02d_nt%03d_tau%03d_run%03d.mat" $gl ${ex[$kk]} ${nt[$kk]} $tl $ii)"
		qsub -l h_vmem=8G -l short -N lohi-notrans $filepath/run_contrastGLM_sim_ind.sh $MATLAB_DIR \
		     $respath/$fn $mu $sd ${ex[$kk]} ${nt[$kk]} $w $lags $pad $alpha $dummycode $baseline $jj $ll
		echo "submit $fn $gl ${ex[$kk]} ${nt[$kk]} $tl $ii with 8G ram"
		
	    done
	done
    done
done


# to submit a test run:
qsub -l h_vmem=8G -l short -N lohi-notrans $filepath/run_contrastGLM_sim_ind.sh $MATLAB_DIR \
     "test.mat" $mu $sd ${ex[$kk]} ${nt[$kk]} $w $lags $pad $alpha $dummycode $baseline $jj $ll

watch 'qstat | wc -l ; qstat'

rsync -avhP \
      angelonc@cubic-login.uphs.upenn.edu:/cbica/home/angelonc/comp_space/contrast_glm/_simulation/_res/*tau*.mat \
      /Users/chris/chris-lab/projects/contrast_glm/_simulation/_res/
