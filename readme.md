Related publication: Angeloni et al. (2023) Dynamics of cortical contrast adaptation predict perception of signals in noise.
https://doi.org/10.1101/2021.08.11.455845

## Description

This code implements a Poisson GLM which estimates the multiplicative relationship between stimulus contrast and stimulus driven neural activity over time. Solving the function

$$\ln \lambda = \beta_{0} + x\beta_{1} + x \cdot C \beta_{2} + C \beta_{3}$$

where 
- $\lambda$ is the observed firing rate of the neuron.
- $x$ is the neurons linear response to a stimulus (stimulus drive).
- $x \cdot C$ is the multiplicative effect of the stimulus contrast and the stimulus drive (in theory, $C$ may be any time varying signal that is expected to affect neuronal gain, here we used stimulus contrast).
- $\beta$ are coefficients to be fit:
    
    - $\beta_{0}$ is the baseline firing rate (scalar)
    - $\beta_{1}$ is the stimulus drive (scalar; conceptually, how strongly the neuron is driven by changes in the stimulus)
    - $\beta_{2}$ is the multiplicative relationship between $C$ and the stimulus drive -- **we interpret this to be the gain of the neuron** (spline coefficients were used to fit this)
    - $\beta_{3}$ is the additive affect of contrast (spline)


There are three main conceptual variants of the model:

1. Symmetric time courses for low and high contrast. Uses the same set of splines to fit $\beta_{2}$ and $\beta_{3}$
2. Asymmetric gain control in each contrast. Uses two sets of splines to fit $\beta_{2}$ by dummy coding each set of coefficients for low and high contrast.
3. Asymmetric gain control and contrast sensitivity in each contrast. Uses four sets of splines to fit $\beta_{2}$ and $\beta_{3}$ by dummy coding each set of coefficients for low and high contrast. This is the model presented in the paper.

## Implementation

Practically, the model was fit in two steps:
1) To determine the stimulus drive $x$, we first fit the spectrotemporal receptive field of the neuron. The stimulus drive was then defined as the convolution of the STRF with the stimulus spectrogram.
2) Next, we fit the coefficients from the function above.

This entire fitting procedure was accomplished using 10-fold cross-validation, and fitting was performed using Poisson regression with `glmfit` in Matlab.

Fits to individual neurons were run on a computing cluster with precompiled Matlab code, see
`contrastGLM_acute.m` for compilation code (lines 34-36). The bash script called `submit_acute_run.sh` runs different versions of the model on the data stored on the cluster. If you wish to run the analysis from scratch, the raw data is provided on [Dryad](link) under `_data/_glm`.

The fitting procedure takes some time (approx 30 minutes per neuron). As such, we have provided the 
fits on [Dryad](link) in `_data/_glm/_res`. You may download the files, and point to the location of this directory in line 5 of `run_gcglm.m` to recreate the figures from the paper. For more details, see https://github.com/geffenlab/contrast_behavior/blob/main/readme.md.
