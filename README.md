# NaPuCco: Group-level inference of fMRI data through non-parametric combination

**A tool for estimating the statistical Power while varying the population size via NPC**

A project developed by Matteo Lionello and Luca Cecchetti, Social and Affective NEuroscience, IMT Lucca, Italy.
The code is part of an on-going project planned as an extensions of the work proposed by [Winkler A. M. et al. (2016). *Non-parametric combination and related permutation tests for neuroimaging*](https://doi.org/10.1002/hbm.23115), aiming at validating power analysis when non-parametric combination is applied between subjects.


## Table of contents
1. [General workflow](#schema)
2. [OOP basic implementation](#oop)
3. [Main Functions](#func)
4. [How to use](#howto)
      - [Step 1: Generating volumes](#genvol)
      - [Step 2: Simulating Experiments](#simexp)
      - [A note on 'volumes' and 'subjects' terms and uses](#volsubj)
5. [Pipeline examples](#examples)
6. [Folder structure](#foldstruct)

## General workflow <a name="schema"/>

This project aims at introducing a power analysis when performing non-parametric combination between subjects in an f-MRI setup.

**DATA**: N participants, M features(regressors), L voxels, T timestemps, J permutations<br />
regression for each subjects between one voxels Y = V(x, y, z, T), and input features X=(M, T).<br />
R2 collection for each voxel, participants, and permutation [L, N, J+1].<br />
_output_:<br />
P-values calcualtion [L, N, J+1] <br />

*Non-parametric combination*:<br />
P-values (u-values winkler) combining function (e.g. Fisher Methods) across subjects [J+1, L]

**(A) VOXELWISE correction**: concatenation of voxels maxima from each partial null Hypothesis [J]<br />
Comparison with non-permutated voxels [L=1] obtaining p-values familiwise correction.

**(B) CLUSTERBASED correction**:<br />
set cluster forming threhsolds on u-values and get binary matrix (J+1, L)<br />
retrieve coordinates for a 3d volume, and apply topographic clustering rules<br />
max K against null hypothesis.

**OUTPUT**
Clusterbased and voxelwise errors are repeated K times for each experiments, the results (pass or fail) are then averaged across experiments to calculate the power and fpr.

The work is split in two parts:
  - subject generation
  - experiment simulation

## OOP EASY PEASY STRAIGHT FOWARD EXAMPLE: <a name="oop"/>
      % initialize experiment design with 6 regressors and 120 timepoints
      napuccoObj = Napucco(6, 120)
      
      % generate 120 volumes with 27 effected and 80K non-effected voxels. The effected volumes are generrated by sampling their R2 target from a N(0.07, 0.02) distribution
      napuccoObj = napuccoObj.generate_volumes(120, 'nb_effvx', 27, 'nb_noneffvx', 80000, 'r2_target', [0.07, 0.02]);
      
      % Voxelwise and Clusterbased correction are calulated for 100 (default) experiments, each combining (NPC Fisher) 5 to 80 participants generated from the generated ones
      napuccoObj = napuccoObj.run_experiments('both', [])
      
      % FPR and Power analysis is averaged across the results
      Napucco.compute_result(napuccoObj)
      
## Main functions: <a name="func"/>
  - **vc_correction**: core function used for both subjects generation and experiment simulation.
  - **subj_sampling**: function called after generating subjects, it samples the subjects generated and it calls vc_correction for running the simulated experiments.

## How to use: <a name="howto"/>
      
### Step 1: Generating volumes <a name="genvol"/>
Volumes carrying effected and non-effected voxels can be generated calling vc_correction. Fo

  - **_vc_correction(nb_subjects, nb_features, nb_outfeat, nb_timepoints, r2_target)_**: it will generate nb_subjects volumes with nb_features regressors and nb_timepoints timepoints. the volume is composed by nb_outfeat voxels. For each volume, a R2 target effect for all its voxel is sampled from $\mathcal{N}$(r2mu, r2std) when **r2_target=[r2mu, r2std]**. [DEPRECATED: If r2_target is given empty, volumes will be still considered as effected but composed by voxels non-effected; see below to build non-effected volumes.]
  
  - **_vc_correction(nb_subjects, nb_features, nb_outfeat, nb_timepoints, r2_target, [] 'greater_dist', 1)_** with a large **nb_subjects**>400 when **greater_dist** is set to 1
  the volumes will be generated to fit a broad distribution of R2. By doing in this way, it will not be necessary to generate different set of volumes to investigate different different R2. The distribution is given in the function **getR2_fromdist()** inside *generators/create_subject.m_
  
  - By default the variability of the R2 in each volume is given by the variability passed by the parameter **_vc_correction(..., 'threshold', 0.002)_**.
  Variance within each volume can be controlled via **_vc_correction(..., 'vox_variability', 0)_**.
  
  - **_vc_correction(nb_subjects, nb_features, 0, [], 'sec_cluster_voxels', nb_outfeat)_** in this case the volume will be composed by 0 effected voxels and nb_outfeat of non-effected voxels (voxels not carrying effects). The generated volume will be stored in _subj_XXX.m_ file in the variable datanull.
  
  - **_vc_correction(nb_subjects, nb_features, nb_outfeat_eff, nb_timepoints, r2_target, 'sec_cluster_voxels', nb_outfeat_noneff)_** will generate nb_subjects volumes files , _subj_{1...nb_subjects}.m_, containing both variable _data_ and _datanull_
  
  - **_vc_correction(... 'numb_permutations', 2000)_** it will generate 2000 timepoints shuffling synchronized with the other volumes and other voxels within the same generation. The shuffling schema is stored in _hyparams.mat_ file

**IF GENERATION IS EARLY INTERRUPTED** it can be resumed via **_resume_generation(path)_**

### Step 2: Simulating Experiments <a name="simexp"/>
Experiments results include Power and Fpr analysis for both Voxelwise and Cluster based correction. The result of each experiment is saved in 
  
  - **_subj_subsampling(correction, homepath, nb_subj, nb_rep, r2)_** It applies the given **correction** to **nb_subj** combined subjects randomly sampled from **homepath** and it repeats the extraction and the correction **nb_rep** times. The results (pass or fail) are saved in _grouping_ folder. If **r2** dist is not empty, such as **r2=[r2_mean, r2_std]**, effected volumes are sampled in order to match the distribution $\mathcal{N}$(**r2_mean, r2_std**).
  
  - **_subj_subsampling(..., 'fpr', 0)_** it computes only the power. Non-effected volumes path is needed to be passed as second item in homepath string array.
    
  - **_subj_subsampling(..., 'power', 0)_** it computes only the fpr. *homepath* is assumed to contain only the non-effected volume path. Otherwise, the first item in *homepath* array will be ignored.

 - **_subj_subsampling(...,  'subj_range' = [5: 5: nb_subj])_** the output will match the combinations of subjects grouped according to the given interval.

 - **_subj_subsampling(...,  'noneff_range_vx', [500:500:10000])_** The experiments will be simulated while varying the size of non-effected volumes (this is a mandatory step to identify the number of non-effected voxels needed for statistics to converge).
 
 -  **_subj_subsampling(...,  'padding_range', [10:1:20])_** Same as before but varying the padding width of non-effected voxel in the case of power in clusterbased. The size of non-effected volumes will be checked according to the size of effected volumes and the amount of available non-effected voxels, and calcualted in _utils/resizecubicpadding.m_

 - **_subj_subsampling(...,  'noneff_range_clus', [10:1:20])_** Same as before but varying the padding width of non-effected voxel in the case of FPR in clusterbased.

  - **_subj_subsampling(..., 'stack_mem', 600)_** When iterating the volumes sampling, up to 600 volumes are kept (based on a ranking queue with respect their usage) in the RAM for the sampling process without need to reload them at every repetition of the experiment.

  - **_compute_mean_res(path)_** it generates the average results for all the experiments simulated with **subj_subsampling**. The path must in the format 'data/inXX_tXX/results/XXXXX' and to contain 'grouping/'.

**IF EXPERIMENT IS EARLY INTERRUPTED** it can be resumed via **_subj_subsampling(..., 'resfolder', path)_** while keeping the other arguemnts explictly identical to the other experiment repetitions. path must in the format 'data/inXXX_tXXXX/results/XXXXXX'.

### A note on 'volumes' and 'subjects' terms and uses <a name="volsubj"/>
The terms subjects and volumes can lead to some missleading. subjects are created by concatenating effected with non-effected volumes. In the case these volumes have been generated separatly they will not share the same permutation schema, making the permutation not synchornized between voxels (although still synch between subjects). **TODO** In order to generate a new volume with a specific permutation schema, please use **_vc_correction(..., 'prev_settings', prev_settings)_**.

## Pipeline examples: <a name="examples"/>
Here follows a few initial examples of how to combine the different part of the code. Further usages can be found in *./scripts/* folder.

    voxelwise correction:
    vc_correction('none', 120, 6, 1, 120, [0.079, 0.017])              # generate effected voxel (1 effected voxels for 120 subjects - mean 0.079, std 0.017)
    vc_correction('none', 400, 6, 0, 120, [], ..                       # generate non-effected voxels (8000 non-eff. for 400 subjects)
      'sec_cluster_voxels', 8000)
    outfolder = subj_subsampling('voxelwise', ...                      # sample 80 volumes from eff and non-eff ones, and calcuate the outcome of 150 experiments
      'data/in006_t0120/subjects/r2_079_017_230101_080000', ...
      80, 150, [], 'subj_range' = [5: 5: 80], ...                      # r2 is not needed to be specified as the pewvious generated folder will be used to be randomly sampled
      'nullfold', 'data/in006_t0120/subjects/nulldst_230101_080002')   # in the case the effected voxels were generated with 'grtdst' arguments, it must be selected the r2 range
    compute_mean_res(outfolder)                                        # average the individual results to get the power of the voxelwise correction
    
    clusterbased correction:
    vc_correction('none', 120, 6, 343, 120, [0.079, 0.017])            # generate effected volumes (343 effected voxels for 120 subjects - mean 0.079, std 0.017)
    outfolder = subj_subsampling('clusterbased',...
      'data/in006_t0120/subjects/r2_079_017_230102_090000', ...
      80, 150, [], 'subj_range' = [5: 5: 80], ...                      # r2 similarly treated as above
      'nullfold', 'data/in006_t0120/subjects/nulldst_230101_080002')   # the non-effected volumes generated early will be used
    compute_mean_res(outfolder)

    Oneshot generation and experiment:
    vc_correction('voxelwise', 120, 6, 1, 120, [0.079, 0.017], 'sec_cluster_voxels', 8000)
    
## Folder structure <a name="foldstruct"/>

    ./napucco/
    ├── vc_correction.m
    ├── subj_sampling.m
    ├── resume_generation.m
    ├── ...
    ├── scripts/
    │   └── ...
    ├── generators/
    │   ├── create_subjects.m
    │   └── ...
    ├── fwc/
    │   └── ...
    ├── data/
    │   ├── inYYY_tYYYY/                                       # This is the home folder for all the simulation sharing the same number of regressors and number of timepoints
    │   │   ├── results/                                       # experiments
    │   │   │   ├── r2_muR2_stdR2_{K,V}DATE_TIME/              # experiments results for given R2 target sampled from population, cluster or voxel wise correction
    │   │   │   │   ├── hyparams.mat
    │   │   │   │   ├── opts.mat
    │   │   │   │   ├── grouping/                              # this folder contains all the individual outcome (pass or fail) of all the experiment repetition {results_XXX.mat}
    │   │   │   │   └── groupingavg/                           # this folder contains power.mat file containing the fpr or power results for the set of experiments
    │   │   │   └── ...
    │   │   └── subjects/                                      # volumes
    │   │       ├── {nulldst,grtdst,r2_muR2_stdR2}_DATE_TIME/  # collection of generated volumes sampled from the same R2 distribution 
    │   │       │   ├── hyparams.mat
    │   │       │   ├── dist_settings.mat
    │   │       │   ├── histR2{null}.csv                      # R2 information storage for volumes extraction during experiments
    │   │       │   ├── subj_XXXX.mat                         # R2 and p-values genereated for each subject/volume.
    │   │       │   └── ...
    │   │       └── ...
    │   └── ...
    └── ...
