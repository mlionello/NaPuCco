# Estimating the statistical Power while varying the population size under NPC

A project developed by Matteo Lionello and Luca Cecchetti, Social and Affective Neuroscience (SANe), IMT Lucca, Italy.
The code is part of an on-going project planned as an application of the work proposed by [Winkler A. M. et al. (2016). *Non-parametric combination and related permutation tests for neuroimaging*](https://doi.org/10.1002/hbm.23115) to validate non-parametric combination across participants in a one-sample multi-column encoding fmri study setup, at the same time, by providing a tool for power analysis.

## Table of contents
1. [General workflow](#schema)
2. [Main Functions](#func)
3. [How to use](#howto)
4. [Pipeline examples](#examples)
5. [Folder structure](#foldstruct)

## General workflow <a name="schema"/>

This project aims at introducing a power analysis when performing non-parametric combination between subjects in an f-MRI setup.
The work is split in three parts:
  - volumes generation
  - experiment simulation
  - experiment-wise analysis
    
### GLM-Volume Generator:​

    Given a design matrix X #timepoints by #features​;
    given a R2target value sampled from a given distribution:
      xi ~ N(R2_µ, R2_std);
    for each volume i, generate voxel j such that:​
    sample
      y ~ N(0, stdVx) + U(beta1, beta2) * X,
    by sampling random beta until: ​xi - regress(y, X) < thr​
    (to generate non-effected volumes no “< thr“ constraints).​
    For each voxel a R2 null-distribution is calculated via timepoint-shuffling.​
    Then, for each voxel, get the parametric u-value via:​
      1 - Beta(R2ij | #features/2, (#timepoints - #features - 1) / 2 )​


### Subject Creation:​

      Voxelwise: concatenating 1 effected voxel with non-effected volumes (true positive); 
            only non-effected volumes (false positives);​
      Clusterbased: 3d-wrap effected and 3d-padding with non-effected (true positive);
            3d-wrap non-effected (false positives).

### Experiment Simulation:​
      
      U-values combination via fisher:​
            combined_subs = - 2 * Σ{i = 1..N}(log(uvalues_subi))​

#### Voxelwise:​

    Null-distribution from max statistic across all permutations.
    Each experiment passes if: p_value < alpha for the effected voxel (true positive);
    and any p_value < alpha for the non-effected volume (false positive).​

#### Clusterbased:​

    Get maximum cluster size from those voxels passing the cluster forming threshold:​
    1 - χ2(comb_sub | #sub * 2) < cft​
    The experiments passes if: p_value < alpha for the mixed volume (true positive) and any p_value < alpha for the non-effected 3d-warp (false positive).​
​
### Power estimation

### FPR estimation

## Main functions: <a name="func"/>
  - **volume_generator**: function used to generate effected and non-effected volumes of voxels.
  - **subj_sampling**: function called after generating volumes,
    it concatenates sampled volumes into subjects and it calls vc_correction for running the simulated experiments.

## How to use: <a name="howto"/>

Please, visit the [Wiki Documentation](https://github.com/mlionello/NPC/wiki)

## Pipeline examples: <a name="examples"/>
Here follows a few initial examples of how to combine the different part of the code. Further usages can be found in the wiki documentation.

    voxelwise correction:
    generate_volumes(300, 6, 120, 1, [0.079, 0.017], ... 
        'greater_dist', 1 )                                           # generate effected voxel (1 effected voxels for 300 subjects - mean 0.079, std 0.017)
    generate_volumes(300, 6, 120, 0, [], ...
        'nb_noneffvx', 90000)                                         # generate non-effected voxels (90000 non-eff. for 300 subjects)
    generate_volumes(300, 6, 120, 27, [[0.079, 0.017], ...            # combination of the previous two commands
        'nb_noneffvx', 90000, ...
        'greater_dist', 1 )

    outfolder = subj_subsampling('voxelwise', ...                      # sample 80 volumes from eff and non-eff ones, and calcuate the outcome of 150 experiments
      ['data/in006_t0120/subjects/r2_079_017_230101_080000', ...
      'data/in006_t0120/subjects/noneffdst_230101_080002'], ...
      80, 150, [], 'subj_range' = [5: 5: 80]);                         # r2 is not needed to be specified as the pewvious generated folder will be used to be randomly sampled
                                                                       # in the case the effected voxels were generated with 'grtdst' arguments, it must be selected the r2 range
    compute_mean_res(outfolder)                                        # average the individual results to get the power of the voxelwise correction
    
    clusterbased correction:
    generate_volumes(120, 6, 4^3, 120, [0.079, 0.017]);                # generate effected volumes (4-voxel-width cube of effected voxels for 120 subjects - mean 0.079, std 0.017)
    outfolder = subj_subsampling('clusterbased', ...
      ['data/in006_t0120/subjects/r2_079_017_230102_090000', ...
      'data/in006_t0120/subjects/noneffdst_230101_080002'], ...
      80, 150, [], 'subj_range' = [5: 5: 80], 
       'nb_noneffvx', get_nb_null_vx(4, 20));                          # r2 similarly treated as above
                                                                       # the non-effected volumes generated early will be used to pad the effected cube with 20 voxels each side
    compute_mean_res(outfolder)
    
## Folder structure <a name="foldstruct"/> (OUTDATED)

    ./
    ├── Napucco.m
    ├── utils/
    │   └── ...
    ├── generators/
    │   ├── create_volumes.m
    |   ├── resume_generation.m
    │   └── ...
    ├── core_functions/
    │   ├── subj_sampling.m
    │   ├── vc_correction.m
    │   └── ...
    ├── data/
    │   ├── inYYY_tYYYY/                                       # This is the home folder for all the simulation sharing the same number of regressors and number of timepoints
    │   │   ├── results/                                       # experiments
    │   │   │   ├── r2_muR2_stdR2_{K,V,B}DATE_TIME/              # experiments results for given R2 target sampled from population, cluster or voxel wise correction
    │   │   │   │   ├── hyparams.mat
    │   │   │   │   ├── opts.mat
    │   │   │   │   └── {voxelwise,clusterbased}_correction/                              # this folder contains all the individual outcome (pass or fail) of all the experiment repetition {results_XXX.mat}
    │   │   │   │       ├── grouping/                              # this folder contains all the individual outcome (pass or fail) of all the experiment repetition {results_XXX.mat}
    │   │   │   │       └── groupingavg/                           # this folder contains power.mat file containing the fpr or power results for the set of experiments
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


