# SHMCMC
 Subspace Hamiltonian Markov Chain Monte Carlo

SHMCMC accompanies the paper: Efficient X-ray and Neutron Data Analysis with Subspace Hamiltonian Markov Chain Monte Carlo

Installation instructions:
1. Install Automatic Differentiation toolbox for Matlab (ADiMat), available at https://www.sc.informatik.tu-darmstadt.de/res/sw/adimat/
2. Add path ~/shmcmc (HMCMC functions) to Matlab search path
3. Examples scripts

3a. Demo_hmcmc_mcmc: compare the performance of Hamiltonian MCMC and Metropolis-MCMC (Figure 1 in the paper)

3b. SiOx_SAXS_SHMCMC: SAXS data (for Figure 2)
  shmcmc_SAXS.m: script to run SHMCMC
  shmcmc_SAXS_result_analysis.m: script to analyze and plot the results.
3c. SiOx_SAXS_HMCMC: SAXS data (for Figure S1)
  hmcmc_SAXS.m: script to run HMCMC
  hmcmc_SAXS_result_analysis.m: script to analyze and plot the results.
3d. XWFH_SHMCMC: XWFH data (for Figure 3)
  shmcmc_XWFH.m: script to run SHMCMC
  shmcmc_XWFH_analysis_result.m: script to analyze and plot the results.
3e. XWFH_HMCMC: XWFH data (for Figure 4)
  hmcmc_XWFH.m: script to run HMCMC
  hmcmc_XWFH_analysis_result.m: script to analyze and plot the results.
3f. PSS_NRef_SHMCMC: neutron reflectivity (for Figure 5)
  shmcmc_NREF.m: script to run SHMCMC
  shmcmc_NREF_result_analysis.m: script to analyze and plot the results.
3e. PSS_NRef_HMCMC: neutron reflectivity (for Figure S2)
  shmcmc_NREF.m: script to run HMCMC
  shmcmc_NREF_result_analysis.m: script to analyze and plot the results.
