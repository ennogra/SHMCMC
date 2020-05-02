# SHMCMC
 Subspace Hamiltonian Markov Chain Monte Carlo

Matlab codes accompaning the paper: An Efficient Inverse Method for X-ray and Neutron Data with Subspace Hamiltonian Markov Chain Monte Carlo

Installation instructions:
1. Install Automatic Differentiation toolbox for Matlab (ADiMat), available at https://www.sc.informatik.tu-darmstadt.de/res/sw/adimat/
2. Download multiESS.m, available at https://github.com/lacerbi/multiESS
3. Add path ~/shmcmc (HMCMC functions) to Matlab search path

Example script folders
1. Demo_hmcmc_mcmc: compare the performance of Hamiltonian MCMC and Metropolis-MCMC (Figure 2 in the paper)
2. SiOx_SAXS_SHMCMC: SAXS data (for Figure 3)
3. SiOx_SAXS_HMCMC: SAXS data (for Figure S1)
4. XWFH_SHMCMC: XWFH data (for Figure 4)
5. XWFH_HMCMC: XWFH data (for Figure 5)
6. PSS_NRef_SHMCMC: neutron reflectivity (for Figure 6)
7. PSS_NRef_HMCMC: neutron reflectivity (for Figure S2)

Animations for the comparison of HMCMC and SHMCMC
1. XWFH_movie
