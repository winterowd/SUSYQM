SUSYQM
======
This code uses the Hybrid Monte Carlo algorithm to study supersymmetry in Quantum Mechanics. Correlations functions
of the primary field \Phi, the auxiliary field F, and the fermion bilinears are computed "on the fly". For more
details see arXiv:1302.2361 and arXiv:1311.3487.

Compilation requires GSL for the random number generator and for computing gaussian distributed random numbers.
Make target is SUSYQM_hmc i.e.

>> make SUSYQM_hmc

All simulation parameters are taken from standard input
Usage: ./SUSYQM_hmc N lambdaR warms trajecs meas step_size num_steps seed
N: (integer) size of lattice
lambdaR: (float) renormalized coupling (referred to as g in arXiv:1311.3487)
warms: (integer) number of warmup trajectories to perform
trajecs: (integer) number of actual trajectories to perform
meas: (integer) number of trajectories before a measurement occurs
step_size: (float) step size for integrating EOM
num_steps: (float) number of updating steps per Molecular Dynamics trajectory (code uses leap frog integration scheme)
seed: (integer) seed for random number generator

Code saves correlation functions to text files. For example, for <F>, it's values are sent to the file 'vevf1.out'.
