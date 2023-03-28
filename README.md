# Tube-Robust-MPC-with-Uncertainty-Quantification
This is the MATLAB code for tube robust MPC with uncertainty quantification.
## Packages for running the code
To run the code you need to install:

**CasADi**: https://web.casadi.org/ (please find the version that fits your computer and MATLAB);

**MPT**: https://www.mpt3.org/Main/Installation (you can choose either Automatic installation or Manual installation; when you have installed MPT you will automatically install YALMIP);

**SDPT3-4.0**: https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/
(Installation of SDPT3-4.0 and configuration can be a bit comprehensive, you may need a few Mex files. Then you need to open MATLAB in the directory SDPT3-4.0, then in MATLAB command window, type:

*Installmex(1)*

to activate the SDP solver. But it is not necessary to run the code, this is only needed to compute $P_s$ in line 161 in the file **offline_parameters_computation.m**. You can avoid running this file and only use the provided parameters, such that you do not need to install the solver.)

## The file *offline_parameters_computation.m*
This is the first file you need to run to get all parameters computed offline, it may take some time, and the results will be automatically saved in **parameters.mat**. The parameters have been provided in the code so if 
you do not want to update those values you can avoid running this file.

## The *Functions* folder
**compute_mrpi_set.m**:

A function to compute the RPI set, implemented according to:
Rakovic, Sasa V., et al. "Invariant approximations of the minimal robust positively invariant set." *IEEE Transactions on automatic control* 50.3 (2005): 406-410.

**ComputeFeasibleRegion.m**

A function to compute the feasible region corresponding to different disturbance sets, e.g., $\mathbb{W}$, $\hat{\mathbb{W}}_k^*$, ${\mathbb{W}}_{\rm true}$, and so on.

**InitialSetComputation.m**

A function to initiate the uncertainty quantification of the set $\mathbb{W}_{\rm true}$ based on $\mathbb{W}$, using initial information set
$\mathcal{I}_0^w$.

**ModelingCar.m**

A function for modeling the EV and LV

**NominalRobustMPC.m**

A function for the conventional Robust MPC controller

**UQRobustMPC.m**

A function for the uncertainty quantification Robust MPC controller

## The **Cases** folder

**Case_1.m**

The file to get Fig. 1, corresponding results are saved in *Results/Results_1.mat*, and figures are produced by *Figures/Fig_Case_1.m*.

**Case_2.m**

The file to get Fig. 2, corresponding results are saved in *Results/Results_2.mat*, and figures are produced by *Figures/Fig_Case_2.m*.

**Case_3.m**

The file to get Fig. 3-4, corresponding results are saved in *Results/Results_3_large.mat* (i.e., $|\mathcal{I}_0^w| = 20000$) and *Results/Results_3_small.mat* (i.e., $|\mathcal{I}_0^w| = 100$), and figures are produced by *Figures/Fig_Case_3.m*.

**Case_4.m**

The file to get Fig. 5, results are saved in *Results/Results_4.mat*, and figures are produced by *Figures/Fig_Case_4.m*.

**Case_5.m**

The file to get Table I, results are saved in *Results/Results_5.mat*, e.g., *Results_5_10.mat* indicates the results with $|\mathcal{I}_0^w| = 10$, and so on.

## The **Results** folder

The folder saves the data by running the files in the folder **Cases** by the author. If you run the files in the folder **Cases**, the data will be saved in the main directory then you need to manually move it to the **Results** folder

## The **Figures** folder

This folder saves the files for making the figures used in the article, you can directly run them to reproduce the figures.

## Some implementation details
(1) It is not necessary to run **offline_parameters_computation.m** if you do not want to update the parameter values, but you need to run **run_first.m** to add the path of folders.

(2) In the article, the horizon $\nu_k$ in (20) is updated according to Algorithm, but this will change the number of constraints of MPC. In our code, we implemented Algorithm 1 and found that $\nu_k$ is almost equal to $\nu_s$. Therefore, we use $\nu_s$ to replace $\nu_k$ (see line 136 in the file *offline_parameters_computation.m*). In practical applications, we can make $\nu_k$ long enough such that the condition on $\nu_k$ such that the condition above Algorithm 1 will be satisfied. For example, a suggestion can be $\nu_k = 2\nu_s$.
