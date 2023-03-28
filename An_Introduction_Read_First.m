Introduction:
To run the code you need to install:
CasADi: https://web.casadi.org/ (please find the version fit your computer and MATLAB)
MPT: https://www.mpt3.org/Main/Installation (you can choose either Automatic installation or Manual installation, when you have installed MPT you will automatic install YALMIP)
SDPT3-4.0: https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/
(Note: installing SDPT3-4.0 and configurating it is a bit comprehensive, you may need a few Mex files.
Then you need to open MATLAB in the directory SDPT3-4.0, then in Matlab command window, type:
>> Installmex(1)
to activate the SDP solver. But it is not necessary to run the code, this is only needed to compute Ps in line 161 in 
the file offline_parameters_computation.m. You can avoid running this file and only use the provided parameters, such
that you do not need to install the solver.)


offline_parameters_computation.m:
This is the first file you need to run to get all parameters computed offline, it may take some time, 
and the results will be automatically saved in the file 'parameters.mat'. The parameters have been provided so if 
you do not want to update the values you can avoid running this file.

% You can check the functions in the folder 'Functions', but you do not need
% to run the following functions
compute_mrpi_set.m
A function to compute the RPI set, implemented according to:
Rakovic, Sasa V., et al. "Invariant approximations of the minimal robust positively invariant set." IEEE Transactions on automatic control 50.3 (2005): 406-410.

ComputeFeasibleRegion.m
A function to compute the feasible region

InitialSetComputation.m
A function to initiate the uncertainty quantification of the set W_true based on W, using initial information set
\mathcal{I}_0^w

ModelingCar.m
A function for modeling the behaviour of EV and LV

NominalRobustMPC.m
A function for the conventional Robust MPC controller

UQRobustMPC.m
A function for the uncertainty quantification Robust MPC controller

% You can run the files in the Cases folder to investigate
Case_1.m:
The file to get Fig. 1, results are saved in Results/Results_1.mat, and figures are produced by Figures/Fig_Case_1.m

Case_2.m:
The file to get Fig. 2, results are saved in Results/Results_2.mat, and figures are produced by Figures/Fig_Case_2.m

Case_3.m:
The file to get Fig. 3-4, results are saved in Results/Results_3_large.mat (i.e., |I_0^w| = 20000) and Results/Results_3_small.mat (i.e., |I_0^w| = 100), and figures are produced by Fig_Case_3.m

Case_4.m:
The file to get Fig. 5, results are saved in Results/Results_4.mat, and figures are produced by Figures/Fig_Case_4.m

Case_5.m:
The file to get Table, results are saved in Results/Results_5.mat, e.g., Results_5_10.mat indicates the results with |I_0^w| = 10, and so on



