# Learning-based Rigid Tube Model Predictive Control
This repository contains code for the article
```
@inproceedings{gao2024learning,
  title={Learning-based Rigid Tube Model Predictive Control},
  author={Yulong Gao, Shuhao Yan, Jian Zhou, Mark Cannon, Alessandro Abate, and Karl H. Johansson},
  booktitle={Learning for Dynamics and Control Conference},
  year={2023},
  pages={}
} 
```
which has been accepted for presentation and publication at the 6th Annual Learning for Dynamics & Control Conference (L4DC), 2024.
## Packages for running the code
To run the code you need to install:

**CasADi**: https://web.casadi.org/;

**MPT**: https://www.mpt3.org/Main/Installation; (The installation will automatically install Yalmip, which is also necessary for running the code.)

## Introduction to the files

## offline_parameters_computation.m
Calculate and define all parameters and save the results in **parameters.mat**. Since the parameters have been provided in the repository, if you do not want to update those values you can avoid running this file.

## Functions
### MRPISet.m:

Compute the minimum robust positive invariant set $\mathbb{S}$.

### ComputeFeasibleRegion.m:

Compute the feasible region corresponding to different disturbance sets, e.g., $\mathbb{W}$, $\mathbb{W}_{\rm true}$, and $\hat{\mathbb{W}}_k^{\star}$.

### InitialSetComputation.m:

Learn the initial uncertainty set $\hat{\mathbb{W}}^{\star}_0$ based on $\mathbb{W}$, using initial information set
$\mathcal{I}_0^w$.

### ModelingCar.m:

Modeling of the ego vehicle (EV) and the leading vehicle (LV).

### NominalRobustMPC.m:

Conventional Robust MPC controller.

### UQRobustMPC.m:

The proposed uncertainty quantification-based Robust MPC controller.

### ScenarioMPC.m:
The scenario MPC controller.

## Cases

### Case_1_Feasible_Region.m
Compute the feasible region of UQ-RMPC and nominal RMPC, results are saved in **Results/Results_1.mat**, and figures are produced by **Results/Fig_Case_1.m**.

### Case_2_MC_Different_Initial_InformationSet.m
Monte-Carlo simulation learning the set $\mathbb{W}_{\rm true}$ with different sizes of initial information set $\mathcal{I}_0^w$, results are saved in **Results/Results_2.mat**, and figures are produced by **Results/Fig_Case_2.m**.

### Case_3_Online_UQRMPC_Different_Initial_InformationSet.m
Online evaluation of UQ-RMPC with different initial information set $\mathcal{I}_0^w$, results are saved in **Results/Results_3_large.mat** ($|\mathcal{I}_0^w| = 20000$) and **Results/Results_3_small.mat** ($|\mathcal{I}_0^w| = 100$), and figures are produced by **Results/Fig_Case_3.m**.

### Case_4_Online_UQRMPC_Long_Simulation_Step.m
Monte-Carlo simulation of UQ-RMPC when the simulation time is long enough, results are saved in **Results/Results_4.mat**, and figures are produced by **Results/Fig_Case_4.m**.

### Case_5_Feasibility_Evaluation_UQRPC.m
Monte-Carlo simulation to evaluate the feasibility of UQ-RMPC with different initial information set $\mathcal{I}_0^w$, results are saved in **Results/Results_5.mat**, e.g., **Results_5_10.mat** indicates the results with $|\mathcal{I}_0^w| = 10$, and so on.

### Case_6_Compare_With_SCMPC.m
Comparing the computation time with Scenario MPC.

## Some implementation details
(1) It is not necessary to run **offline_parameters_computation.m** if you do not want to update the parameter values, but you need to run **run_first.m** to add the path of folders.

(2) In the article, the horizon $\nu_k$ is updated according to Algorithm 1, but this will change the number of constraints of MPC. In our code, we implemented Algorithm 1 and found that $\nu_k$ is almost equal to $\nu_s$. Therefore, we use $\nu_s$ to replace $\nu_k$. In practical applications, we can make $\nu_k$ long enough such that the condition on $\nu_k$ will be satisfied. For example, a suggestion can be $\nu_k = 2\nu_s$.
