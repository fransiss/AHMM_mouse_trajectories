# AHMM_mouse_trajectories


Python code to run Autoregressive Hidden Markov Model (AHMM) on longitudinal data of weight incorporating age, sex, and APOE genotype to identify emergent weight trajectories and phenotypes. 
The data of the hAPOE mouse colony used in the analyses is reported in the /data subfolder.

Pipeline: 
1- runAHMM_1.m to run 10 times the AHMM 
2 - analyzeTrajectories_2.m to select the model with max trace and analyze the trajectories based on the transition matrix
3 - results are saved in the Result folder
