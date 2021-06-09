This folder consists of the codes for the algorithm used to perform the simulation studies in the manuscript.

You will need the following R files to run the simulations. In this README file, I will explain each R file and then outline how to run them in the correct order.

###########################################################

0. nodes.RDS: contains all the nodes of the triangulation used in the simulations.
   tri.RDS: contains the information about all the triangles formulated based on the nodes information saved in "nodes.RDS".

1. DataGen_simple.R: the R code used to generate data for the simple scenarios of the simulation study presented in the manuscript.

   DataGen_complex.R: the R code used to generate data for the complex scenarios of the simulation study presented in the manuscript.

2. Estimate_FPCs.R: the R code to apply FPCA on the dataset.

3. Fitting.R: Bernstein polynomials over a given triangulation are first generated and then the optimization problem corresponding to the functional quantile regression is formulated into a format that can be solved by a MATLAB package called "quadprog". The required inputs of the MATLAB package "quadprog" are calculated and saved as .RDS files.

4. cv_Fitting.R: the R code used to perform the cross validation for tuning parameters, which is similar to "Fitting.R".

5. cv-R-MATLAB.R: the R code used to generate multiple MATLAB scripts that need to be submitted to MATLAB/2018b for the first round cross validation. Note that for the first round cross validation, there is only penalty term included into the model and we aim to select the "optimal" value for the corresponding tuning parameter from a given set of candidate values.

6. pre_compare.R: the R code used to summarize the results of the first round cross validation obtained from "cv-R-MATLAB.R". The output is the "optimal" value for the tuning parameter of the first penalty.

7. cv-R-MATLAB_2.R: the R code used to generate MATLAB scripts for the second round cross validation, which is similar to "cv-R-MATLAB.R".

8. cv_compare.R: the R code used to summarize the results of the second round cross-validation obtained from "cv-R-MATLAB_2.R". The output is the "optimal" value for the tuning parameter of the second penalty.

9. results_check.R: the R code used to summarize the simulation results and compare the performance of the proposed method and the conventional method.

##########################################################

To reproduce the analysis, please follow the steps below.

Note that the following steps only present how to obtain one replicate of the simulation (simulate one time).

Regarding the simulation results in the manuscript, we repeat the following steps for 100 times under each simulation setting by using high performance computing clusters.

We recommend to perform the following steps on high performance computing clusters due to three reasons:

  A. The computational cost for each replication on a standard 2018 MacBook Pro is more than one hour.

  B. The numerical errors of quadratic programming in MATLAB/2018b on a laptop (2018 MacBook Pro) can significantly impact the accuracy of the estimation.
     We compared the estimation results obtained from a 2018 MacBook Pro with the results obtained from the computing clusters based on some simulated data, and we observe that the estimation obtained from a laptop is much less accurate than that from the high performance computing clusters.

  C. On the clusters, we can use parallel computing to perform the simulations.

We treat each simulation as a single computation job and in each step, and submit 100 computation jobs to the clusters. More specifically, we first create 100 folders for 100 replicates. Then we submit 100 jobs to the clusters to perform step 1 for each folder. Next, we submit another 100 jobs to the clusters to perform step 2 for each folder. Last, we submit 100 jobs to perform step 3 to obtain 100 simulation results.

Step 1: 

a. In R, set the current folder as the working directory, and issue

   ```
   source('DataGen_simple.R')

   ```
or

   ```
   source('DataGen_complex.R')

   ```

   Then issue

   ```
   source('Estimate_FPCs.R')

   ```

b. Next issue
 
   ```
   source('cv-R-MATLAB.R')

   ```

   This code creates 10 subfolders named "fold_1", "fold_2", ..., "fold_10", and each of them contains five MATLAB files named "sim1.m" , "sim2.m", ..., "sim5.m".

   The number of subfolders depends on the number of folds for the cross validation. In the manuscript, we use 10-fold cross validation and therefore, there will be 10 subfolders. The number of MATLAB files in each of the 10 subfolders depends on the number of candidate values for the tuning of the first penalty. In my code, I choose five candidate values.

c. Execute "sim1.m" , "sim2.m", ..., "sim5.m" in MATLAB/2018b for all the subfolders "fold_1", "fold_2", ..., "fold_10" and save the outputs of "sim1.m" , "sim2.m", ..., "sim5.m" as txt files named "output1.txt", ..., "output5.txt" respectively.

   Then there should be five txt files named "output1.txt", ..., "output5.txt" (because I used five canadidate values) in each of the 10 subfolders "fold_1", "fold_2", ..., "fold_10".

Step 2:

a. In R, issue
  
   ```
   source('cv-R-MATLAB_2.R')
 
   ```
   This code creates 10 subfolders named "second_fold_1", "second_fold_2", ..., "second_fold_10", and each of them contains three MATLAB files named "sim1.m" , "sim2.m" and "sim3.m". The number of MATLAB files in each of the 10 subfolders depends on the number of candidate values for the tuning parameter of the second penalty. In my code, I choose three candidate values for it.

b. Execute "sim1.m" , "sim2.m", and "sim3.m" in MATLAB/2018b for each of the 10 subfolders "second_fold_1", "second_fold_2", ..., "second_fold_10" and save the outputs of "sim1.m" , "sim2.m" and "sim3.m" as txt files named "output1.txt", "output2.txt" and "output3.txt" respectively.

   Then there should be three txt files named "output1.txt", "output2.txt" and "output3.txt" in each of the 10 subfolders "second_fold_1", "second_fold_2", ..., "second_fold_10".

Step 3: 

a. In R, issue

   ```
   source('cv_compare.R')
 
   ```
   This step generates a file named "sim.m". 

b. Execute "sim.m" in MATLAB/2018b and save the output as "output.txt".

   Suppose we have 100 folders named "simulation1", "simulation2", ..., "simulation100", and each folder contains one replicate of the simulation. Then issue the following code in R to get the summary of 100 replicates of the simulation.

   ```
   source('results_check.R')
 
   ```
