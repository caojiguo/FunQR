This folder consists of the codes for the application of our algorithm on the Berkeley growth data in the manuscript.

You will need the following R files to implement the algorithm. I will first explain each R file and then outline how to run them in order.

###########################################################

0. nodes.RDS: contains all the nodes of the triangulation used in the data analysis.

   tri.RDS: contains the information about all the triangles formulated based on the nodes information saved in "nodes.RDS".

1. rd-FPCA-growth.R: the R code for data cleaning and applying FPCA on the dataset.

2. rd-Fitting.R: the R code used to generate Bernstein polynomials over a given triangulation and formulate the optimization problem corresponding to the functional quantile regression into a format that can be solved by a MATLAB/2018b package called "quadprog". The required inputs of the package "quadprog" are calculated and saved as RDS files.

3. cv_Fitting.R: the R code used to perform the cross validation for tuning parameters, which is similar to "rd-Fitting".

4. cv-R-MATLAB.R: the R code used to generate multiple MATLAB scripts that need to be submitted to MATLAB/2018b for the first round cross validation. Note that for the first round cross validation, there is only penalty term included into the model and we aim to select the "optimal" value for the corresponding tuning parameter from a given set of candidate values.

5. pre_compare.R: the R code used to summarize the results of the first round cross validation obtained from "cv-R-MATLAB.R". The output is the "optimal" value for the tuning parameter of the first penalty.

6. cv-R-MATLAB_2.R: the R code used to generate MATLAB scripts for the second round cross validation, which is similar to "cv-R-MATLAB.R".

7. cv_compare.R: the R code used to summarize the results of the second round cross-validation obtained from "cv-R-MATLAB_2.R". The output is the "optimal" value for the tuning parameter of the second penalty.

8. growth-plots.R: the R code used to generate plots for the estimated slope function obtained from the proposed method in the manuscript based on the Berkeley growth data. Specifically, it can generate a plot of the estimated slope function for five different quantiles. It also can create a heat map for the estimated slope function as the overall visualization.

9. growth-plots_2.R: the R code used to generate plots for the comparison of the estimated slope functions obtained from the proposed method and the conventional method to illustrate the advantages of our method. Specifically, it can generate two pairs of the fitted quantile functions obtained from the proposed method and the conventional method. It also can create a heat map for the estimated slope function obtained from the conventional method.

##########################################################

To reproduce the analysis, please follow the steps below.

We perform the following steps on the high performance computing clusters due to two reasons:

  A. The computational cost is very high because of the large scale of matrices in the optimization.

  B. The numerical errors of quadratic programming in MATLAB/2018b on a laptop (2018 MacBook Pro) can significantly impact the accuracy of the estimation. We compared the estimation results obtained from a 2018 MacBook Pro with the results obtained from the computing clusters based on some simulated data, and we observe that the estimation obtained from a laptop is much less accurate than that from the high performance computing clusters.

Because of the issues explained above, we also include the fitted results obtained by ourselves in a sub-folder named "results_of_model_fitting". If the user only wants to generate all the figures we used in the manuscript, please do the following steps:

a. In R, set the sub-folder "results_of_model_fitting" as the working directory, and issue
   
   ```
   source('growth-plots.R')

   source('growth-plots_2.R')

   ```


If the user wants to reproduce all the model fitting results and then generate the figures, please do the following steps:

Step 1: 

a. In R, set the current folder as the working directory, and issue

   ```
   source('rd-FPCA-growth.R')

   ```
   Then issue

   ```
   source('cv-R-MATLAB.R')

   ```

   This code creates 10 subfolders named "fold_1", "fold_2", ..., "fold_10", and each of them contains five MATLAB files named "sim1.m" , "sim2.m", ..., "sim5.m".

   The number of subfolders depends on the number of folds for the cross validation. In the manuscript, we use 10-fold cross validation and therefore, there will be 10 subfolders. The number of MATLAB files in each of the 10 subfolders depends on the number of candidate values for the tuning of the first penalty. In my code, I choose five candidate values.

b. Execute "sim1.m" , ..., "sim5.m" in MATLAB/2018b for each of the 10 subfolders "second_fold_1", "second_fold_2", ..., "second_fold_10" and save the outputs of "sim1.m", ..., "sim5.m" as txt files named "output1.txt", ...,"output5.txt" respectively.

   Then there should be five txt files named "output1.txt", ..., "output5.txt" in each of the 10 subfolders "fold_1", "fold_2", ..., "fold_10".

Step 2:

a. In R, issue
  
   ```
   source('cv-R-MATLAB_2.R')

   ```

   This code creates 10 subfolders named "second_fold_1", "second_fold_2", ..., "second_fold_10", and each of them contains five MATLAB files named "sim1.m" , ..., "sim5.m". The number of MATLAB files in each of the 10 subfolders depends on the number of candidate values for the tuning parameter of the second penalty. In my code, I choose five candidate values for it.

b. Execute "sim1.m" , ..., "sim5.m" in MATLAB/2018b for each of the 10 subfolders "second_fold_1", "second_fold_2", ..., "second_fold_10" and save the outputs of "sim1.m", ..., "sim5.m" as txt files named "output1.txt", ...,"output5.txt" respectively.

   Then there should be five txt files named "output1.txt", ..., "output5.txt" in each of the 10 subfolders "second_fold_1", "second_fold_2", ..., "second_fold_10".

Step 3: 

a. In R, issue

   ```
   source('cv_compare.R')
 
   ```
   This step generates a file named "sim.m". 

b. Execute "sim.m" in MATLAB/2018b and save the output as "output.txt". Then in R, issue

   ```
   source('growth-plots.R')

   source('growth-plots_2.R')

   ```
   to generate the plots presented in the manuscript.
