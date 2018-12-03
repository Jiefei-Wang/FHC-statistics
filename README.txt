Title: A novel exact equal local levels higher criticism test
Paper authors: Jiefei Wang, Jeffrey C. Miecznikowski
Code author: Jiefei Wang (jwang96@buffalo.edu)
Software:
R:
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

MATLAB:
9.0.0.341360 (R2016a)
February 11, 2016


Description:
since the HC and FHC statistic requires a symbolic computation to obtain their p-values, MATLAB is used to compute them and R is used to run the simulations. In case where critical values and p-values are required in R, we first run the code in MATLAB and save the result on disk and then load them from R.

For the MATLAB code, each folder has its own README file

For the R code:
simulation1:
Compare the exact power of HC and FHC statistic under a beta-uniform mixture distribution(see section 4.1). The power of KS and asymptotic FHC are also obtained by the simulation

simulateion2:
Compare the power of Exact/permute HC/FHC statictic for a case-control simulation study(see section 4.2). The sample are simulated from a normal distribution, the correlation structure are described in the paper. The t-test is used to calculate the p-value for each variable and the HC and FHC statistic are computed based on the p-values obtained from the t-test. 


case_study:
Apply FHC,HC,KS to the data in Yang et al.(2012)(see section 5)

