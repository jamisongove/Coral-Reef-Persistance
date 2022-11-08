# Coral-Reef-Persistence
Code and data for Gove and Williams et al. (in review)

This is a repository for code and data that are used within Gove and Williams et al. (in review)

Software: 
GAMM_Disturbance.R requires R

All other code requires matlab.

Contents: 

PreDisturbance.m: Runs all analysis on CORAL.PERSISTENCE.RAW.csv that produces Figure 2a,b,c in our submitted manuscript. The output from this script that produces Figure 2c is FIG2c_PERCENT.DIFFERENCES.xlsx

Disturbance.m: Runs all analysis on CORAL.PERSISTENCE.RAW.csv that produces the dataset GAMM_MHW.csv. GAMM_MHW.csv is the input for the R script (GAMM_Disturbance.R) that produces the results in Fig. 3d

PostDisturbance.m: Runs all anayis on CORAL.PERSISTENCE.RAW.csv that produces the data set OLR_REEF.BUILDER.COVER.csv, which is the input for OLR_Probability_Scenario.m

OLR_Probability_Scenario.m: Takes OLR_REEF.BUILDER.COVER.csv and runs the following: RunOrdinalRegression.m and deltaProbability_2variables. The outputs are the probability matrices for acheiveing HIGH, MODERATE, and LOW reef-builder cover. This also produces figure 4. 


plot_corr_marix.m is called within Predisturbance.m, Disturbance.m, and PostDisturbance.m 

