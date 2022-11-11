# Coral-Reef-Persistence
This repository is for the reviewers of Gove and Williams et al. (in review)

All data and code used within Gove and Williams et al. (in review) are below. 

Software: 
GAMM_Disturbance.R requires R

All other code requires matlab.

CONTENTS:

DATA FILES

CORAL.PERSISTENCE.RAW.csv: Benthic and fish survey data from 2003 - 2019 along with the associated land-sea human impacts and environmental factors for each reef for all years from 2000 - 2019. This data set is called by numerous scripts below to produce the results presented in Figures 2, 3, and 4. 

CORAL.REEF.TRAJETORY.csv: Positive and negative trajectory reefs along with the associated land-sea human impacts and environmental factors presented in Figures 2a,b, and used to calculate percent difference in Figure 2c. 

FIG2c.PERCENT.DIFFERENCE.xlsx: Percent differences used in Figure 2c

GAMM.MHW.csv: Coral change (percent difference) between 2014/2015 and 2016 along with the associated land-sea human impacts and environmental factors for all reefs used in the GAMM analysis presented in Fig. 3.  

OLR_REEF.BUILDER.COVER.csv: Total cover of reef-building organisms in 2019 along with the associated land-sea human impacts and environmental factors for all reefs used in the Ordinal Logistic Regression analysis presented in Fig 4. 

HUMAN.IMPACTS_ENVIRONMENTAL.FACTORS_100mGRID.csv: High reslution (100m) data for all land-sea human impacts and environmental factors used in Figure 1 and supplementary figures. Data are in annual time steps and primarly extend from 2000 - 2019. Some factors extend further back in time. All data in this file are from the 10 m contour, but these data are available from 0 - 20 m depth upon request. 

CODE

PreDisturbance.m: Runs all analysis on CORAL.PERSISTENCE.RAW.csv that produces Figure 2a,b,c in our submitted manuscript. The output from this script is the data file CORAL.REEF.TRAJECTORY that produces Figures 2a,b,c. See also FIG2c_PERCENT.DIFFERENCES.xlsx which is produced by this script and are the percent differences in Figure 2c. 

Disturbance.m: Runs all analysis on CORAL.PERSISTENCE.RAW.csv that produces the dataset GAMM_MHW.csv. GAMM_MHW.csv is the input for the R script (GAMM_Disturbance.R) that produces the results in Figure 3d and Supplementary Figure 2. 

GAMM_Disturbance.R: Runs the GAMM analysis used to produce results for Figure 3d and Figure S2. 

PostDisturbance.m: Runs all analysis on CORAL.PERSISTENCE.RAW.csv that produces the data set OLR_REEF.BUILDER.COVER.csv, which is the input for OLR_Probability_Scenario.m

OLR_Probability_Scenario.m: Takes OLR_REEF.BUILDER.COVER.csv and runs the following: RunOrdinalRegression.m and deltaProbability_2variables. The outputs are the probability matrices for achieving HIGH, MODERATE, and LOW reef-builder cover. This also produces figure 4. 

plot_corr_matrix.m is called within Predisturbance.m, Disturbance.m, and PostDisturbance.m and produced the correlation matrices show in the supplementary information. 

