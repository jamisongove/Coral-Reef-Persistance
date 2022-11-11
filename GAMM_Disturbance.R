rm(list=ls())# cleans the workspace
gc() # cleans the workspace

library(gamm4) # loads the gamm4 package (used for running gamm)
library(foreach) # loads the foreach package (used for parrallelized looping)
library(doMC)   
library(MuMIn)
library(ggplot2) 
library(gmt)
library(tidyverse)
registerDoMC(12) # WARNING!! defines the number of cores to be used. To know how many cores you can use, use detectCores().


### Here, set the path of your working directory. Rename this to your directory
setwd("~/Documents/...")
### Get the dataset
data=read.csv('GAMM_MHW.csv',header = T)
datalist=read.csv('GAMM_MHW.csv',header = T)
### Create the matrix of Euclidean distances between locations
distance_island_matrix=matrix(NA,length(data$Lat),length(data$Lat))   # create the empty matrix in which the distances between reefs are going to be recorded
for(isl1 in 1:length(data$Lat)){    
  for(isl2 in isl1:length(data$Lat)){
    distance_island_matrix[isl1,isl2]=geodist(data$Lat[isl1],data$Long[isl1], data$Lat[isl2],data$Long[isl2], units="km")   # calculate the distance between 2 reefs
  }}
distance_island_matrix=as.dist(t(distance_island_matrix))   # create a 'dist' type object for hclust to process
full=hclust(distance_island_matrix,method="complete")   # performs a hierarchical cluster analysis using a set of dissimilarities for the n objects being clustered.

###############
nameRV='RelChange' # defines response variable
h_cut=10      # defines the distance (in km) to form the random factor in the GAMM based on the hierarchical cluster analysis
###############
GROUP=cutree(full, h=h_cut) # cuts the hierarchical clustering analysis tree
GROUP  <- as.factor(GROUP)
datalist=cbind(data,GROUP) # column bind the dataset and the previsouly determined grouping


#################################### Distance v Group Comparision ############################################################################################################
require(foreach)
distance_exploration=seq(1,50,length.out=1000)
results=as.data.frame(foreach(h_cut=distance_exploration,.combine=rbind)%do%{    # defines the distance (in km) to form the random factor in the GAMM based on the hierarchical cluster analysis
  ###############
  c(h_cut,length(unique(cutree(full, h=h_cut))))
})
colnames(results)=c("Distance (km)", "Number of groups")
plot(results,las=1,type='l',col='red',main='Distance v Groups')
#############################################################################################################################################################

########################################################################################
################################## GAMM models #########################################
########################################################################################
base_k=4

# Define the response variables
colnames(datalist)[1]

# All models for a given group
name='RelChange'

mod<-uGamm(RelChange ~ s(Depth, k = base_k, bs = "cr")  # make sure this list of predictors matches EXACTLY the list in the data set
           + s(GearRegRank, k = base_k, bs = "cr")
           + s(total_biomass, k = base_k, bs = "cr")
           + s(grazers, k = base_k, bs = "cr")
           + s(scrapers,  k = base_k, bs = "cr")
           + s(DHW_MAX_2015, k = base_k, bs = "cr")
           + s(OSDS_Eff, k = base_k, bs = "cr")
           + s(Nuts, k = base_k, bs = "cr")
           + s(Imperv, k = base_k, bs = "cr")
           + s(Sediment, k = base_k, bs = "cr")
           + s(WPow975pct, k = base_k, bs = "cr")
           + s(CHL_MHW, k = base_k, bs = "cr")
           + s(RainMax3dS, k = base_k, bs = "cr"),
           data = datalist, random=~(1|GROUP), family="gaussian", lme4 = T)
 
# Run all model combinations and output model evaluation metrics
# m.lim defines the number of predictors that can be included in any single candidate model (here set to 5 to reduce overfitting)
M.set.ori<-dredge(mod,beta=FALSE, rank="AICc", extra = alist(AIC, "R^2", "adjR^2"), m.lim=c(1,5))
M.set=as.data.frame(M.set.ori)

###############################3#################################################################################

# Filter models by AICc delta <2
subtab=M.set[which(M.set$delta<=2),]
subtab

############################ CUMULATIVE AKAIKE WEIGHT TABLE ############################

###GET TABLE OF TOP (AICc Delta <=2) MODEL RESULTS. 
stab <- subtab %>% 
  mutate(across(contains('s('), ~ ifelse(.x=='+' & !is.na(.x), T, F))) %>% 
  rename_with(~sub('s\\(([A-z]*.*), k.*','\\1',.x), contains('s('))

stab$vars <- apply(stab, 1, function(x) paste(names(x[x==T]), collapse='; '))

# print compact
mod_tab <- stab %>% select(-`(Intercept)`, -where(~ is.logical(.x)))

varlist <- unique(unlist(lapply(mod_tab$vars, strsplit,'; ')))

sapply(varlist, function(v) {
  
  this.weight <- sum(mod_tab$weight[grepl(v, mod_tab$vars)])
  
})
View(mod_tab)
########################## Average across top-ranking models ##########################

avg=model.avg(get.models(M.set.ori,subset = delta<=2)) 

############################ FITTED FUNCTIONS ############################

plot_pred <- function(model){
  
  datamedian=(as.data.frame(lapply(apply(datalist[,4:16],2, median, na.rm = TRUE), rep, nrow(datalist)))) # DOUBLE CHECK THIS NUMBER - Be sure to define predictor variables correctly!
  
  varnames <- strsplit(mod_tab$vars[rownames(mod_tab)==model], '; ')[[1]]
  
  plot_frame <- lapply(varnames, function(name){
    #####
    datamedian[name]=datalist[name]
    #####
    
    M<-predict(avg, datamedian,se.fit=T)  # Note if only one top-ranking model then replace "avg" with "modcr"
    fit.vals<-M$fit
    se.vals<-M$se
    
    lcl<-fit.vals-(1.645*se.vals)
    ucl<-fit.vals+(1.645*se.vals)
    
    tmp<-data.frame(fit.vals, lcl, ucl,datalist[name],name)
    names(tmp)<-c('RelChange', "lcl", "ucl",'x','variable')
    tmp
  }) %>% bind_rows()
  
  long_data <- data %>% pivot_longer(cols=4:(ncol(datalist)-1), names_to = 'variable', values_to = 'xx') %>% 
    filter(variable %in% varnames) %>%
    rename(Change = RelChange)
  
  plot_frame %>% inner_join(long_data) %>%
    group_by(variable) %>% 
    ### Plot a defined fitted function using ggplot 
    ggplot()+
    geom_line(aes(x = x, y = RelChange), size =2, colour='red') +
    geom_ribbon(aes(x = x, y = RelChange, ymin=lcl,ymax=ucl), alpha=0.3, fill = 'red')+
    geom_point(aes(x = xx, y = Change), color="blue") +
    xlab('')+
    ylab('Coral Change(% difference)')+
    coord_cartesian(ylim=c(-95, 15)) +
    theme_bw(base_size=8)+
    facet_wrap(~variable, scales='free_x') + 
    theme(legend.position="right")
  ############################
}

plot_pred(model ='3618') #define the model you want to visualize - in this case the top model. 
########################################################################################
########################################################################################
########################################################################################



