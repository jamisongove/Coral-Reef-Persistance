%% REEF TRAJECTORIES PRE-DISTURBANCE
%Hawaii Time Series Analysis
%Gove & Williams et al., Nature (in review)
%Updated 7-NOV-2022
clear all close all

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 395);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["SurveyID", "SiteID", "Name", "Day", "Month", "Year", "Depth_survey", "Lat", "Long", "Coral", "Macro", "Substrate", "Turf", "Calg", "Total", "MPA_25m_Name", "MPA_25m_Estab_Yr", "MPA_25m_GearRegSummary", "MPA_25m_GearRegRank", "GRID_ID_100m", "OSDS_Eff_2000", "OSDS_Eff_2001", "OSDS_Eff_2002", "OSDS_Eff_2003", "OSDS_Eff_2004", "OSDS_Eff_2005", "OSDS_Eff_2006", "OSDS_Eff_2007", "OSDS_Eff_2008", "OSDS_Eff_2009", "OSDS_Eff_2010", "OSDS_Eff_2011", "OSDS_Eff_2012", "OSDS_Eff_2013", "OSDS_Eff_2014", "OSDS_Eff_2015", "OSDS_Eff_2016", "OSDS_Eff_2017", "OSDS_Eff_2018", "OSDS_Eff_2019", "OSDS_N_2000", "OSDS_N_2001", "OSDS_N_2002", "OSDS_N_2003", "OSDS_N_2004", "OSDS_N_2005", "OSDS_N_2006", "OSDS_N_2007", "OSDS_N_2008", "OSDS_N_2009", "OSDS_N_2010", "OSDS_N_2011", "OSDS_N_2012", "OSDS_N_2013", "OSDS_N_2014", "OSDS_N_2015", "OSDS_N_2016", "OSDS_N_2017", "OSDS_N_2018", "OSDS_N_2019", "Golf_N_2000", "Golf_N_2001", "Golf_N_2002", "Golf_N_2003", "Golf_N_2004", "Golf_N_2005", "Golf_N_2006", "Golf_N_2007", "Golf_N_2008", "Golf_N_2009", "Golf_N_2010", "Golf_N_2011", "Golf_N_2012", "Golf_N_2013", "Golf_N_2014", "Golf_N_2015", "Golf_N_2016", "Golf_N_2017", "Golf_N_2018", "Golf_N_2019", "Imperv_1992", "Imperv_1993", "Imperv_1994", "Imperv_1995", "Imperv_1996", "Imperv_1997", "Imperv_1998", "Imperv_1999", "Imperv_2000", "Imperv_2001", "Imperv_2002", "Imperv_2003", "Imperv_2004", "Imperv_2005", "Imperv_2006", "Imperv_2007", "Imperv_2008", "Imperv_2009", "Imperv_2010", "Imperv_2011", "Imperv_2012", "Imperv_2013", "Imperv_2014", "Imperv_2015", "Imperv_2016", "Imperv_2017", "Imperv_2018", "Imperv_2019", "Sediment_1990", "Sediment_1991", "Sediment_1992", "Sediment_1993", "Sediment_1994", "Sediment_1995", "Sediment_1996", "Sediment_1997", "Sediment_1998", "Sediment_1999", "Sediment_2000", "Sediment_2001", "Sediment_2002", "Sediment_2003", "Sediment_2004", "Sediment_2005", "Sediment_2006", "Sediment_2007", "Sediment_2008", "Sediment_2009", "Sediment_2010", "Sediment_2011", "Sediment_2012", "Sediment_2013", "Sediment_2014", "Sediment_2015", "Sediment_2016", "Sediment_2017", "Sediment_2018", "Sediment_2019", "HPop15km_2000", "HPop15km_2001", "HPop15km_2002", "HPop15km_2003", "HPop15km_2004", "HPop15km_2005", "HPop15km_2006", "HPop15km_2007", "HPop15km_2008", "HPop15km_2009", "HPop15km_2010", "HPop15km_2011", "HPop15km_2012", "HPop15km_2013", "HPop15km_2014", "HPop15km_2015", "HPop15km_2016", "HPop15km_2017", "HPop15km_2018", "HPop15km_2019", "HPop15km_2020", "RainAnnSum_1990", "RainAnnSum_1991", "RainAnnSum_1992", "RainAnnSum_1993", "RainAnnSum_1994", "RainAnnSum_1995", "RainAnnSum_1996", "RainAnnSum_1997", "RainAnnSum_1998", "RainAnnSum_1999", "RainAnnSum_2000", "RainAnnSum_2001", "RainAnnSum_2002", "RainAnnSum_2003", "RainAnnSum_2004", "RainAnnSum_2005", "RainAnnSum_2006", "RainAnnSum_2007", "RainAnnSum_2008", "RainAnnSum_2009", "RainAnnSum_2010", "RainAnnSum_2011", "RainAnnSum_2012", "RainAnnSum_2013", "RainAnnSum_2014", "RainAnnSum_2015", "RainAnnSum_2016", "RainAnnSum_2017", "RainAnnSum_2018", "RainAnnSum_2019", "RainMax3dS_1990", "RainMax3dS_1991", "RainMax3dS_1992", "RainMax3dS_1993", "RainMax3dS_1994", "RainMax3dS_1995", "RainMax3dS_1996", "RainMax3dS_1997", "RainMax3dS_1998", "RainMax3dS_1999", "RainMax3dS_2000", "RainMax3dS_2001", "RainMax3dS_2002", "RainMax3dS_2003", "RainMax3dS_2004", "RainMax3dS_2005", "RainMax3dS_2006", "RainMax3dS_2007", "RainMax3dS_2008", "RainMax3dS_2009", "RainMax3dS_2010", "RainMax3dS_2011", "RainMax3dS_2012", "RainMax3dS_2013", "RainMax3dS_2014", "RainMax3dS_2015", "RainMax3dS_2016", "RainMax3dS_2017", "RainMax3dS_2018", "RainMax3dS_2019", "WPow975pct_1998", "WPow975pct_1999", "WPow975pct_2000", "WPow975pct_2001", "WPow975pct_2002", "WPow975pct_2003", "WPow975pct_2004", "WPow975pct_2005", "WPow975pct_2006", "WPow975pct_2007", "WPow975pct_2008", "WPow975pct_2009", "WPow975pct_2010", "WPow975pct_2011", "WPow975pct_2012", "WPow975pct_2013", "WPow975pct_2014", "WPow975pct_2015", "WPow975pct_2016", "WPow975pct_2017", "WPow975pct_2018", "WPow975pct_2019", "DHW_MAX_2014", "DHW_MAX_2015", "SST_MAX_2015", "SSTA_MAX_2015", "HS_MAX_2015", "DHW_MAX_2019", "SST_MEAN_2000", "SST_MEAN_2001", "SST_MEAN_2002", "SST_MEAN_2003", "SST_MEAN_2004", "SST_MEAN_2005", "SST_MEAN_2006", "SST_MEAN_2007", "SST_MEAN_2008", "SST_MEAN_2009", "SST_MEAN_2010", "SST_MEAN_2011", "SST_MEAN_2012", "SST_MEAN_2013", "SST_MEAN_2014", "SST_MEAN_2015", "SST_MEAN_2016", "SST_MEAN_2017", "SST_MEAN_2018", "SST_STD_2000", "SST_STD_2001", "SST_STD_2002", "SST_STD_2003", "SST_STD_2004", "SST_STD_2005", "SST_STD_2006", "SST_STD_2007", "SST_STD_2008", "SST_STD_2009", "SST_STD_2010", "SST_STD_2011", "SST_STD_2012", "SST_STD_2013", "SST_STD_2014", "SST_STD_2015", "SST_STD_2016", "SST_STD_2017", "SST_STD_2018", "CHL_MHW", "CHL_2016", "CHL_2017", "CHL_2018", "CHL_2019", "PAR_MHW", "PAR_2016", "PAR_2017", "PAR_2018", "PAR_2019", "total_biomass_2000", "total_biomass_2001", "total_biomass_2002", "total_biomass_2003", "total_biomass_2004", "total_biomass_2005", "total_biomass_2006", "total_biomass_2007", "total_biomass_2008", "total_biomass_2009", "total_biomass_2010", "total_biomass_2011", "total_biomass_2012", "total_biomass_2013", "total_biomass_2014", "total_biomass_2015", "total_biomass_2016", "total_biomass_2017", "total_biomass_2018", "total_biomass_2019", "herbivores_2000", "herbivores_2001", "herbivores_2002", "herbivores_2003", "herbivores_2004", "herbivores_2005", "herbivores_2006", "herbivores_2007", "herbivores_2008", "herbivores_2009", "herbivores_2010", "herbivores_2011", "herbivores_2012", "herbivores_2013", "herbivores_2014", "herbivores_2015", "herbivores_2016", "herbivores_2017", "herbivores_2018", "herbivores_2019", "browsers_2000", "browsers_2001", "browsers_2002", "browsers_2003", "browsers_2004", "browsers_2005", "browsers_2006", "browsers_2007", "browsers_2008", "browsers_2009", "browsers_2010", "browsers_2011", "browsers_2012", "browsers_2013", "browsers_2014", "browsers_2015", "browsers_2016", "browsers_2017", "browsers_2018", "browsers_2019", "grazers_2000", "grazers_2001", "grazers_2002", "grazers_2003", "grazers_2004", "grazers_2005", "grazers_2006", "grazers_2007", "grazers_2008", "grazers_2009", "grazers_2010", "grazers_2011", "grazers_2012", "grazers_2013", "grazers_2014", "grazers_2015", "grazers_2016", "grazers_2017", "grazers_2018", "grazers_2019", "scrapers_2000", "scrapers_2001", "scrapers_2002", "scrapers_2003", "scrapers_2004", "scrapers_2005", "scrapers_2006", "scrapers_2007", "scrapers_2008", "scrapers_2009", "scrapers_2010", "scrapers_2011", "scrapers_2012", "scrapers_2013", "scrapers_2014", "scrapers_2015", "scrapers_2016", "scrapers_2017", "scrapers_2018", "scrapers_2019"];
opts.VariableTypes = ["string", "string", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["SurveyID", "SiteID"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["SurveyID", "SiteID", "Name", "MPA_25m_Name", "MPA_25m_GearRegSummary"], "EmptyFieldRule", "auto");

% Import the data
drivers_raw = readtable("/Users/jgove/Documents/jamisongove/EOD/KONA_IEA/PaperI/data/Analysis/Final_Data/Revision/Submission/CORAL.PERSISTENCE.RAW.csv", opts);

%% DATA ANALYSIS %% 
%Renaming variable
drivers = drivers_raw;

drivers = drivers_raw; 
[A I] = sort(drivers.Year); drivers = drivers(I,:);
[A I] = sort(drivers.SiteID); drivers = drivers(I,:);
predictornames = drivers.Properties.VariableNames; %Names of columns in drivers_raw

%% Below is a series of steps to get rid of all the data not included in this part of the analysis 

%Take out everything but DAR Sites - only want longest time series. 
ind = drivers.Name == {'TNC'}; drivers{ind,:} = missing;
ind = (contains(drivers.SiteID,'Temporary')); drivers{ind,:} = missing;
ind = contains(drivers.SiteID,'96'); drivers{ind,:} = missing;
ind = contains(drivers.SiteID,'97'); drivers{ind,:} = missing;
ind = contains(drivers.SiteID,'98'); drivers{ind,:} = missing;
ind = drivers.Name == {'NPS'}; drivers{ind,:} = missing;
%Also, only want before the MHW, so get rid of anythying beyond that. 
ind = drivers.Year >=2015; drivers.Year(ind) = missing;clear ind
ind = ~ismissing(drivers.Year);
drivers = drivers(ind,:);

%% CHANGE THROUGH TIME
drivers_avg = drivers; %Not necessary but easier to backtrack if something goes wrong.  
%Sort by year and site ID
[A I] = sort(drivers_avg.Year); drivers_avg = drivers_avg(I,:); clear A I
[A I] = sort(drivers_avg.SiteID); drivers_avg = drivers_avg(I,:); 

%Section to setup the matrices we wish to generate. Not all of them are
%needed for this analysis. 
vars = drivers_avg.Properties.VariableNames;
PredictorNames = {'SiteID','GridID','Lat', 'Long','Depth','Year_Start','Year_End',...
    'Years','GearRegRank','CoralStart','CoralEnd','Coral_Change','Turf_Change','Macro_Change','CCA_Change',...
    'total_biomass','herbivores','browsers','grazers','scrapers',...
    'HPop15km','OSDS_Eff','OSDS_N','Golf_N','Imperv','Sediment',...
    'RainAnnSum','RainMax3dS','WPow975pct','SST_MEAN','SST_STD'};

varTypes = {'double'} %need this as many times as size varnames
varTypes = repmat(varTypes,1,size(PredictorNames,2));
site = unique(drivers_avg.SiteID);
Benthic_Change = table('Size',[size(site,1) size(PredictorNames,2)],'VariableTypes',varTypes,'VariableNames',PredictorNames);
Benthic_Change{:,:} = missing; 
%Also need to variables for coral and coral change
CoralNames = {'SiteID','GridID','Lat', 'Long','Depth','Year_2002','Year_2003','Year_2004','Year_2005','Year_2006',...
    'Year_2007','Year_2008','Year_2009','Year_2010','Year_2011','Year_2012','Year_2013','Year_2014'};
varTypes = {'double'} %need this as many times as size varnames
varTypes = repmat(varTypes,1,size(CoralNames,2));
site = unique(drivers_avg.SiteID);
Coral_Cover = table('Size',[size(site,1) size(CoralNames,2)],'VariableTypes',varTypes,'VariableNames',CoralNames);
Coral_Change = table('Size',[size(site,1) size(CoralNames,2)],'VariableTypes',varTypes,'VariableNames',CoralNames);
Coral_Cover{:,:} = missing; Coral_Change{:,:} = missing; 
%now generate empty matrices for the other fuctional groups
Turf_Change = Coral_Change; CCA_Change = Coral_Change; Macro_Change = Coral_Change; 

%Headers to identify columns we wish to average. 
fish_headers = {'total_biomass','herbivores','browsers','grazers','scrapers'};
human_headers = {'HPop15km','OSDS_Eff','OSDS_N','Golf_N','Imperv',};
naturalmean_headers = {'RainAnnSum','RainMax3dS','SST_MEAN','SST_STD'};
%The difference is how we want to average them. The above will be a
%long-term average, the below will be based on peak events.  
peak_headers = {'Sediment','WPow975pct'};

for i = 1:size(site,1);
    ind = find(drivers_avg.SiteID == site(i));%get site 1.
    [yearmx Imx] = max(drivers_avg.Year(ind));
    [yearmn Imn] = min(drivers_avg.Year(ind));
Benthic_Change.SiteID(i) = site(i);
Benthic_Change.GridID(i) = drivers_avg.GRID_ID_100m(ind(1));
Benthic_Change.Lat(i) = drivers_avg.Lat(ind(1));
Benthic_Change.Long(i) = drivers_avg.Long(ind(1));
Benthic_Change.Depth(i) = drivers_avg.Depth_survey(ind(1));
Benthic_Change.Year_Start(i) = yearmn;
Benthic_Change.Year_End(i) = yearmx;
Benthic_Change.Years(i) = yearmx - yearmn;
Benthic_Change.CoralStart(i) = round(drivers_avg.Coral(ind(Imn)),1);% starting coral cover
Benthic_Change.CoralEnd(i) = round(drivers_avg.Coral(ind(Imx)),1);% starting coral cover
Benthic_Change.GearRegRank(i) = drivers_avg.MPA_25m_GearRegRank(ind(1));
%now calculating functional groups
Benthic_Change.Coral_Change(i) = round(drivers_avg.Coral(ind(Imx)),1) - round(drivers_avg.Coral(ind(Imn)),1);
Benthic_Change.Turf_Change(i) = round(drivers_avg.Turf(ind(Imx)),1) - round(drivers_avg.Turf(ind(Imn)),1);
Benthic_Change.Macro_Change(i) = round(drivers_avg.Macro(ind(Imx)),1) - round(drivers_avg.Macro(ind(Imn)),1);
Benthic_Change.CCA_Change(i) = round(drivers_avg.Calg(ind(Imx)),1)- round(drivers_avg.Calg(ind(Imn)),1);


%% FISH <><><><><><
% now calculating fish data. 
for j = 1:length(fish_headers)
  
col2avg = [find(contains(vars,join([fish_headers(j),mat2str(yearmn)],'_'))):1:find(contains(vars,join([fish_headers(j),mat2str(yearmx)],'_')))];
col2insert = find(contains(PredictorNames,fish_headers(j)));
fishdata = [drivers_avg{ind(1),col2avg}]; % Separates data out
fishavg = nanmean(fishdata); 
Benthic_Change{i,col2insert} = fishavg;
end

%% LAND-BASED DRIVERS

%find columns to average over 
for k = 1:length(human_headers)
    yr2start = 2000;%IMPORTANT - makes the starting year the first 3 years before the start of the fish/benthic data set
    col2avg = [find(contains(vars,join([human_headers(k),mat2str(yr2start)],'_'))):1:find(contains(vars,join([human_headers(k),mat2str(yearmx)],'_')))];
    col2insert = find(contains(PredictorNames,human_headers(k)));
    humandata = [drivers_avg{ind(1),col2avg}]; % Separates data out
    humanavg = (nanmean(humandata)); 
    Benthic_Change{i,col2insert} = humanavg;
end

%% ENVIRONMENTAL DRIVERS
for m = 1:length(naturalmean_headers)
    yr2start = 2000;%
    col2avg = [find(contains(vars,join([naturalmean_headers(m),mat2str(yr2start)],'_'))):1:find(contains(vars,join([naturalmean_headers(m),mat2str(yearmx)],'_')))];
    col2insert = find(contains(PredictorNames,naturalmean_headers(m)));
    naturaldata = [drivers_avg{ind(1),col2avg}]; % Separates data out
    naturalavg = nanmean(naturaldata); 
    %AVERAGE ALL YEARS
    Benthic_Change{i,col2insert} = naturalavg; %
end
clear m 
%% PEAK DRIVERS 
for m = 1:length(peak_headers)
    yr2start = 2000;
    col2avg = [find(contains(vars,join([peak_headers(m),mat2str(yr2start)],'_'))):1:find(contains(vars,join([peak_headers(m),mat2str(yearmx)],'_')))];
    col2insert = find(contains(PredictorNames,peak_headers(m)));
    peakdata = [drivers_avg{ind(1),col2avg}]; % Separates data out
    %Finds the top 3 events and takes the average. 
    peakmax = mean(maxk(peakdata,5)); 
    Benthic_Change{i,col2insert} = (peakmax);
end

%Fill in data for Coral/Turf/CCA/Macro Change and Coral Cover
%Same metadata values go into Change and Coral_over
Coral_Change(i,1) = Benthic_Change(i,1);Coral_Cover(i,1) = Benthic_Change(i,1);
Coral_Change{i,2:5} = Benthic_Change{i,2:5}; Coral_Cover{i,2:5} = Benthic_Change{i,2:5};
Turf_Change(i,1) = Benthic_Change(i,1); Turf_Change{i,2:5} = Benthic_Change{i,2:5};
CCA_Change(i,1) = Benthic_Change(i,1); CCA_Change{i,2:5} = Benthic_Change{i,2:5};
Macro_Change(i,1) = Benthic_Change(i,1); Macro_Change{i,2:5} = Benthic_Change{i,2:5};

yr = (drivers_avg.Year(ind));
for z = 1:size(yr,1)
    col = find(contains(CoralNames,mat2str(yr(z)))); %Gets column to insert
    Coral_Cover{i,col} = round(drivers_avg.Coral(ind(z)),1);
    Coral_Change{i,col} = round(drivers_avg.Coral(ind(z)),1) - round(drivers_avg.Coral(ind(1)),1); %Removes the first year of rthe data to calculate 'change'
    Turf_Change{i,col} =  round(drivers_avg.Turf(ind(z)),1) -  round(drivers_avg.Turf(ind(1)),1);
     CCA_Change{i,col} = round(drivers_avg.Calg(ind(z)),1) - round(drivers_avg.Calg(ind(1)),1);
    Macro_Change{i,col} = round(drivers_avg.Macro(ind(z)),1) - round(drivers_avg.Macro(ind(1)),1);
end

clear ind yearmx Imx yearmn Imn col* j k m human natural
end
%% Import CHL data
%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 9);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:I27";

% Specify column names and types
opts.VariableNames = ["SiteID", "Name", "Day", "Month", "Year", "Depth_survey", "Lat", "Long", "CHL_CLIM_M"];
opts.VariableTypes = ["double", "categorical", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Name", "EmptyFieldRule", "auto");

% Import the data
CHL = readtable("/Users/jgove/Documents/jamisongove/EOD/KONA_IEA/PaperI/data/Analysis/Final_Data/Revision/Submission/CHL_PRE.DISTURBANCE.xlsx", opts, "UseExcel", false);

%% Clear temporary variables
clear opts

%% Continue Analysis
[A I] = sort(CHL.SiteID); CHL = CHL(I,:); %Sorts by SiteID
%remove sites 96, 96 , 98 as they are not included in the above due to
%length of time series; 
ind = CHL.SiteID == 96; CHL(ind,:) = [];
ind = CHL.SiteID == 97; CHL(ind,:) = [];
ind = CHL.SiteID == 98; CHL(ind,:) = [];
%Now sort Benthic Change - this will align CHL and Benthic Change
[B J] = sort(Benthic_Change.SiteID);
Benthic_Change = Benthic_Change(J,:); 

%ADD IN CHL METRICS 
Benthic_Change = [Benthic_Change CHL(:,end)]
%clean up workspace

clearvars -except drivers WESTHAWAIIALLDRIVERS Benthic_Change Coral_Change Coral_Cover
%% CLEANING UP
%Now make a series of adustments for each predictor. 
BC = Benthic_Change; 
BC.Sediment = BC.Sediment*1000; %Changes to tons/ha

%Adjust driver resoltion to appropriate sig digits. 
BC.total_biomass = round(BC.total_biomass,2);
BC.herbivores = round(BC.herbivores,2);
BC.grazers = round(BC.grazers,2);
BC.scrapers = round(BC.scrapers,2);
BC.browsers =  round(BC.browsers,2);
BC.HPop15km = round(BC.HPop15km,1);
BC.OSDS_Eff = round(BC.OSDS_Eff,1);
BC.OSDS_N = round(BC.OSDS_N,1);
BC.Golf_N= round(BC.Golf_N,1);
BC.Nuts = BC.Golf_N + BC.OSDS_N; %create total nutrients variable. 
BC.Imperv = round(BC.Imperv,1);
BC.Sediment = round(BC.Sediment,1);
BC.RainAnnSum = round(BC.RainAnnSum,1);
BC.RainMax3dS = round(BC.RainMax3dS,1);
BC.WPow975pct = round(BC.WPow975pct,1);
BC.SST_STD = round(BC.SST_STD,2);
BC.SST_MEAN= round(BC.SST_MEAN,2);
BC.Depth = round(BC.Depth,1);

%% VISUALIZE DISTRIBUTIONS AND ID OUTLIERS 
%This produces a plot for each predictor. Simply revove the % if you wish
%to do this. 
% vars = BC.Properties.VariableNames;
%vars2plot = {'Coral_Change','CoralStart','CoralEnd','Depth','GearRegRank','total_biomass',...
%'herbivores','browsers','grazers','scrapers','HPop15km','OSDS_Eff','Nuts','Imperv','Sediment',...
%'RainAnnSum','RainMax3dS','WPow975pct','SST_MEAN','SST_STD','CHL_CLIM_M'};
 
% for i = 1:numel(vars2plot)
%    col = matches(vars,vars2plot(i)); %find the column of data to plot. 
%    figure('Renderer', 'painters', 'Position', [500 500 900 400])
%    subplot(2,1,1);  histogram(BC{:,col},8,'normalization','probability');
%    hold on
%    pos_thresh = nanmedian(BC{:,col}) + 3*std(BC{:,col},'omitnan');
%    %xline(pos_thresh,'r'); 
%      title(vars2plot(i))
%      ylabel('Probability')
%     subplot(2,1,2);
%    plot(BC{:,col},BC.Coral_Change,'k.','markersize',30);hold on
%    xline(pos_thresh,'r'); )
%     indbelow = BC{:,col} <pos_thresh; %indices below threshold
%     A = max(BC{indbelow,col}); %gets max value below the threshold. 
%     B = A*1.25; %Creates line at 25% mark 
%     xline(B,'g','linewidth',2); %plots line 
%    hold on
%    ylabel('Coral Change')
%  
% end


%% GET RID OF OUTLIERS (see Methods for outlier ID approach)
ind = BC.RainMax3dS >22800; BC.RainMax3dS(ind)= NaN;
ind = BC.Imperv >36000; BC.Imperv(ind)= NaN; 
ind = BC.Sediment >19500; BC.Sediment(ind) = NaN; 
ind = BC.OSDS_Eff >1800000; BC.OSDS_Eff(ind)= NaN;
ind = BC.Nuts > 150; BC.Nuts(ind) = NaN; 
ind = BC.grazers> 35; BC.grazers(ind) = NaN; 
ind = BC.total_biomass> 200; BC.total_biomass(ind) = NaN; 

clearvars -except drivers BC Benthic_Change Coral_Change Coral_Cover

%% PLOTTING 
%%%%%%%%%%%%%%%%%FIGURE 2A%%%%%%%%%%%%%%%%%%%
X = drivers.Coral; %define variable. 
%Figures
figure; hold on
ind = drivers.Year ==2003; X03 = X(ind); 
ind = drivers.Year ==2007; X07 = X(ind); 
ind = drivers.Year ==2011; X11 = X(ind); 
ind = drivers.Year ==2014; X14 = X(ind);

[fi,xi] = ksdensity(X03,'Bandwidth',7); fi = fi*10; 
p= plot(xi,fi,'linewidth',2); 

[fi,xi] = ksdensity(X07,'Bandwidth',7);fi = fi*10;
p= plot(xi,fi,'linewidth',2); 

[fi,xi] = ksdensity(X11,'Bandwidth',7);fi = fi*10;
p= plot(xi,fi,'linewidth',2); 

[fi,xi] = ksdensity(X14,'Bandwidth',7);fi = fi*10;
p= plot(xi,fi,'linewidth',2); 

set(gca,'xtick',[0:5:80])
xlim([0 80])
ylim([0 0.32]); set(gca,'ytick',[0:0.05:0.3]); 

ylabel('Relative Distribution')
xlabel('Coral Cover (%)')
legend('2003','2007','2011','2014')

clearvars -except drivers BC Benthic_Change Coral_Change Coral_Cover

%%%%%%%%%%%%%%FIGURE 2B %%%%%%%%%%%%%%
[A I] = sort(Coral_Change.SiteID); 
Coral_Change = Coral_Change(I,:); 
years = 2002:2014;

X = BC;
Y = Coral_Change;

figure('Renderer', 'painters', 'units','centimeters','Position', [50 50 12 8]); hold on
map = brewermap(100,'RdBu');
for i = 1:size(Y,1)
idx = 6:18; %colums with data)
a = find(~isnan(Y{i,6:18})); %not all columns have data. 
if X.Coral_Change(i) <=-3
p1 = plot(years(a),Y{i,idx(a)},'-','color',map(12,:),'linewidth',5);
p2 = plot(years(a),Y{i,idx(a)},'.','color',map(12,:),'markersize',20);
p1.Color(4) = 0.2; p2.Color(4) = 0.1;
elseif X.Coral_Change(i) >=3
p1 = plot(years(a),Y{i,idx(a)},'-','color',map(90,:),'linewidth',5);
p2 = plot(years(a),Y{i,idx(a)},'.','color',map(90,:),'markersize',20);
p1.Color(4) = 0.2; p2.Color(4) = 0.1;
end
end
set(gca,'xtick',[2003:1:2014]); 
set(gca,'ytick',[-20:5:20])
axis tight
grid on

ylabel('Coral Cover Change (%)')
xlabel('Year')
ylim([-21 20])

clearvars -except drivers BC Benthic_Change Coral_Change Coral_Cover


%%%%%%%%%%%%%%%%%%%% FIGURE 2C PERCENT DIFFERENCE %%%%%%%%%%%

matrix2avg = BC; 
matrix_headers = matrix2avg.Properties.VariableNames;
drivers2avg = {'browsers','scrapers','HPop15km','grazers','herbivores','total_biomass','WPow975pct',...
'OSDS_Eff','Nuts','Imperv','RainMax3dS','RainAnnSum','Sediment',...
'SST_MEAN','SST_STD','CHL_CLIM_M','GearRegRank','Depth'};

%get positive and negative trajectory sites. 
ind_pos = matrix2avg.Coral_Change >3;
ind_neg = matrix2avg.Coral_Change <-3;
Delta = [];

% PERCENT DIFFERENCE
for i = 1:numel(drivers2avg)
    col = matches(matrix_headers,drivers2avg(i)); %find the column of data
    YY = jackknife(@nanmean,matrix2avg{ind_pos,col}); %get all positive data
     NN = jackknife(@nanmean,matrix2avg{ind_neg,col});%get all negative data
 
if mean(YY) - mean(NN) > 0   
P = [mean(YY)-mean(NN)]/[(mean(YY)+mean(NN))/2];
Pmin = [min(YY)-max(NN)]/[(min(YY)+max(NN))/2];
Pmax = [max(YY)-min(NN)]/[(max(YY)+min(NN))/2];

Delta(i,1) = P;
Delta(i,2) = P-Pmin; %Lower bound
Delta(i,3) = Pmax-P; %Upper bound
clear YY NN
else
    %Ensure values are negative, also makes 'lower bound' as the maximum
    %negative value, so min and max switch compared to postive values. 
P = [mean(NN)-mean(YY)]/[(mean(NN)+mean(YY))/2];
Pmin = [(min(NN)-max(YY))]/[(min(NN)+max(YY))/2];
Pmax = [(max(NN)-min(YY))]/[(max(NN)+min(YY))/2];
Delta(i,1) = P*-1;
Delta(i,2) = (Pmax-P)*-1; %More negative bound
Delta(i,3) = (P-Pmin)*-1; %Less negative bound 
end
end
clearvars -except BC Benthic_Change drivers Delta matrix2avg drivers2avg
% 

%PLOTTING: This order is custom based on visualizing the results. 
plot_order = {'browsers','scrapers','HPop15km','grazers','herbivores','total_biomass','WPow975pct',...
'OSDS_Eff','Nuts','Imperv','RainMax3dS','RainAnnSum','Sediment',...
'SST_MEAN','SST_STD','CHL_CLIM_M','GearRegRank','Depth'};
%Changes things to cm so you can better size graphics. 
figure('Renderer', 'painters', 'units','centimeters','Position', [50 50 18 10]); hold on

Delta = Delta.*100; %Turns it into a percentage
for i = 1:numel(plot_order);
    ind = matches(drivers2avg,plot_order(i));
    e = errorbar(i,Delta(ind,1),[Delta(ind,2)],[Delta(ind,3)],'vertical','o','Markersize',10,'MarkerEdgeColor','k');        
        e.LineWidth = 3;
        e.CapSize = 0;
        e.Color = 'r';
end
set(gca,'Xtick',1:numel(plot_order));        
set(gca,'XtickLabel',plot_order) 
ylim([-125 125])
set(gca,'ytick',[-125:25:125])
xlabel('Land-Sea Driver'); ylabel('Percent Difference','FontSize',20)
title('3% Coral Change Cutoff')

% Calculate cutoff for 'minimal difference' between trajectories.  
X = Delta(:,1); %All mean differences. 
ZZ = bootstrp(10000,@median,X); %Bootstrap
cutoff = std(ZZ)*2; %Calculate 2x SD;
%Draw cutoff lines. 
yline(cutoff,':k','linewidth',2); yline(-cutoff,':k','linewidth',2)
%Note that figure refining was all done in AI. 
