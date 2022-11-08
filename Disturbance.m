%% CORAL RESPONSE TO THE 2015 MARINE HEATWAVE
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


%% Clear temporary variables
clear opts
%% DATA ANALYSIS %% 
%Make sure to rename it: drivers_raw
drivers = drivers_raw;
[A I] = sort(drivers.Year); drivers = drivers(I,:);
[A I] = sort(drivers.SiteID); drivers = drivers(I,:);
predictornames = drivers.Properties.VariableNames; %Names of columns in drivers_raw

%% PRE-CLEANING
%We only want data that is from 2014 to 2016 
ind = drivers.Year <=2013; drivers.Year(ind) = missing;clear ind
ind = drivers.Year >=2017; drivers.Year(ind) = missing;clear ind
ind = ~ismissing(drivers.Year); 
drivers = drivers(ind,:);
%Site 96 does not have fish data...remove. 
ind = (drivers.SiteID == '96');
drivers{ind,:} = missing;
ind = ~ismissing(drivers.Year); 
drivers = drivers(ind,:);

%% CHANGE THROUGH TIME 
%Section builds the variable to populate
drivers_avg = drivers;
vars = drivers_avg.Properties.VariableNames;
PredictorNames = {'SiteID','GridID','Lat', 'Long','Depth','Year_Start','Year_End','Years','GearRegRank',...
    'CoralStart','CoralEnd','CCAStart','Coral_Change','Turf_Change','Macro_Change','CCA_Change',...
    'total_biomass','herbivores','browsers','grazers','scrapers',...
    'HPop15km','OSDS_Eff','OSDS_N','Golf_N','Imperv','Sediment','RainAnnSum','RainMax3dS',...
    'WPow975pct','SST_MEAN','SST_STD','DHW_MAX_2015','SST_MAX_2015','SSTA_MAX_2015','HS_MAX_2015',...
    'CHL_MHW','PAR_MHW'};

varTypes = {'double'} %need this as many times as size varnames
varTypes = repmat(varTypes,1,size(PredictorNames,2));
site = unique(drivers_avg.SiteID);
Benthic_Change = table('Size',[size(site,1) size(PredictorNames,2)],'VariableTypes',varTypes,'VariableNames',PredictorNames);
Benthic_Change{:,:} = missing; 
%Also need to variables for coral and coral change
CoralNames = {'SiteID','GridID','Lat', 'Long','Depth','Year_2015','Year_2016'};
varTypes = {'double'} ;%need this as many times as size varnames
varTypes = repmat(varTypes,1,size(CoralNames,2));
site = unique(drivers_avg.SiteID);
Coral_Cover = table('Size',[size(site,1) size(CoralNames,2)],'VariableTypes',varTypes,'VariableNames',CoralNames);
Coral_Change = table('Size',[size(site,1) size(CoralNames,2)],'VariableTypes',varTypes,'VariableNames',CoralNames);
Coral_Cover{:,:} = missing; Coral_Change{:,:} = missing;
Benthic_Change = convertvars(Benthic_Change,{'SiteID'},'categorical');%Need to change these columns to categorical. 
Coral_Change = convertvars(Coral_Change,{'SiteID'},'categorical');%Need to change these columns to categorical. 
Coral_Cover = convertvars(Coral_Cover,{'SiteID'},'categorical');%Need to change these columns to categorical. 

%% Calculating the change over time for each predictor

fish_headers = {'total_biomass','herbivores','browsers','grazers','scrapers'};
human_headers = {'HPop15km','Imperv','OSDS_Eff','OSDS_N','Golf_N'};
natural_headers = {'RainAnnSum','RainMax3dS','WPow975pct'};
other_headers = {'Sediment'};
resilience_headers  = {'SST_MEAN','SST_STD'};
for i = 1:size(site,1);
    ind = find(drivers_avg.SiteID == site(i));
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
Benthic_Change.GearRegRank(i) = drivers_avg.MPA_25m_GearRegRank(ind(1));
Benthic_Change.DHW_MAX_2015(i) = drivers_avg.DHW_MAX_2015(ind(1));
Benthic_Change.SST_MAX_2015(i) = drivers_avg.SST_MAX_2015(ind(1));
Benthic_Change.SSTA_MAX_2015(i) = drivers_avg.SSTA_MAX_2015(ind(1));
Benthic_Change.HS_MAX_2015(i) = drivers_avg.HS_MAX_2015(ind(1));
Benthic_Change.CHL_MHW(i) = drivers_avg.CHL_MHW(ind(1));
Benthic_Change.PAR_MHW(i) = drivers_avg.PAR_MHW(ind(1));
Benthic_Change.CoralStart(i) = round(drivers_avg.Coral(ind(Imn)),1);% starting coral cover
Benthic_Change.CoralEnd(i) = round(drivers_avg.Coral(ind(Imx)),1);% starting coral cover
Benthic_Change.CCAStart(i) = round(drivers_avg.Calg(ind(Imx)),1);% starting coral cover
%now calculating functional groups
Benthic_Change.Coral_Change(i) = round(drivers_avg.Coral(ind(Imx)),1) - round(drivers_avg.Coral(ind(Imn)),1);
Benthic_Change.Turf_Change(i) = round(drivers_avg.Turf(ind(Imx)),1) - round(drivers_avg.Turf(ind(Imn)),1);
Benthic_Change.Macro_Change(i) = round(drivers_avg.Macro(ind(Imx)),1) - round(drivers_avg.Macro(ind(Imn)),1);
Benthic_Change.CCA_Change(i) = round(drivers_avg.Calg(ind(Imx)),1)- round(drivers_avg.Calg(ind(Imn)),1);

%% FISH <><><><><><
% now calculating fish data. 
for j = 1:length(fish_headers);
 
col2avg = [find(contains(vars,join([fish_headers(j),mat2str(yearmn)],'_'))) find(contains(vars,join([fish_headers(j),mat2str(yearmx)],'_')))];
col2insert = find(contains(PredictorNames,fish_headers(j)));
fishdata = [drivers_avg{ind(1),col2avg}]; % Separates data out
fishavg = nanmean(fishdata); 
Benthic_Change{i,col2insert} = fishavg;
end

%% HUMAN DRIVERS

for k = 1:length(human_headers);
    yr2start = 2012; %Starts 3 years before MWH 
    col2avg = [find(contains(vars,join([human_headers(k),mat2str(yr2start)],'_'))):1:find(contains(vars,join([human_headers(k),mat2str(yearmx)],'_')))];
    col2insert = find(contains(PredictorNames,human_headers(k)));
    humandata = [drivers_avg{ind(1),col2avg}]; % Separates data out, Only need 1 row as all rows will be the same data
    humanavg = nanmean(humandata);
    Benthic_Change{i,col2insert} = humanavg;
end

%% NATURAL DRIVERS
for m = 1:length(natural_headers)
yr2start = 2012;%Starts 3 years before MWH 
    col2avg = [find(contains(vars,join([natural_headers(m),mat2str(yr2start)],'_'))):1:find(contains(vars,join([natural_headers(m),mat2str(yearmx)],'_')))];
    col2insert = find(contains(PredictorNames,natural_headers(m)));
    naturaldata = [drivers_avg{ind(1),col2avg}]; % Separates data out, Only need 1 row as all rows will be the same data
    naturalavg = nanmean(naturaldata);
    Benthic_Change{i,col2insert} = naturalavg;
end
clear m

%% SEDIMENT
for m = 1:length(other_headers);
    yr2start = 2006;
    col2avg = [find(contains(vars,join([other_headers(m),mat2str(yr2start)],'_'))):1:find(contains(vars,join([other_headers(m),mat2str(yearmx)],'_')))];
    col2insert = find(contains(PredictorNames,other_headers(m)));
    otherdata = [drivers_avg{ind(1),col2avg}]; % Separates data out, Only need 1 row as all rows will be the same data
    othermax= mean(maxk(otherdata,3));
    Benthic_Change{i,col2insert} = othermax;
end
clear m
 %% RESILIENCE DRIVERS
resilience_headers  = {'SST_MEAN','SST_STD'};
for m = 1:length(resilience_headers);
    yr2start = 2000;
    yr2end = 2014; 
    col2avg = [find(contains(vars,join([resilience_headers(m),mat2str(yr2start)],'_'))):1:find(contains(vars,join([resilience_headers(m),mat2str(yr2end)],'_')))];
    col2insert = find(contains(PredictorNames,resilience_headers(m)));
    resiliencedata = [drivers_avg{ind(1),col2avg}]; % Separates data out, Only need 1 row as all rows will be the same data
    resilienceavg = nanmean(resiliencedata);
    Benthic_Change{i,col2insert} = resilienceavg;
end
%Fill in data for Coral Change and Coral Cover
%Same metadata values go into Coral_Change and Coral_over
Coral_Change(i,1) = Benthic_Change(i,1);Coral_Cover(i,1) = Benthic_Change(i,1);
Coral_Change{i,2:5} = Benthic_Change{i,2:5}; Coral_Cover{i,2:5} = Benthic_Change{i,2:5};
yr = drivers_avg.Year(ind);
for z = 1:size(yr,1)
    col = find(contains(CoralNames,mat2str(yr(z)))); %Gets column to insert
    Coral_Cover{i,col} = drivers_avg.Coral(ind(z));
    Coral_Change{i,col} = drivers_avg.Coral(ind(z)) - drivers_avg.Coral(ind(1)); %Removes the first year of rthe data to calculate 'change'
end
clear ind yearmx Imx yearmn Imn col* j k m human natural
end
%clean up workspace
clearvars -except Benthic_Change Coral_Cover drivers drivers_avg drivers_raw 

%Now make a series of adustments for each predictor. 
BC = Benthic_Change; %Get new variable to preserve all work above
BC.Sediment = BC.Sediment*1000; %Changes sediment to tons/ha
%Need to adjust driver resolution to sig digits that make sense
%First fish
BC.total_biomass = round(BC.total_biomass,2);
BC.herbivores = round(BC.herbivores,2);
BC.grazers = round(BC.grazers,2);
BC.scrapers = round(BC.scrapers,2);
BC.browsers =  round(BC.browsers,2);
%Humans
BC.HPop15km = round(BC.HPop15km,1);
BC.OSDS_Eff = round(BC.OSDS_Eff);
BC.OSDS_N = round(BC.OSDS_N,1);
BC.Golf_N= round(BC.Golf_N,1);
BC.Nuts = BC.OSDS_N+BC.Golf_N; %ADDING NEW VARIABLE 
BC.Imperv = round(BC.Imperv,1);
BC.Sediment = round(BC.Sediment,1);
BC.RainAnnSum = round(BC.RainAnnSum,1);
BC.RainMax3dS = round(BC.RainMax3dS,1);
BC.WPow975pct = round(BC.WPow975pct,2);
BC.DHW_MAX_2015 = round(BC.DHW_MAX_2015,2);
BC.SST_STD = round(BC.SST_STD,3);
BC.CHL_MHW= round(BC.CHL_MHW,3);
BC.Depth = round(BC.Depth,1);
BC.RelChange = [(BC.CoralEnd-BC.CoralStart)./BC.CoralStart*100];%Calculates 

%Remove unnecessary values and shift variables around. 
BC = removevars(BC, {'SiteID','GridID','Year_Start','Year_End','Years',...
    'OSDS_N','Golf_N','CoralStart','CoralEnd','Macro_Change','Turf_Change','CCA_Change','CCAStart',});
BC  = removevars(BC, {'SST_MAX_2015','SSTA_MAX_2015','HS_MAX_2015'}); %None of these are needed. 
BC = movevars(BC, 'Coral_Change', 'Before', 'Lat');
BC = movevars(BC, 'RelChange', 'Before', 'Lat');
BC = movevars(BC, 'Nuts', 'Before', 'Imperv');

%% VISUALIZE DISTRIBUTIONS AND OUTLIERS 
%This produces a plot for each predictor. Simply revove the % if you wish
%to do this. 
% vars = BC.Properties.VariableNames;
% vars2plot = vars;
% for i = 1:numel(vars); 
%    col = matches(vars,vars2plot(i)); %find the column of data to plot. 
%    figure('Renderer', 'painters', 'Position', [500 500 600 300])
%    subplot(2,1,1);  histogram(BC{:,col},6,'normalization','probability');
%    hold on
%    pos_thresh = nanmedian(BC{:,col}) + 3*std(BC{:,col},'omitnan');
%    xline(pos_thresh,'r'); 
%      title(vars2plot(i))
%      ylabel('Probability')
%     subplot(2,1,2);
%    plot(BC{:,col},BC.RelChange,'r.','markersize',40);hold on
%    xline(pos_thresh,'r'); 
%    hold on
%    indbelow = BC{:,col} <pos_thresh; %indices below threshold
% A = max(BC{indbelow,col}); %gets max value below the threshold. 
% B = A*1.25; %Additional threshold for inclusion
% xline(B,'g','linewidth',2); 
% xline(pos_thresh,'r');
%    ylabel('Coral Change (% difference)'); 
%    xlabel(vars2plot(i)); title(vars2plot(i))
% end
%% REMOVE OUTLIERS AND TRANSFORM
%See methods for approach or uncomment the above and plot all
%distributions. 
ind = BC.Sediment > 13000; BC.Sediment(ind) = NaN; 
ind = BC.Imperv > 31000; BC.Imperv(ind) = NaN; 
ind = BC.Nuts > 140; BC.Nuts(ind) = NaN;
ind = BC.OSDS_Eff >1800000; BC.OSDS_Eff(ind) = NaN;
BC.OSDS_Eff = sqrt(BC.OSDS_Eff); 
BC.Nuts= sqrt(BC.Nuts); 
BC.Sediment = nthroot(BC.Sediment,4);
BC.HPop15km = sqrt(BC.HPop15km);
BC.CHL_MHW = sqrt(BC.CHL_MHW); 
BC.RainMax3dS = sqrt(BC.RainMax3dS);
BC.Imperv = sqrt(BC.Imperv);

%FISH
ind = BC.scrapers > 40; BC.scrapers(ind) = NaN;
ind = BC.browsers > 15; BC.browsers(ind) = NaN;
ind = BC.herbivores > 150; BC.herbivores(ind) = NaN;
ind = BC.total_biomass > 800; BC.total_biomass(ind) = NaN;
BC.total_biomass = sqrt(BC.total_biomass); 
BC.herbivores = sqrt(BC.herbivores);
BC.scrapers = sqrt(BC.scrapers); 
BC.grazers = sqrt(BC.grazers);
BC.browsers = sqrt(BC.browsers); 

%% Plot Correlation Matrix (remove variables that are not predictors). 
BCmodel = BC; 
BCmodel = removevars(BCmodel,{'Lat','Long','RelChange','Coral_Change'});
plot_corr_matrix(BCmodel)

%% New varable to export
BCexport = BC;
%Remove correlated variables (see Methods)
BCexport = removevars(BCexport, {'herbivores','browsers','RainAnnSum','SST_STD','SST_MEAN','HPop15km','PAR_MHW','Coral_Change'});
cd '/Users/jgove/Documents/jamisongove/EOD/KONA_IEA/PaperI/data/Analysis/Final_Data/Revision'

writetable(BCexport,'GAMM_MHW.csv');

%%%%%%PLOTTING%%%%%%%

%% Distribution of Coral Cover Change
X = Benthic_Change.Coral_Change; 
%Figures
figure; hold on
histogram(X,'normalization','probability')
[fi,xi] = ksdensity(X,'Bandwidth',4.5); fi = fi*10; 
p= plot(xi,fi,'linewidth',3); 
ylim([0.01 0.35])
xlim([-50 12])
set(gca,'Xtick',[-55:5:10])
set(gca,'Ytick',[0:0.05:0.35])
ylabel('Proportion of Reefs')
xlabel('Coral Cover Change')
%% DHW FIGURE 3b
XX = BC.DHW_MAX_2015; 
figure('Renderer', 'painters', 'units','centimeters','Position', [50 50 12 12])
histogram(XX,'normalization','probability')
[fi,xi] = ksdensity(XX,'Bandwidth',0.38); 
hold on
plot(xi,fi,'r','linewidth',3); hold on
ha1 = area(xi,fi,'FaceColor','r'); ha1.FaceAlpha = 0.2; 
ha1.EdgeColor = 'none';
set(gca,'xtick',[9:0.5:14]);
set(gca,'xticklabel',{'9','','10','','11','','12','','13','','14'})
xlim([9 14])
ylim([0 0.52])    
ylabel('Proportion of Reefs')
xlabel('Degree Heating Weeks (C-weeks)')

%% Figure 3c
figure('Renderer', 'painters', 'units','centimeters','Position', [50 50 8 14])
hold on

for i = 1:size(Benthic_Change,1);
y = [Benthic_Change.CoralStart(i) Benthic_Change.CoralEnd(i)] 
x = [2015 2016];
col = [Benthic_Change.CoralStart(i) Benthic_Change.CoralEnd(i)]   % This is the color, vary with x in this case.
z = [0 0];
surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',6);
plot(x,y,'-','color',[0.75 0.75 0.75],'linewidth',2)
end
scatter((ones(size(Coral_Cover,1),1)*2015),Benthic_Change.CoralStart,150,Benthic_Change.CoralStart,'filled')
scatter((ones(size(Coral_Cover,1),1)*2015),Benthic_Change.CoralStart,150,[0.7 0.7 0.7])
scatter((ones(size(Coral_Cover,1),1)*2016),Benthic_Change.CoralEnd,150,Benthic_Change.CoralEnd,'filled')
scatter((ones(size(Coral_Cover,1),1)*2016),Benthic_Change.CoralEnd,150,[0.7 0.7 0.7])
xlim([2014.95 2016.05])
set(gca,'xtick',[2015 2016])
set(gca,'ytick',[0:2.5:65])
ylim([0 65])
ylabel('Coral Cover (%)')
xlabel('Year')
%NOTE Brewer map is a custom color ramp. Download it from here:
%https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
set(gcf,'renderer','Painters');
map = brewermap(80,'RdPu');
colormap(map)
caxis([5 60])

