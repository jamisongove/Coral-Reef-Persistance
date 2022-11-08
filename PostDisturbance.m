%% FOUR YEARS POST-DISTURBANCE
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
%Not necessary, but simply creating a new variable just to preserve the
%above.
drivers = drivers_raw; 
predictornames = drivers.Properties.VariableNames; %Names of columns in drivers

%Only including reefs that were surveyed up to 2019. 
ind = drivers.Year == 2019; %Note this keeps all data from first survey through 2019. 
drivers{~ind,:} = missing;
%Get rid of missing values
ind = ~ismissing(drivers.Year);drivers = drivers(ind,:);
%Not necessary, but simply creating a new variable just to preserve the
%above. 
drivers_avg = drivers; 
[A I] = sort(drivers_avg.Year); drivers_avg = drivers_avg(I,:);
[A I] = sort(drivers_avg.SiteID); drivers_avg = drivers_avg(I,:);

%% Section that finds each driver, takes an average over a set timeframe, and populates a new variable.  
vars = drivers_avg.Properties.VariableNames;
PredictorNames = {'SiteID','GridID','Lat', 'Long','Depth','Year_Start','Year_End','Years','GearRegRank',...
    'CoralStart','CoralEnd','Coral_Change','Turf_Change','Macro_Change','CCAStart','CCAEnd','CCA_Change',...
    'TCStart','TCEnd','total_biomass','herbivores','browsers','grazers','scrapers','HPop15km','OSDS_Eff',...
    'OSDS_N','Golf_N','Imperv','Sediment','RainAnnSum','RainMax3dS','WPow975pct','SST_MEAN','SST_STD',...
    'CHL','PAR'};
varTypes = {'double'} %need this as many times as size varnames
varTypes = repmat(varTypes,1,size(PredictorNames,2));
site = unique(drivers_avg.SiteID); %Gets us all of the sites. 
Benthic_Change = table('Size',[size(site,1) size(PredictorNames,2)],'VariableTypes',varTypes,'VariableNames',PredictorNames);
Benthic_Change{:,:} = missing; %Keeps the new variable without data. 
site = unique(drivers_avg.SiteID);
Benthic_Change = convertvars(Benthic_Change,{'SiteID'},'categorical');%Need to change these columns to categorical. 


%% Calculating the change over time for each predictor
fish_headers = {'total_biomass','herbivores','browsers','grazers','scrapers'};
human_headers = {'HPop15km','OSDS_Eff','OSDS_N','Golf_N','Imperv'};
natural_headers = {'RainAnnSum','RainMax3dS','WPow975pct'};
other_headers = {'Sediment'};
resilience_headers  = {'SST_STD','SST_MEAN'};
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
Benthic_Change.CHL(i) = nanmean([drivers_avg.CHL_2016(ind(1)) drivers_avg.CHL_2017(ind(1)) drivers_avg.CHL_2018(ind(1)) drivers_avg.CHL_2019(ind(1))]);
Benthic_Change.PAR(i) = nanmean([drivers_avg.PAR_2016(ind(1)) drivers_avg.PAR_2017(ind(1)) drivers_avg.PAR_2018(ind(1)) drivers_avg.PAR_2019(ind(1))]);
Benthic_Change.CoralStart(i) = round(drivers_avg.Coral(ind(Imn)),1);% starting coral cover
Benthic_Change.CoralEnd(i) = round(drivers_avg.Coral(ind(Imx)),1);% ending coral cover
Benthic_Change.CCAStart(i) = round(drivers_avg.Calg(ind(Imn)),1);% starting CCA cover
Benthic_Change.CCAEnd(i) = round(drivers_avg.Calg(ind(Imx)),1);% ending coral cover
Benthic_Change.TCStart(i) = [Benthic_Change.CCAStart(i)+Benthic_Change.CoralStart(i)]; %Starting Reef Builder Cover
Benthic_Change.TCEnd(i) = [Benthic_Change.CCAEnd(i)+Benthic_Change.CoralEnd(i)];%Ending Reef Builder Cover
%now calculating functional groups
Benthic_Change.Coral_Change(i) = round(drivers_avg.Coral(ind(Imx)),1) - round(drivers_avg.Coral(ind(Imn)),1);
Benthic_Change.Turf_Change(i) = round(drivers_avg.Turf(ind(Imx)),1) - round(drivers_avg.Turf(ind(Imn)),1);
Benthic_Change.Macro_Change(i) = round(drivers_avg.Macro(ind(Imx)),1) - round(drivers_avg.Macro(ind(Imn)),1);
Benthic_Change.CCA_Change(i) = round(drivers_avg.Calg(ind(Imx)),1)- round(drivers_avg.Calg(ind(Imn)),1);

%% FISH <><><><><><
for j = 1:length(fish_headers)
   year2start = 2016; 
col2avg = [find(contains(vars,join([fish_headers(j),mat2str(year2start)],'_'))):1:find(contains(vars,join([fish_headers(j),mat2str(yearmx)],'_')))];
col2insert = find(contains(PredictorNames,fish_headers(j)));
fishdata = [drivers_avg{ind(1),col2avg}]; % Separates data out
fishavg = nanmean(fishdata); 
Benthic_Change{i,col2insert} = fishavg;end

%% HUMAN DRIVERS
for k = 1:length(human_headers)
    yr2start = 2016; 
    col2avg = [find(contains(vars,join([human_headers(k),mat2str(yr2start)],'_'))):1:find(contains(vars,join([human_headers(k),mat2str(yearmx)],'_')))];
    col2insert = find(contains(PredictorNames,human_headers(k)));
    humandata = [drivers_avg{ind(1),col2avg}]; % Separates data out
    humanavg = mean(humandata);
    Benthic_Change{i,col2insert} = humanavg;
end

%% NATURAL DRIVERS
for m = 1:length(natural_headers)
yr2start = 2016;
    col2avg = [find(contains(vars,join([natural_headers(m),mat2str(yr2start)],'_'))):1:find(contains(vars,join([natural_headers(m),mat2str(yearmx)],'_')))];
    col2insert = find(contains(PredictorNames,natural_headers(m)));
    naturaldata = [drivers_avg{ind(1),col2avg}]; % Separates data out
    naturalavg = nanmean(naturaldata);
    Benthic_Change{i,col2insert} = naturalavg;

end
clear m
%% Sediment
for m = 1:length(other_headers)
yr2start = 2006;
    col2avg = [find(contains(vars,join([other_headers(m),mat2str(yr2start)],'_'))):1:find(contains(vars,join([other_headers(m),mat2str(yearmx)],'_')))];
    col2insert = find(contains(PredictorNames,other_headers(m)));
    %window = yearmx - yearmn;
    otherdata = [drivers_avg{ind(1),col2avg}]; % Separates data out
    othermax= mean(maxk(otherdata,3)); %finds meean of the top 3 events over this time window. 
    Benthic_Change{i,col2insert} = othermax;
end
clear m
 %% RESELIENCE DRIVERS
for m = 1:length(resilience_headers)
    yr2start = 2000;
    yr2end = 2018; 
    col2avg = [find(contains(vars,join([resilience_headers(m),mat2str(yr2start)],'_'))):1:find(contains(vars,join([resilience_headers(m),mat2str(yr2end)],'_')))];
    col2insert = find(contains(PredictorNames,resilience_headers(m)));
    resiliencedata = [drivers_avg{ind(1),col2avg}]; % Separates data out
    resilienceavg = nanmean(resiliencedata);
    Benthic_Change{i,col2insert} = resilienceavg;
end
clear ind yearmx Imx yearmn Imn col* j k m human natural
end
%clean up
clearvars -except BC Benthic_Change drivers drivers_avg

%Create new variable to preserve the old
BC = Benthic_Change; 
BC.Sediment = BC.Sediment*1000; %Changes to tons/ha
%Need to adjust driver resolution to sig digits that make sense
%First fish, rounding to hundreths
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
BC.Nuts = BC.OSDS_N+BC.Golf_N; %Create total nutrient variable, already rounded
BC.Imperv = round(BC.Imperv,1);
BC.Sediment = round(BC.Sediment,1);
BC.RainAnnSum = round(BC.RainAnnSum,1);
BC.RainMax3dS = round(BC.RainMax3dS,1);
BC.WPow975pct = round(BC.WPow975pct,2);
BC.SST_STD = round(BC.SST_STD,3);
BC.Depth = round(BC.Depth,1);

%% REMOVING VARIABLES
BC = removevars(BC, {'Lat','Long','CoralStart','CoralEnd',...
    'Coral_Change','TCStart','CCAStart','CCAEnd','CCA_Change','SiteID','GridID','Year_Start',...
    'Year_End','Years','Turf_Change','Macro_Change','OSDS_N','Golf_N'});
% 
%% VISUALIZE DISTRIBUTIONS AND ID OUTLIERS 
%This produces a plot for each predictor. Simply revove the % if you wish
%to do this. 

% vars = BC.Properties.VariableNames;
% thresh = round(prctile(BC.TCEnd,[25 75]))
% for i = 1:numel(vars); 
%      col = matches(vars,vars(i));
%       figure('Renderer', 'painters', 'Position', [500 500 800 600])
%       subplot(2,1,1);  histogram(BC{:,col},6,'normalization','probability');
%    hold on
%       pos_thresh = nanmedian(BC{:,col}) + 3*std(BC{:,col},'omitnan');
%       xline(pos_thresh,'r'); 
%       title(vars(i))
%       subplot(2,1,2);
%       plot(BC{:,col},BC.TCEnd,'k.','markersize',30);hold on
%       yline(thresh(1),'--k','linewidth',1)
%       yline(thresh(2),'--k','linewidth',1)
%       indbelow = BC{:,col} <pos_thresh; %indices below threshold
%       A = max(BC{indbelow,col}); %gets max value below the threshold. 
%       B = A*1.25; 
%       xline(B,'g','linewidth',2); 
%       xline(pos_thresh,'r'); %xline(neg_thresh,'r');
%       ylabel('Total Calcification')
%       title(vars(i))
% end

%% GET RID OF OUTLIERS (see Methods for outlier ID approach). 
ind = BC.CHL > 0.2; BC.CHL(ind) = NaN;
ind = BC.RainMax3dS > 25000; BC.RainMax3dS(ind) = NaN; 
ind = BC.Sediment >30000; BC.Sediment(ind) = NaN; 
ind = BC.Imperv > 30000; BC.Imperv(ind)= NaN;
%transformations
BC.OSDS_Eff = sqrt(BC.OSDS_Eff); 
BC.Nuts = sqrt(BC.Nuts);
BC.Sediment = sqrt(BC.Sediment); 
BC.HPop15km =  sqrt(BC.HPop15km); 
%FISH
ind = BC.scrapers > 35; BC.scrapers(ind) = NaN; 
ind = BC.herbivores> 110; BC.herbivores(ind) = NaN;
ind = BC.total_biomass> 200; BC.total_biomass(ind) = NaN;
BC.total_biomass = sqrt(BC.total_biomass); 
%% Plot Correlation Matrix (remove variables that are not predictors). 
BCplot = BC; 
BCplot= removevars(BCplot,{'TCEnd'});%remove response variable 
plot_corr_matrix(BCplot)
%Remove variables after assessing correlation.
BCmodel =  BC; 
BCmodel= removevars(BCmodel,{'herbivores','HPop15km','SST_STD','SST_MEAN','CHL','RainAnnSum','browsers'});

%% CALCULATE PERCENTILES AND SAVE 
% Establishsed Percentiles and Categorical variable. 
thresh = round(prctile(BCmodel.TCEnd,[25 75]));
indlow = BCmodel.TCEnd <= thresh(1,1); BCmodel.CAT(indlow) = 3;
indmod = BCmodel.TCEnd > thresh(1,1) & BCmodel.TCEnd < thresh(1,2);BCmodel.CAT(indmod) = 2;
indhigh = BCmodel.TCEnd >=thresh(1,2); BCmodel.CAT(indhigh) = 1;
%Convert CAT to a categorical 
BCmodel = convertvars(BCmodel,{'CAT'},'categorical'); OLRdata = BCmodel;
cd '/Users/jgove/Documents/jamisongove/EOD/KONA_IEA/PaperI/data/Analysis/Final_Data/Revision'
writetable(BCmodel,'OLR_REEF.BUILDER.COVER.csv')


%% FIGURE 4A
figure; hold on
[fi,xi] = ksdensity(BC.TCEnd,'bandwidth',6);
fi = fi*10; 
xi2 = [min(xi):0.05:max(xi)]; 
fi2 = interp1(xi,fi,xi2); 
hold on
map = brewermap(15,'RdPu');
indL = find(xi2 <= thresh(1,1)); 
indM = find(xi2 >= thresh(1,1) & xi2 <= thresh(1,2));
indH = find(xi2 >= thresh(1,2)); 
clrs = [map(4,:);map(8,:);map(12,:)];
ha1 = area(xi2(indL),fi2(indL),'FaceColor',clrs(1,:)); ha1.FaceAlpha = 0.5; 
%ha2 = area(dx1,PH(ind(2),:),'FaceColor',[1 1 1]);
ha2 = area(xi2(indM),fi2(indM),'FaceColor',clrs(2,:)); ha2.FaceAlpha = 0.5; 
ha3 = area(xi2(indH),fi2(indH),'FaceColor',clrs(3,:)); ha3.FaceAlpha = 0.5; 
ha1.EdgeColor = [1 1 1]; ha2.EdgeColor = [1 1 1]; ha3.EdgeColor = [1 1 1]
p= plot(xi,fi,'linewidth',3,'color',[0.7 0.7 0.7]); xlim([0 53]); ylim([0.05 (max(fi)+0.01)])
xline(thresh(1),'--k','linewidth',1.5)
xline(thresh(2),'--k','linewidth',1.5)



