%% Script  modified with permission from Aryan Safaie
%Originally published in Safaie et al., 2018 Nature Communications 
%(DOI:10.1038/s41467-018-04074-2)
%Jamison Gove NOV 2022

%% LOAD DATA FOR OLR
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 14);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Depth", "GearRegRank", "TCEnd", "total_biomass", "grazers", "scrapers", "OSDS_Eff", "Imperv", "Sediment", "RainMax3dS", "WPow975pct", "PAR", "Nuts", "CAT"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
OLRdata = readtable("/Users/jgove/Documents/jamisongove/EOD/KONA_IEA/PaperI/data/Analysis/Final_Data/Revision/Submission/OLR_REEF.BUILDER.COVER.csv", opts);
%% RUN OLR ANALYSIS
RunOrdinalRegression(OLRdata); 
%% SECTION FINDS TOP MODEL AND PLOTS PROBABILITIES
%Select Top Models: Determined by AICc 
indtop = find((dAIC<=2)); %find all dAIC <=0

%Create new variables of the top models
[AICtop] = AIC(indtop);
R2top = R2(indtop);
BVtop = BV(indtop,:); 
SEtop = SE(indtop,:);
PVtop = PV(indtop,:); 
VARStop = VARS(indtop,:);

%Clean up workspace
clearvars AIC BV PV SE VARS R2 

%Order Models by AIC from low to high then rearrange all other variables. 
[A i] = sort(AICtop); 
AICtop = AICtop(i); R2top = R2top(i); BVtop = BVtop(i,:);
SEtop = SEtop(i,:); PVtop = PVtop(i,:);  VARStop = VARStop(i,:); 

%Assess models and note which models have significant predictors (p < 0.05).   
sigpred1 = find(PVtop(1,3:end) < 0.05); 
isempty(sigpred1)
    
%If output is 0, then
mi = 1;
%If output is 1, then there are no significant predictors...move on to next
%model. 

% Get Data: We need the actual data to run the scenario analysis
%Get variable choices from the dataset 
choices = OLRdata.Properties.VariableNames;

vars = VARStop(mi,:);%Get the variables that are asscoiated with the scenario model  
ivs = find(ismember(choices,vars(3:end))); %find the data columns - omitting the first 2 as these are the intercepts 
X = OLRdata(:,ivs); %pullout the data 

%Get Betas from scenario model 
Betas = BVtop(mi,:); %Selects predictor Betas

%% Probability Scenario Analysis
    %Inputs 
        %X: matrix of native-unit covariates
        % Betas: MLE estimates of parameters for logit model [intercepts + coeffs]
        % indIV: index column of predictors we're running in the
        % scenario
        % indIV: index column of covariate we want to change w/in X
        % flag: if =1, x1 starts at its mean; if =0, x1 starts at its min
    %OUTPUTS
        % dx: amount, in native-units, change in the desired IV
        % PH, PM, PL: Probability Curves for High, Moderate, and Low
        % ReefBuilder categories. 
        
%ind1 and ind2 represent the two variables for the scenario. In this case, the two variables...
%are Scrapers and Wastewater Pollution; indIV1 = 1; indIV2 = 2; 
%
indIV1 = 1; indIV2 = 2; 
flag = 0; %flag: if =1, x1 starts at its mean; if =0, x1 starts at its min. Suggest 0 
    [dx1, dx3, PH, PM, PL] = deltaProbability_2variables(X{:,:}, Betas, indIV1,indIV2, flag);
    
%% 

%% CONTOUR PLOT
%COMBINED LOW, MODERATE, HIGH; 
figure('Renderer', 'painters', 'units','centimeters','Position', [50 50 12 12])
t = tiledlayout(1,1);
ax1 = axes(t); %start with the first axis
clev = ([0:0.01:1])%range of values to plot
[cs1, hc1] = contourf(ax1,PM,clev); hold on %filled contour for Probability Moderate/Low
set(hc1(:), 'edgecolor', 'none');

map = brewermap(15,'RdPu');%Custome color ramp from Brewer Map
clrs = [map(4,:);map(8,:);map(12,:)]; %These are colors for low, moderate, high (respectively)
cmap = customColormap([map(4,:);1 1 1;map(8,:)],100,'linear');%Create a divergent color ramp
colormap(ax1,cmap)%set colors 
caxis(ax1,[0 1])%set color axis

steps_large = [0:0.1:1];%set probability contours to draw
[c1,h1] = contour(ax1,PM,steps_large,'color',[0.6 0.6 0.6],'linewidth',1);hold on
steps_small = [0.05:0.1:1];%set probability contours to draw
[c11,h11] = contour(ax1,PM,steps_small,'--','color',[0.6 0.6 0.6],'linewidth',0.5);hold on
%draw 0.5 line thicker to show the transition. 
[c50,h50] = contour(ax1,PM,[0.5 0.5],'color',[0.5 0.5 0.5],'linewidth',3);hold on

%Do all the same for High, but on ax2

% Only want to ploy PH that is >=0.5
PHplot = PH; %rename to preserve original data
ind= PH<=0.49; PHplot(ind) = NaN; 
ax2 = axes(t);
clev2 = ([0.5:0.01:1]);
[cs2, hc2] = contourf(ax2,PHplot,clev2); hold on
set(hc2(:), 'edgecolor', 'none');
cmap = customColormap([1 1 1;map(12,:)],50,'linear');%Create a divergent color ramp
cmap = cmap(4:end,:)%We don't want to plot all white on the lowest value, but really light purple
colormap(ax2,cmap); caxis(ax2,[0.5 1]);

steps_large = [0:0.1:1];%set probability contours to draw
steps_small = [0.05:0.1:1];%set probability contours to draw
[c2,h2] = contour(ax2,PHplot,steps_small,'--','color',[0.6 0.6 0.6],'linewidth',0.5);hold on
steps = [0.5:0.05:1];%contour lines to plot
[c22,h22] = contour(ax2,PHplot,steps_large,'color',[0.6 0.6 0.6],'linewidth',1);hold on
%Thick line along the 50:50 line
[c50H,h50H] = contour(ax2,PHplot,[0.5 0.5],'color',[0.5 0.5 0.5],'linewidth',3);hold on

ax2.Visible = 'off';

% link the limits of axes
linkaxes(t.Children)
%Set colorbars and revisit color axis 
cb1 = colorbar(ax1);
cb1.Layout.Tile = 'south';
cb1.Label.String = 'Low - Moderate';
cb1.Ticks = [0:0.1:1];
cb1.TickLabels = {'1','0.9','0.8','0.7','0.6','0.5','0.6','0.7','0.8','0.9','1'}

cb2 = colorbar(ax2);
cb2.Layout.Tile = 'south';
cb2.Label.String = 'High';
cb2.Ticks = [0.5:0.1:1];
cb2.TickLabels = {'0.5','0.6','0.7','0.8','0.9','1'}
%SET X AND Y AXES and LIMITS

fishaxis = round(dx1*10);%converts to kg/ha and rounds to whole number. 
fishind = [0:25:350];  
%Gets dx3 indices closest to these values. 
[val,xticks] = min(abs(fishaxis-fishind)); 
set(ax1,'Xtick',xticks);
set(ax1,'XtickLabel',fishind);
xlabel(ax1,'Scraper Biomass (kg/ha)');

%WW has been square root transformed, so take this into account.
WWind = [0:100:1000]; 
WWvals = WWind.^2; %Get y-labels correct. 
%Gets dx3 indices closest to these values. 
[val,yticks] = min(abs(dx3-WWind)); 
set(ax1,'Ytick',yticks)
WWlabel = WWvals;
set(ax1,'YtickLabel',WWlabel)
ylabel(ax1,'Wastewater Pollution (L/ha)')

set(ax1,'xlim',[1 xticks(end-2)]); %0 value is at 1; sets upper at 325/kg ha
set(ax1,'ylim',[1 max(yticks)]); %0 value is at 1; sets upper at; sets bounds to 1 million L/ha, 

ax = gcf; 
cd '/Users/jgove/Documents/jamisongove/EOD/KONA_IEA/PaperI/figures/Manuscript/SUBMITTED/REVISION'

%% MARKING MANAGMENT SCENARIOS
%STARTING FISH LOW WW HIGH
fish = 30;
ww = sqrt(600000)
[val,fishlow] = min(abs(fishaxis-fish));
[val,wwhigh] = min(abs(dx3-ww))

ScenS_PL = PL(wwhigh,fishlow);
ScenS_PM = PM(wwhigh,fishlow);
ScenS_PH = PH(wwhigh,fishlow);
ScenarioS = [ScenS_PL ScenS_PM ScenS_PH]

%SCENARIO B: HIGH FISH HIGH WW
fish = 250;
ww = sqrt(600000);
[val,fishhigh] = min(abs(fishaxis-fish));
[val,wwhigh] = min(abs(dx3-ww));
ScenA_PL = PL(wwhigh,fishhigh);
ScenA_PM = PM(wwhigh,fishhigh);
ScenA_PH = PH(wwhigh,fishhigh);
ScenarioA = [ScenA_PL ScenA_PM ScenA_PH];
%SCENARIO B: LOW FISH LOW WW
fish = 30;
ww = sqrt(2500);
[val,fishlow] = min(abs(fishaxis-fish));
[val,wwlow] = min(abs(dx3-ww));
ScenB_PL = PL(wwlow,fishlow);
ScenB_PM = PM(wwlow,fishlow);
ScenB_PH = PH(wwlow,fishlow);
ScenarioB = [ScenB_PL ScenB_PM ScenB_PH];

%SCEARIO C: FISH HIGH WW LOW
fish = 250;
ww = sqrt(2500)
[val,fishlow] = min(abs(fishaxis-fish));
[val,wwhigh] = min(abs(dx3-ww));
ScenC_PL = PL(wwhigh,fishlow);
ScenC_PM = PM(wwhigh,fishlow);
ScenC_PH = PH(wwhigh,fishlow);
ScenarioC = [ScenC_PL ScenC_PM ScenC_PH];

AllScenarios = [ScenarioS;ScenarioA;ScenarioB;ScenarioC];

%% ADD IN SCENARIO LOCATRIONS. 
%STARTING FISH LOW WW HIGH
fish = 30;
[val,fishlow] = min(abs(fishaxis-fish));
xline(fishlow,'g');
fish = 250;
[val,fishhigh] = min(abs(fishaxis-fish));
xline(fishhigh,'g');

ww = sqrt(600000);
[val,wwhigh] = min(abs(dx3-ww));
yline(wwhigh,'g');

ww = sqrt(2500);
[val,wwlow] = min(abs(dx3-ww));
yline(wwlow,'g');


 