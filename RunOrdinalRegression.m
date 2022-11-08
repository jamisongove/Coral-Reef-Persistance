%% Compute all OLR permutations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RunOrdinalRegression(BC); 

clc

DepVar = BC.CAT; 
BC = removevars(BC, {'TCEnd'});
BC = removevars(BC, {'CAT'});

% % Define some cell arrays, vectors, and matrices for storage:
%First standardize the variables
BCnorm = normalize(BC{:,:}); %Note that zscores is the default and this handles NaNs. This is the same as Value = [X - mean/std]
V = [BC.Properties.VariableNames];
VN = {}; % cell array of variable names
BV = []; % cell array of logit parameters
devs = []; % vector of sum of deviance residuals for each logit model
Residuals = []; % array of all N residuals from each mnrfit run
AICs = []; % AIC
AICcs = []; % AIC correction
MPR2 = []; % McFadden's Pseudo R2
SE = []; % cell array of standard errors for each parameter for each logit model
SigVars = []; % cell array of names of significant variables for each logit model
Tvals = []; % cell array of students t values for each parameter for each logit model
PV = []; % cell array of p values for each parameter of each logit model
Rd = []; % 2 column vector of the average and sum (= devs) of deviance residuals for each logit model
Brants = []; % Results from Brant's Wald test: Pval(1) Pval(2) h(1) h(2) Wstat(1) Wstat(2) df
LR = []; % vector of results (h,pValue,stat,cValue) from likelihood ratio tests for each non-proportional odds test
GroupNums = []; % Number of models in each group of IVs
numpreds = 4; %number of predicto
modelruns = nchoosek([1:size(BC,2)],numpreds); %Set up to run all possible combination of predictors but limited to 6 predictors to reduce overfitting
%modelruns = ff2n(size(BC,2)); modelruns = logical(modelruns(2:end,:)); 

% Now for each of these combos, run an mnrfit
  for jp = 1:size(modelruns,1)
            %GNcounter = GNcounter+1;
            W = BCnorm(:,(modelruns(jp,:)));  % W = temporary matrix of I for current logit run
            % Define some constants
            n = size(W,1);  % number of observations
            k = numel(unique(DepVar)); % number of probability categories
            p = size(W,2); % number of IV categories
            DF = n*(k - 1) - (k - 1 + p); % degrees of freedom
            VN(jp,:) = [{'C1'},{'C2'}, V(modelruns(jp,:))]; %Variables in the model 
        
            [B dev stats] = mnrfit(W,DepVar, 'model','ordinal','link','logit','interactions','off');
            % also run the Null (intercept only) model
            [Bnull devnull statsnull] = mnrfit([],DepVar, 'model','ordinal','link','logit','interactions','off');
            BV(jp,:) = B';
            rss1 = sum(stats.residd);
            rss2 = sum((stats.residd).^2);
            NumParam = numel(B);
            aictemp = 2*NumParam + dev; % via McCullagh and Nelder 1990
            devs = [devs; dev];
            Residuals{end+1} = stats.residd;
            AICs = [AICs; aictemp];
            aiccorrection = (2*NumParam*(NumParam+1)) /(n - NumParam - 1);
            AICcs = [AICcs; (aictemp + aiccorrection)];
            % And compute McFadden's Pseudo R2
            MPR2 = [MPR2; 1 - (dev/devnull)];           
            SE(jp,:) = stats.se';    
            Tvals(jp,:) = stats.t';
            PV(jp,:) = stats.p';
            Rd(jp,:) = [nanmean(stats.residd) nansum(stats.residd)];        
  end
% 3. SUMMED EFFCTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AIC = AICcs; % Use weighted AIC
[AIC0] = min(AIC); % find min AIC
dAIC = AIC-AIC0;
assignin('base','AIC',AIC); assignin('base','dAIC',dAIC)

indtop = find((dAIC<=2)); %find all dAIC <=0
%Reduce all models to top and get rid of intercetps (i.e. first 2 indicies)
k = numel(unique(DepVar))-1;%this is the number of intercept in the model output, which we want to get rid of. 

R2top = MPR2(indtop); assignin('base','R2',MPR2)
BVtop = BV(indtop,k+1:end); assignin('base','BV',BV)
SEtop = SE(indtop,k+1:end); assignin('base','SE',SE)
PVtop = PV(indtop,k+1:end); assignin('base','PV',PV)
VARStop = VN(indtop,k+1:end); assignin('base','VARS',VN)
% 
 
