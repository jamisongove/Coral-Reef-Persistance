function [dx1, dx3, PH, PM, PL] = deltaProbability_2variables(X, Betas, indIV1,indIV2, flag)
% function to calculate the change in probability of bleaching categories
% based on a dx change in independent variable indexed by indIV
% Inputs:

    % X: matrix of native-unit covariates
    % Betas: MLE estimates of parameters for logit model [intercepts + coeffs]
    % dx: amount, in native-units, change in the desired IV
    % indIV: index column of covariate we wanna change w/in X
    % flag: if =1, x1 starts at its mean; if =0, x1 starts at its min

    % Aryan Safaie, 04/27/2017
    % Modified for J.Gove
% Parse out intercepts.
C1 = Betas(1);
C2 = Betas(2);

% Which covariate are we changing?
Bs = Betas(3:end);%The first two are interecepts; this will need to change based on number of ordinal inputs
Beta = Bs(indIV1);
mu = nanmean(X(:,indIV1));
sig = nanstd(X(:,indIV1));
if flag == 0
    x1 = min(X(:,indIV1)); % base minimum value of covariate
    z1 = (x1 - mu)/(sig);
    z2 = z1 + linspace(0,6,300); % %Customize based on disrubtion of data and extent
elseif flag == 1
    x1 = nanmean(X(:,indIV1));
    z1 = (x1 - mu)/(sig);
    z2 = z1 + linspace(-6,6,100); % z2 corresponds to a -6:6 standardized unit range
end
% dx = x2 - x1, dz = z2 - z1
dx1 = sig*(z2-z1);
dx1 = dx1(:);
dz1 = z2 - z1;
dz1 = dz1(:);

%Now for 2nd Driver
Beta = Bs(indIV2);
mu = nanmean(X(:,indIV2));
sig = nanstd(X(:,indIV2));
if flag == 0
    x3 = min(X(:,indIV2)); % base minimum value of covariate
    z3 = (x3 - mu)/(sig);
    z4 = z3 + linspace(0,6,300); %%Customize based on disrubtion of data and extent
elseif flag == 1
    x3 = nanmean(X(:,indIV2));
    z3 = (x3 - mu)/(sig);
    z4 = z3 + linspace(-6,6,100); % z2 corresponds to a -6:6 standardized unit range
end
% dx = x2 - x1, dz = z2 - z1
dx3 = sig*(z4-z3);
dx3 = dx3(:);
dz3 = z4 - z3;
dz3 = dz3(:);

% Which covariates are we NOT changing?
betaL  = Bs; betaL([indIV1 indIV2]) = [];
XL = X; XL(:,[indIV1 indIV2]) = [];
mus = nanmean(XL);
sigs = nanstd(XL);

    zL = (mus - mus)./sigs; % ie keep the other variables at their mean vals, zL should be zero

BLeft = sum(zL(:).*betaL(:));

%%%% CALCUALTE HIGH, MODERATE, LOW PROBABILITIES  %%%%%%%%
%Note that for the below, it holds each row in z4 and calculates all
%column values (i.e. all z2) over each row in z4. So the y axis will be z4
%and x axis will be z2
C = C1;
PH = ones(numel(z2),numel(z4))*NaN; 
for i = 1:numel(z2)
  for j = 1:numel(z4)
den1 = 1 + (exp(-C)) * exp(-Bs(indIV1)*z1) * exp(-Bs(indIV2)*z3) * (exp(-BLeft));
den2 = 1 + (exp(-C)) * exp(-Bs(indIV2)*z4(i)) * exp(-Bs(indIV1)*z2(j)) * (exp(-BLeft));
dP1 = (1./(den2)) - (1./(den1));   
P1 = 1./den1 + dP1; 
PH(i,j) = P1; 
  end
end

C = C2;
PM = ones(numel(z2),numel(z4))*NaN; 
for i = 1:numel(z2)
  for j = 1:numel(z4)
den1 = 1 + (exp(-C)) * exp(-Bs(indIV1)*z1) * exp(-Bs(indIV2)*z3) * (exp(-BLeft));
den2 = 1 + (exp(-C)) * exp(-Bs(indIV2)*z4(i)) * exp(-Bs(indIV1)*z2(j)) * (exp(-BLeft));
dP2 = (1./(den2)) - (1./(den1));   
P2 = 1./den1 + dP2; 
PM(i,j) = P2; 
  end
end

PL = 1 - PM; 

end

