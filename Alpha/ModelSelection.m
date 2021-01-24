%% MODEL SELECTION STRATEGY

% the idea is that you have to find the model which is associated with the
% smallest AIC and BIC.

clc;clear

%% Load MI and LV data
saveFolderMat = '/Volumes/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/group_level/Lateralization_indices/'; %Mac
load([saveFolderMat 'MI_all_dt.mat'])
load([saveFolderMat 'LV_all.mat'])

LV([23,28],:) = [];
X_LV = LV;
HLM  = modulationIdx(:,4); 
structs = {'Th', 'CN', 'Pu', 'GP', 'Hpc', 'Amg', 'Acb'};

%% reate all unique combinations  of 'n' numbers, from 1 to n
%this so you can just try all combinations of regressors to include in the model and
%for all possible number of regressors (meaning, 3,4,5,6, etc).
nReg = cell(7,1);

for n_regr = 1:7
nReg{n_regr} = unique(sort(combnk(1:7,n_regr)), 'rows', 'stable'); % contains all unique combinations of n numbers from 1:7
end
%% 3 regr models
% Preallocation
LME = cell(7,13);
str_lbl = cell(7,13);
AIC = zeros(7,13);
BIC = zeros(7,13);

formula{1} = 'HLM ~ 1 + LV';
formula{2} = 'HLM ~ 1 + LV_1 + LV_2'; 
formula{3} = 'HLM ~ 1 + LV_1 + LV_2 + LV_3';
formula{4} = 'HLM ~ 1 + LV_1 + LV_2 + LV_3 + LV_4';
formula{5} = 'HLM ~ 1 + LV_1 + LV_2 + LV_3 + LV_4 + LV_5';
formula{6} = 'HLM ~ 1 + LV_1 + LV_2 + LV_3 + LV_4 + LV_5 + LV_6';
formula{7} = 'HLM ~ 1 + LV_1 + LV_2 + LV_3 + LV_4 + LV_5 + LV_6 + LV_7';

for n_regr = 1:7
    for j = 1:size(nReg{n_regr,:},1)
        
        designtmp = table(HLM, X_LV(:,nReg{n_regr}(j,1:n_regr)),'VariableNames',{'HLM','LV'});
        designtbl = splitvars(designtmp);
        
        LME{n_regr,j} = fitlme(designtbl, formula{n_regr} ,'DummyVarCoding','effects', 'CheckHessian',true, 'FitMethod','REML');
        str_lbl{n_regr,j} = {structs{nReg{n_regr}(j,1:n_regr)}};
        AIC(n_regr,j) = LME{n_regr,j}.ModelCriterion.AIC;
        BIC(n_regr,j) = LME{n_regr,j}.ModelCriterion.BIC;
    end
end

AIC(AIC==0)=nan;  %NaN all the zeros, to make sure they don't interfer with calculations
BIC(BIC==0)=nan;  

%% find smallest coefficients for each condition (best AIC or best BIC), eventually could check best average both...
best_AICs = zeros(7,2); %nr. of row => nr. of regressors, column 1 => AIC, column 2 => idx
best_BICs = zeros(7,2); %nr. of row => nr. of regressors, column 1 => BIC, column 2 => idx
AIC_best_lbl = cell(7,7);
BIC_best_lbl = cell(7,7);

for n_regr = 1:7
    best_AICs(n_regr,1) = min(AIC(n_regr,:));
    best_AICs(n_regr,2) = find(AIC(n_regr,:)==min(AIC(n_regr,:)));
    best_BICs(n_regr,1) = min(BIC(n_regr,:));
    best_BICs(n_regr,2) = find(BIC(n_regr,:)==min(BIC(n_regr,:)));
    AIC_best_lbl(n_regr,1:n_regr) = str_lbl{n_regr,AIC(n_regr,:)==min(AIC(n_regr,:))}; %LABELS ASSOCIATED WITH THE WINNING MODEL CONTAINING n regressors
    BIC_best_lbl(n_regr,1:n_regr) = str_lbl{n_regr,BIC(n_regr,:)==min(BIC(n_regr,:))};
end

%best model regarding AIC/BIC
avg_AICs_BICs   = mean([best_AICs(:,1),best_BICs(:,1)],2);
best_n_reg(1,1) = find(best_AICs(:,1)==min(best_AICs(:,1))); %smallest AIC
best_n_reg(1,2) = find(best_BICs(:,1)==min(best_BICs(:,1))); %smallest BIC
best_n_reg(1,3) = find(avg_AICs_BICs==min(avg_AICs_BICs));   %smallest in the average of AIC and BIC

%% FIND MODEL WITH LOWEST AIC, HIGHEST R2, HIGHEST LOGLIKELIHOOD
% another thing you can do is to use not only AIC and BIC but also Log
% likelihood and R squared as criteria to decide which one is the best
% model, see below for examples on how to extract them

%Put the best loglikelihoods and Rsquared in a matrix
LogLikelihoods = zeros(7,1);
Rsquared       = zeros(7,2); % columns => ordinary and adjusted R-square

for n_regr=1:7
    %find the labels (name of sub structures) corresponding to best model with n regressors
    LogLikelihoods(n_regr,1) = LME{n_regr,best_AICs(n_regr,2)}.ModelCriterion.LogLikelihood;
    Rsquared(n_regr,:)       = cell2mat(struct2cell(LME{n_regr,best_AICs(n_regr,2)}.Rsquared))';
end
%% consider the model associated with the best coefficients among all of them
%and fit linear regression model (instead of fitlme) for plotting purposes 

%mdldef is chosen model - in your case with three regressors 
%define which structures showed the best performance
mdldef =fitlm(X_LV(:, nReg{best_n_reg(1,3)}(best_AICs(best_n_reg(1,3),2),:)),HLM,...
    'PredictorVars',AIC_best_lbl(best_n_reg(1,3),1:best_n_reg(1,3)));

%% Plot all the Beta estimates of the best model

colmp = bone(10);
fsz   = 20;
xvars = X_LV(:, nReg{best_n_reg(1,3)}(best_AICs(best_n_reg(1,3),2),:));

figure();clf;

for n_est = 1:size(xvars,2)
    bar(n_est,mdldef.Coefficients.Estimate(n_est+1),'facecolor',colmp(n_est*3,:));
    hold all
    errorbar(n_est,mdldef.Coefficients.Estimate(n_est+1),mdldef.Coefficients.SE(n_est+1),'k')
%     errorbar(n,mdl.Coefficients.Estimate(n+1),CI(n+1,1),CI(n+1,2),'k'); title('Predictors Betas');
    set(gca,'XTick',1:3,'FontSize', fsz,'XTickLabel', structs, 'XTickLabelRotation', 45)
%  ylim([-3 1]);
end

txt = num2str(mdldef.Rsquared.Adjusted, 'R(adjusted)= %.3f');
text(.5, 1, txt, 'FontSize', 10);

ylabel('Beta Coefficient', 'FontSize', fsz)


%% added variable plots 
%this is basically the partial regression, so the correlation between Th. CN and GP and HLM, 
% for example, after correcting for the influence of the other variables in the model (so in your case Th, CN)
for plts = 1:best_n_reg(1,3)
figure();clf; 
plotAdded(mdldef, plts+1); %GP (because 4 is GP's index in the model (the first is the intercept always))
end

