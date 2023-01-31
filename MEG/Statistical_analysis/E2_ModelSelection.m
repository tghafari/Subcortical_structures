function [] = E2_ModelSelection(predictor_data,badSubs)

% MODEL SELECTION STRATEGY
% the idea is that you have to find the model which is associated with the
% smallest AIC and BIC.

% predictor_data = which data to be predicted from LVs: 1.1 -> alpha (all conditions) 1.2 -> alpha(load conditions separated)
%                                                       1.3 -> alpha(conditions 3 & 4) 1.4 -> HLM difference between load
%                                                       2.1-> RFT(group) 2.2=> RFT(individual) 
%                                                       3.1-> behavioral(all conditions) 3.2-> behavioral(conditions separately)
% baSubs = code of bad subjects ([23,28] for now)
predictor_data = 1.1;
badSubs = [23,28];
%% Load MI and LV data
% saveFolderMat = '/Volumes/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/group_level/Lateralization_indices/'; %Mac
saveFolderMat = 'Z:\Programming\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices\'; %Windows
load([saveFolderMat 'LV_all.mat'])
LV(badSubs,:) = [];
X_LV = LV;
structs = {'Thalamus', 'Caudate n.', 'Putamen', 'Globus Pallidus', 'Hippocampus', 'Amygdala', 'n. Accumbens'};

if predictor_data == 1.1
    load([saveFolderMat 'MI_all_dt_sym']) %alpha
    HLM = modulationIdx(:,4); 
    predVar = HLM;
elseif predictor_data == 1.2
    load([saveFolderMat 'MI_all_dt_sym_load.mat']) %separated load conditions
    HLM(:,1:4)  = modulationIdx(:,4,:); 
    predVar_all = HLM;
elseif predictor_data == 1.3
    load([saveFolderMat 'MI_all_dt_sym_salient.mat']) %load conditions 3 & 4 (salient distractor)
    HLM = modulationIdx(:,4); 
    predVar = HLM;
elseif predictor_data == 1.4
    load([saveFolderMat 'HLM_diff.mat']) %Low load - high load in salient distractor
    HLM = HLM_diff_load(:,2);  
    predVar = HLM;
elseif predictor_data == 2.1
    load([saveFolderMat 'RFT_MI_all_dt.mat']) %group roi
    HLM  = modulationIdx(:,4); 
    predVar = HLM;
elseif predictor_data == 2.2
    load([saveFolderMat 'RFT_MI_indiv_dt']); %indiv roi
    HLM  = modulationIdx(:,4); 
    predVar = HLM;
elseif predictor_data == 3.1
    load([saveFolderMat 'BAP_RT_acc_IES_allsubs.mat']) %behavioral indices
    RT = BAP_RT_acc(:,1); acc = BAP_RT_acc(:,2); IES = BAP_RT_acc(:,3);
    predVar = input('which variable do you want to be predicted for LV in the GLM? (RT/acc/IES)'); 
elseif predictor_data == 3.2
    load([saveFolderMat 'BAP_RT_acc_IES_allsubs_load.mat']) %behavioral indices
    RT = BAP_RT_acc(:,:,1); acc = BAP_RT_acc(:,:,2); IES = BAP_RT_acc(:,:,3);
    predVar_all = input('which variable do you want to be predicted for LV in the GLM? (RT/acc/IES)'); 
end

%% Create all unique combinations  of 'n' numbers, from 1 to n
%this so you can just try all combinations of regressors to include in the model and
%for all possible number of regressors (meaning, 3,4,5,6, etc).
nReg = cell(7,1);

for n_regr = 1:7
nReg{n_regr} = nchoosek(1:7,n_regr); % contains all combinations of n numbers from 1:7
end
%% n regr models
% % if  predictor_data == 1.2 -> should be done manually
beta_est = cell(1,2); %only needed for load separate
mdl_pVal = cell(1,2); %model p-values
mdl_fVal = cell(1,2); %model's f values
% for ld = 1:4 %1:4 %comment this "for loop" if load conditions are appended
%     predVar = predVar_all(:,ld);
%     predVar = predVar_all(:,setxor(1:4,ld)); %jacknife method
%     predVar = mean(predVar,2); %average over HLM in three remaining conditions
% Preallocation
LME = cell(7,13);
str_lbl = cell(7,13);
AIC = zeros(7,13);
BIC = zeros(7,13);
logLikelihood = zeros(7,13);
Rsqrd_adjst = zeros(7,13);

formula{1} = 'predVar ~ 1 + LV';
formula{2} = 'predVar ~ 1 + LV_1 + LV_2'; 
formula{3} = 'predVar ~ 1 + LV_1 + LV_2 + LV_3';
formula{4} = 'predVar ~ 1 + LV_1 + LV_2 + LV_3 + LV_4';
formula{5} = 'predVar ~ 1 + LV_1 + LV_2 + LV_3 + LV_4 + LV_5';
formula{6} = 'predVar ~ 1 + LV_1 + LV_2 + LV_3 + LV_4 + LV_5 + LV_6';
formula{7} = 'predVar ~ 1 + LV_1 + LV_2 + LV_3 + LV_4 + LV_5 + LV_6 + LV_7';

for n_regr = 1:7
%     2 %Th and CN
%`    7   %all structures
%     1:7
%     3 %Th, CN and GP
    for j = 1:size(nReg{n_regr,:},1)
%         1 %Th and CN
%         1:size(nReg{n_regr,:},1) % make all structure and all regressors
%         2 %Th, CN and GP
      
        designtmp = table(predVar, X_LV(:,nReg{n_regr}(j,1:n_regr)),'VariableNames',{'predVar','LV'});
        designtbl = splitvars(designtmp);
        
        LME{n_regr,j} = fitlme(designtbl, formula{n_regr} ,'DummyVarCoding','effects', 'CheckHessian',true, 'FitMethod','REML');
        str_lbl{n_regr,j} = {structs{nReg{n_regr}(j,1:n_regr)}};
        AIC(n_regr,j) = LME{n_regr,j}.ModelCriterion.AIC;
        BIC(n_regr,j) = LME{n_regr,j}.ModelCriterion.BIC;
        logLikelihood(n_regr,j) = LME{n_regr,j}.LogLikelihood;
        Rsqrd_adjst(n_regr,j) = LME{n_regr,j}.Rsquared.Adjusted;
    end
end

AIC(AIC==0)=nan;  %NaN all the zeros, to make sure they don't interfer with calculations
BIC(BIC==0)=nan;  
logLikelihood(logLikelihood==0)=nan;  
Rsqrd_adjst(Rsqrd_adjst==0)=nan;  

%% find smallest coefficients for each condition (best AIC or best BIC), eventually could check best average both...
best_AICs = zeros(7,2); %nr. of row => nr. of regressors, column 1 => AIC, column 2 => idx
best_BICs = zeros(7,2); %nr. of row => nr. of regressors, column 1 => BIC, column 2 => idx
AIC_best_lbl = cell(7,7);
BIC_best_lbl = cell(7,7);

for n_regr = 1:7
%     2 %Th and CN
%     1:7
%     3 %Th, CNd an GP
    best_AICs(n_regr,1) = min(AIC(n_regr,:));
    best_AICs(n_regr,2) = find(AIC(n_regr,:)==min(AIC(n_regr,:)));
    best_BICs(n_regr,1) = min(BIC(n_regr,:));
    best_BICs(n_regr,2) = find(BIC(n_regr,:)==min(BIC(n_regr,:)));
    AIC_best_lbl(n_regr,1:n_regr) = str_lbl{n_regr,AIC(n_regr,:)==min(AIC(n_regr,:))}; %LABELS ASSOCIATED WITH THE WINNING MODEL CONTAINING n regressors
    BIC_best_lbl(n_regr,1:n_regr) = str_lbl{n_regr,BIC(n_regr,:)==min(BIC(n_regr,:))};
end

%best model regarding AIC/BIC
avg_AICs_BICs   = mean([best_AICs(:,1),best_BICs(:,1)],2);
AB_best_n_reg(1,1) = find(best_AICs(:,1)==min(best_AICs(:,1))); %smallest AIC
AB_best_n_reg(1,2) = find(best_BICs(:,1)==min(best_BICs(:,1))); %smallest BIC
AB_best_n_reg(1,3) = find(avg_AICs_BICs==min(avg_AICs_BICs));   %smallest in the average of AIC and BIC

%% FIND MODEL WITH LOWEST AIC, HIGHEST R2, HIGHEST LOGLIKELIHOOD
% another thing you can do is to use not only AIC and BIC but also Log
% likelihood and R squared as criteria to decide which one is the best
% model, see below for examples on how to extract them

% Choose the model with the highest Rsqrd & loglikelihood
best_LLs = zeros(7,2); %nr. of row => nr. of regressors, column 1 => AIC, column 2 => idx
best_Rsqrds = zeros(7,2); %nr. of row => nr. of regressors, column 1 => BIC, column 2 => idx
LL_best_lbl = cell(7,7);
Rsqrd_best_lbl = cell(7,7);

for n_regr = 1:7
%     2 %Th and CN
%     1:7
%     3 %Th, CN and GP
    best_LLs(n_regr,1) = max(logLikelihood(n_regr,:));
    best_LLs(n_regr,2) = find(logLikelihood(n_regr,:)==max(logLikelihood(n_regr,:)));
    best_Rsqrds(n_regr,1) = max(Rsqrd_adjst(n_regr,:));
    best_Rsqrds(n_regr,2) = find(Rsqrd_adjst(n_regr,:)==max(Rsqrd_adjst(n_regr,:)));
    LL_best_lbl(n_regr,1:n_regr) = str_lbl{n_regr,logLikelihood(n_regr,:)==max(logLikelihood(n_regr,:))}; %LABELS ASSOCIATED WITH THE WINNING MODEL CONTAINING n regressors
    Rsqrd_best_lbl(n_regr,1:n_regr) = str_lbl{n_regr,Rsqrd_adjst(n_regr,:)==max(Rsqrd_adjst(n_regr,:))};
end

%best model regarding LL/Rsqrd
avg_LLs_Rsqrds  = mean([best_LLs(:,1),best_Rsqrds(:,1)],2);
LR_best_n_reg(1,1) = find(best_LLs(:,1)==max(best_LLs(:,1))); %greatest LR
LR_best_n_reg(1,2) = find(best_Rsqrds(:,1)==max(best_Rsqrds(:,1))); %greatest R-squared
LR_best_n_reg(1,3) = find(avg_LLs_Rsqrds==max(avg_LLs_Rsqrds));   %smallest in the average of AIC and BIC

%% Find the loglikelihood and R_squared associated with the best BIC/AIC model - not of interest atm
% for n_regr=1:7
%     best_LLs(n_regr,1) = LME{n_regr,best_AICs(n_regr,2)}.ModelCriterion.LogLikelihood;
%     best_Rsqrds(n_regr,:) = cell2mat(struct2cell(LME{n_regr,best_AICs(n_regr,2)}.Rsquared))';
% end

%% Choose which best model to use (AIC/BIC or LL/Rsqrd)

sprintf('best AIC/BIC = %d',AB_best_n_reg(1,3))
sprintf('best LL/R_squared = %d',LR_best_n_reg(1,3))

AB = AB_best_n_reg;
LR = LR_best_n_reg;
best_n_reg = input('AB or LR?');
% AB_best_n_reg; %Th, CNd an GP
if best_n_reg == AB
    best_model_selected = best_AICs;
    best_lbl = AIC_best_lbl;
elseif best_n_reg == LR
    best_model_selected = best_Rsqrds;
    best_lbl = Rsqrd_best_lbl;
end


%% consider the model associated with the best coefficients among all of them
%and fit linear regression model (instead of fitlme) for plotting purposes 

%mdldef is chosen model - in your case with three regressors 
%define which structures showed the best performance
mdldef =fitlm(X_LV(:, nReg{best_n_reg(1,3)}(best_model_selected(best_n_reg(1,3),2),:)),predVar,...
    'PredictorVars',best_lbl(best_n_reg(1,3),1:best_n_reg(1,3)));
if predictor_data == 1.2 || predictor_data == 3.2
    beta_est{ld} = mdldef.Coefficients;
    [mdl_pVal{ld},mdl_fVal{ld}] = coefTest(mdldef);
else
    beta_est = mdldef.Coefficients;
    [mdl_pVal,mdl_fVal,df] = coefTest(mdldef);
end

%% Plot all the Beta estimates of the best model

colmp = bone;
fsz   = 40;
num_best_struct = nReg{best_n_reg(1,3)}(best_model_selected(best_n_reg(1,3),2),:);
xvars = X_LV(:,num_best_struct);

figure();clf;

for n_est = 1:size(xvars,2)
    bar(n_est,mdldef.Coefficients.Estimate(n_est+1),'facecolor',colmp(n_est*50,:));
    hold all
    errorbar(n_est,mdldef.Coefficients.Estimate(n_est+1),mdldef.Coefficients.SE(n_est+1),'k', 'LineWidth', 3)
    if predictor_data == 1.2 || predictor_data == 3.2
    title(['Load Condition ' num2str(ld)],'FontSize',4);
    end
end
set(gca,'XTick',1:length(num_best_struct),'FontSize', fsz,'XTickLabel',...
    best_lbl(best_n_reg(1,3),1:best_n_reg(1,3)), 'XTickLabelRotation', 45)
ylim([-3.5, 1.5])
ylabel('Beta Coefficient', 'FontSize', fsz)

txt = sprintf('R_{adjusted}= %.3f', mdldef.Rsquared.Adjusted);
text_x = min(xlim) + .5; text_y = max(ylim)-.5;
tx = text(text_x, text_y, txt);
tx.FontSize = 15;
tx.FontWeight = 'bold';

if predictor_data == 1.2 || predictor_data == 3.2
    txt_pVal = ['\itp\rm\bf', sprintf(' = %.3f', mdl_pVal{ld})];
else
    txt_pVal = ['\itp\rm\bf', sprintf(' = %.3f', mdl_pVal)];
end
t = text(text_x, text_y-.2, txt_pVal);
t.FontSize = 15;
t.FontWeight = 'bold';
%% added variable plots 
%this is basically the partial regression, so the correlation between Th. CN and GP and HLM, 
% for example, after correcting for the influence of the other variables in the model (so in your case Th, CN)
for plts = 1:best_n_reg(1,3)
figure();clf; 
plotAdded(mdldef, plts+1); %GP (because 4 is GP's index in the model (the first is the intercept always))
txt = num2str(mdldef.Coefficients{plts+1,4}, 'p-value= %.3f');
text_x = min(xlim); text_y = min(ylim);
text(text_x, text_y, txt, 'FontSize', 10);
if predictor_data == 1.2 || predictor_data == 3.2
title(['Load Condition ' num2str(ld)],'FontSize',4);
end
end

end %comment this for append == 1 and salient
% save([saveFolderMat filesep 'Beta_estimates_alpha_sal'],'beta_est')

