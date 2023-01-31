function [] = E3_Beta_Comparison_bootstrp(badSubs,slctd_strct,append,n_BS)

% baSubs = code of bad subjects ([23,28] for now)
% slctd_strct = which structures to do comarison on ([1,2,4]; %Th, CN and GP: selected based
%                                                    on best model on all conditions combined)
% append = load condition appended (1) or separately (2) ->only works for separated for now
% n_BS = number of bootstrapping iterations (n_BS = 10000;)

saveFolderMat = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices\'; %Windows

%LV variables
load([saveFolderMat 'LV_all.mat'])
LV(badSubs,:) = [];
reg_LV = LV;
structs = {'Th', 'CN', 'Pu', 'GP', 'Hpc', 'Amg', 'Acb'};

%HLM variables
if append == 1
    load([saveFolderMat 'MI_all_dt_sym_salient.mat']) %load conditions 3 & 4 (salient distractor)
    HLM = modulationIdx(:,4);
    predVar = HLM;
elseif append == 2
    load([saveFolderMat 'MI_all_dt_sym_load.mat']) %separated load conditions
    HLM(:,3:4)  = modulationIdx(:,4,3:4);
    predVar_all = HLM;
end

BS_predVar_all = predVar_all; %HLM (to be prredicted values) 
BS_reg_LV_all = reg_LV(:,slctd_strct); %LV values
N=size(BS_predVar_all,1);
% dist_beta=nan(n_BS+1,3,4); mdl_pVal = nan(n_BS+1,1,4); %one p-value for each model 

for n = 1:n_BS
    for ld = 3:4
    bsSelVec = randi(N,[N,1]);
    BS_slctd_HLM = BS_predVar_all(bsSelVec,:);
    BS_slctd_LV = BS_reg_LV_all(bsSelVec,:);
    
    mdldef = fitlm(BS_slctd_LV,BS_slctd_HLM(:,ld),'PredictorVars',{structs{slctd_strct}});
    dist_beta(n,:,ld) = mdldef.Coefficients.Estimate(2:end,:)';
    mdl_pVal(n,1,ld) = coefTest(mdldef);
    end
end
%add beta estimates of the main model (all subs) to the end of dist_beta
for ld = 3:4
    mdldef = fitlm(BS_reg_LV_all,BS_predVar_all(:,ld),'PredictorVars',{structs{slctd_strct}});
    dist_beta(n_BS+1,:,ld) = mdldef.Coefficients.Estimate(2:end,:)';    
    mdl_pVal(n_BS+1,1,ld) = coefTest(mdldef);
end
save([saveFolderMat filesep 'Beta_load_sal_pVal_BS_2'],'dist_beta','mdl_pVal')
