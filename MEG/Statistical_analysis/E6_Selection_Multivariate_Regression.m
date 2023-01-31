% function [RegressionSummary] = multivariateregression(Y,X)
%Multivariate multiple regression on DVar Y along dimensions d and IVar X along dimensions p. 
%Y is the dependent variable, each row subject, each column is condition. Xmat is matrix containing IV, should be all n x p, where n is subject
%and every column p is different predictor
clear
%Windows
saveFolderMat = 'Z:\Programming\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices\'; 
load([saveFolderMat 'LV_all.mat'])
load([saveFolderMat '/MI_all_dt_sym_load.mat'])

% %Mac
% load('/Users/Tara/Documents/MATLAB/MATLAB-Programs/CHBH-Programs-move to RDS/Perceptual-Load/FieldTrip-Just checking/MI_all_dt_sym_load.mat')
% load('/Users/Tara/Documents/MATLAB/MATLAB-Programs/CHBH-Programs-move to RDS/Perceptual-Load/FieldTrip-Just checking/LV_all.mat')

structs = {'Th', 'CN', 'Pu', 'GP', 'Hpc', 'Amg', 'Acb'};


%% multivariate regression

% Define dependent (HLMs) and independent variables (LV)
Y = squeeze(modulationIdx(:,4,:)); % to get the HLM values of 4 conditions separetely
X = LV(setxor(1:35,[23,28]),:); 

% Set some basic variables
[n,d] = size(Y);
Xcell = cell(1,n);

% Create all unique combinations  of 'n' numbers, from 1 to n
nReg = cell(7,1);

for n_regr = 1:7
nReg{n_regr} = unique(sort(combnk(1:7,n_regr)), 'rows', 'stable'); % contains all unique combinations of n numbers from 1:7
end

% Run the model for all possible combinations
for n_regr = 1:7
    for j = 1:size(nReg{n_regr,:},1)
        Xmat  = [ones(n,1) X(:,nReg{n_regr}(j,1:n_regr))]; %add intercept
        for i = 1:n %prepare the dependent variable for mvregress
            Xcell{i} = [kron([Xmat(i,:)],eye(d))];
        end
        [beta{n_regr,j},Sigma{n_regr,j},E{n_regr,j},CovB{n_regr,j},logL(n_regr,j)] = mvregress(Xcell,Y); %run the model
        str_lbl{n_regr,j} = {structs{nReg{n_regr}(j,1:n_regr)}};
    end
end

% Find the model assiciated with greatest loglikelihood
for n_regr = 1:7
    best_LLs(n_regr,1) = max(logL(n_regr,:));
    best_LLs(n_regr,2) = find(logL(n_regr,:)==max(logL(n_regr,:)));
    LL_best_lbl(n_regr,1:n_regr) = str_lbl{n_regr,logL(n_regr,:)==max(logL(n_regr,:))}; %LABELS 
end

%best model regarding log likelihood
LL_best_n_reg(:,1) = find(best_LLs(:,1)==max(best_LLs(:,1))); %greatest LL
LL_best_n_reg(:,2) = best_LLs(LL_best_n_reg(:,1),2);
sprintf('greatest LogLikelihood is for %d regressors',LL_best_n_reg(1,1))

%% Goodness of fit measures for best model

% Partial Slopes - Standard Errors Bx
PartialStdErr = diag(sqrt(CovB{LL_best_n_reg(1,1),LL_best_n_reg(1,2)})); % mvregress
%t-statistics (Beta / PartialStdErr)
best_beta = beta{LL_best_n_reg(1,1),LL_best_n_reg(1,2)};
tRatio = best_beta ./ PartialStdErr;
%p-values
pVals = 2*(1-tcdf(abs(tRatio),n-2));
% Summary Table of the Regression Model
RegressionSummary = [best_beta, PartialStdErr, tRatio, pVals];

%% Plotting

beta_cmbnd_mdl = [-1.7215;.9410;-.1345;.5543;-.0724;-.1204;.0056]; 
SE_cmbnd_mdl = [1.0015;.3615;.4258;.3334;.1517;.1148;.1372]; %beta estimates from combined full model
betaMat = reshape(best_beta(5:end)',4,[]); %for bar plots
betaMat = [beta_cmbnd_mdl,betaMat'];
STDerrors_neg = reshape(PartialStdErr(5:end)',4,[]); %for errorbars 
STDerrors_neg = [SE_cmbnd_mdl,STDerrors_neg'];
STDerrors_pos = STDerrors_neg; 
STDerrors_pos(betaMat<0) = 0; STDerrors_neg(betaMat>0) = 0; 
xErrorbarMat = repmat([.69,.85,1,1.15,1.3],7,1) + [0;1;2;3;4;5;6]*ones(1,5);
star_locs_x = xErrorbarMat(:,2:end)'; star_locs_x = star_locs_x(pVals(5:end,:) < .05);
star_locs_x = [1.69,star_locs_x']-.05;
star_locs_y = [2.8,2.8,2.8,2.8];
figure
hold on
bar(betaMat)
errorbar(xErrorbarMat,betaMat, STDerrors_neg,STDerrors_pos,'k.','LineWidth',.5);
text(star_locs_x', star_locs_y, '*','FontSize',20);
xticks(1:7)
xticklabels(structs);
xtickangle(45);
ylim([-5 3]);
grid on
ylabel('Beta Estimates');
legend({'Combined','Cond. 1','Cond. 2','Cond. 3','Cond. 4'});

