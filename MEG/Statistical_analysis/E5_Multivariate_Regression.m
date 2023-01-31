% function [RegressionSummary] = multivariateregression(Y,X)
%Multivariate multiple regression on DVar Y along dimensions d and IVar X along dimensions p. 
%Y is the dependent variable, each row subject, each column is condition. Xmat is matrix containing IV, should be all n x p, where n is subject
%and every column p is different predictor
clear
%Windows
saveFolderMat = 'Z:\Programming\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices\'; 
load([saveFolderMat 'LV_all.mat'])
load([saveFolderMat '/MI_all_dt_sym_load.mat'])

%Mac
load('/Users/Tara/Documents/MATLAB/MATLAB-Programs/CHBH-Programs-move to RDS/Perceptual-Load/FieldTrip-Just checking/MI_all_dt_sym_load.mat')
load('/Users/Tara/Documents/MATLAB/MATLAB-Programs/CHBH-Programs-move to RDS/Perceptual-Load/FieldTrip-Just checking/LV_all.mat')

strcts = {'Th', 'CN', 'Pu', 'GP', 'Hpc', 'Amg', 'Acb'};
slctd_strc = [1,2,4];

%% checking the correlation between structures 

%pairwise correlation
[R,P]=corrcoef(LV(setxor(1:35,[23,28]),:));
cntr = 0;
for stri = 1:numel(strcts)-1
    for strj = stri+1:numel(strcts)
        cntr = cntr+1;
        subplot(4,6,cntr)
        scatter(LV(:,stri),LV(:,strj),'filled','SizeData',50)
        xlabel(strcts{stri}); ylabel(strcts{strj});
%         p_txt = sprintf('p-value = %.3f',P(stri,strj));
%         r_txt = sprintf('Rho = %.3f',R(stri,strj));
%         text(.02,.02,r_txt)
%         text(.02,0,p_txt)
    end
end

%varianve inflation factor (VIF)
multicol_vif = vif(LV(setxor(1:35,[23,28]),:));
multicol_vif2 = diag(inv(corrcoef(LV(setxor(1:35,[23,28]),:))))'; %same result with previous line

%% multivariate regression
Y = squeeze(modulationIdx(:,4,:)); % to get the HLM values of 4 conditions separetely
X = LV(setxor(1:35,[23,28]),slctd_strc); 
X_null = LV(setxor(1:35,[23,28]),[3,5,6,7]);

[n,d] = size(Y);
Xmat  = [ones(n,1) X]; %add intercept
Xmat_null = [ones(n,1) X_null];
% p     = size(Xmat,2)-1;

Xcell = cell(1,n);
for i = 1:n
    Xcell{i} = [kron([Xmat(i,:)],eye(d))];
end
[beta,Sigma,E,CovB,logL] = mvregress(Xcell,Y);

% Hypothesis = beta(1) + beta(2) * X(:,1) +  beta(3) * X(:,2) + ...
% beta(4) *  X(:,3);
% 
% % r-squared
% r2 = sum((Hypothesis - mean(Y)).^2) / sum((Y - mean(Y)).^2);
% number of explanatory variables
% Data from the Regression Model
% beta(1); % intercept
% beta(2); % X1
% beta(3); % X2
% beta(4); % X3

% % model Standard Error
% Se2 = sum((Y - beta(1) - beta(2) * X(:,1) -  beta(3) * X(:,2) + ...
%     beta(4) *  X(:,3)).^2) / (n-p-1);
% Se = sqrt(Se2);

%% statistical analysis of selected model
% Partial Slopes - Standard Errors Bx
PartialStdErr = diag(sqrt(CovB)); % mvregress
%t-statistics (Beta / PartialStdErr)
tRatio = beta ./ PartialStdErr;
%p-values
pVals = 2*(1-tcdf(abs(tRatio),n-2));
% Summary Table of the Regression Model
RegressionSummary = [beta, PartialStdErr, tRatio, pVals];

%% plotting

% plot residuals
z = E/chol(Sigma);
figure()
plot(z(:,1),z(:,2),'.')
title('Standardized Residuals')
hold on

% Overlay standard normal contours
z1 = linspace(-5,5);
z2 = linspace(-5,5);
[zx,zy] = meshgrid(z1,z2);
zgrid = [reshape(zx,100^2,1),reshape(zy,100^2,1)];
zn = reshape(mvnpdf(zgrid),100,100);
[c,h] = contour(zx,zy,zn);
clabel(c,h)

% barplot beta estimates (Ole's plot)
beta_cmbnd_mdl = [-2.1912;.9255;.5062]; SE_cmbnd_mdl = [.8000;.3274;.2593]; %beta estimates from combined model
betaMat = reshape(beta(5:end)',4,[]); %for bar plots
betaMat = [beta_cmbnd_mdl,betaMat'];
STDerrors_neg = reshape(PartialStdErr(5:end)',4,[]); %for errorbars 
STDerrors_neg = [SE_cmbnd_mdl,STDerrors_neg'];
STDerrors_pos = STDerrors_neg; 
STDerrors_pos(1,:) = 0; STDerrors_pos(3,end) = 0;
STDerrors_neg(2,:) = 0; STDerrors_neg(3,1:end-1) = 0;
% xErrorbarMat = repmat([.72,.9,1.1,1.28],length(slctd_strc),1) + [0,0,0,0;1,1,1,1;2,2,2,2;3,3,3,3;4,4,4,4;5,5,5,5;6,6,6,6]; %full model
xErrorbarMat = repmat([.69,.85,1,1.15,1.3],length(slctd_strc),1) + [0,0,0,0,0;1,1,1,1,1;2,2,2,2,2]; 
star_locs_x = xErrorbarMat(:,2:end)'; star_locs_x = star_locs_x(pVals(5:end,:) < .05);
star_locs_x = [.69,1.69,star_locs_x']-.02;
star_locs_y = [-5.5,2.5,-5.5,2.5,2.5,2.5];
figure
hold on
bar(betaMat_1')
errorbar(xErrorbarMat,betaMat,STDerrors_neg, STDerrors_pos,'k.','LineWidth',1);
text(star_locs_x', star_locs_y, '*','FontSize',20);
xticks([1,2,3])
xticklabels(strcts(slctd_strc)); 
xtickangle(45);
ylim([-6.5 3.5]);
grid on
ylabel('Beta Estimates');
legend({'Combined','Cond. 1','Cond. 2','Cond. 3','Cond. 4'});

% barplot- remove combined beta and put GP in the middle
betaMat = reshape(beta(5:end)',4,[]); %for bar plots
betaMat = betaMat(:,[1,3,2]);
STDerrors_neg = reshape(PartialStdErr(5:end)',4,[]); %for errorbars 
STDerrors_neg = STDerrors_neg(:,[1,3,2]);
STDerrors_pos = STDerrors_neg; 
STDerrors_pos(:,1) = 0; STDerrors_pos(end,2) = 0;
STDerrors_neg(:,3) = 0; STDerrors_neg(1:end-1,2) = 0;
xErrorbarMat = repmat([.72,.91,1.1,1.28],length(slctd_strc),1) + [0,0,0,0;1,1,1,1;2,2,2,2]; 
pVals = reshape(pVals(5:end)',4,[]); pVals = pVals(:,[1,3,2]); 
star_locs_x = xErrorbarMat'; star_locs_x = star_locs_x(pVals(:) < .05);
star_locs_x = star_locs_x'-.02;
star_locs_y = [-5.5,2.5,2.5,2.5];
figure
hold on
bar(betaMat')
errorbar(xErrorbarMat,betaMat',STDerrors_neg', STDerrors_pos','k.','LineWidth',1);
text(star_locs_x, star_locs_y, '*','FontSize',20);
xticks([1,2,3])
xticklabels(strcts(slctd_strc([1, 3, 2]))); 
xtickangle(45);
ylim([-6 3]);
grid on
ylabel('Beta Estimates');
legend({'Cond. 1','Cond. 2','Cond. 3','Cond. 4'});



%% saving mvregress outputs
pV = array2table(reshape(pVals,4,[]),'VariableNames',{'pVal:Intercept','Th','CN','GP'},'RowNames',{'condition 1','condition 2','condition 3','coniditon 4'});
tR = array2table(reshape(tRatio,4,[]),'VariableNames',{'tRatio:Intercept','Th','CN','GP'},'RowNames',{'condition 1','condition 2','condition 3','coniditon 4'});
stdev = array2table(reshape(PartialStdErr,4,[]),'VariableNames',{'std:Intercept','Th','CN','GP'},'RowNames',{'condition 1','condition 2','condition 3','coniditon 4'});
bt = array2table(reshape(beta,4,[]),'VariableNames',{'beta:Intercept','Th','CN','GP'},'RowNames',{'condition 1','condition 2','condition 3','coniditon 4'});
sig = array2table(reshape(Sigma,4,[]),'VariableNames',{'sigma:Intercept','Th','CN','GP'},'RowNames',{'condition 1','condition 2','condition 3','coniditon 4'});
res = array2table(E,'VariableNames',{'E:Intercept','Th','CN','GP'});
cov=array2table(CovB,'VariableNames',{'Cov:Intercept 1','Intercept 2','Intercept 3','Intercept 4','Th 1','Th 2','Th 3','Th 4','CN 1','CN 2','CN 3','CN 4','GP 1','GP 2','GP 3','GP 4'},'RowNames',{'Intercept 1','Intercept 2','Intercept 3','Intercept 4','Th 1','Th 2','Th 3','Th 4','CN 1','CN 2','CN 3','CN 4','GP 1','GP 2','GP 3','GP 4'});
outputStr = {bt,sig,res,cov,pV,tR,stdev};

