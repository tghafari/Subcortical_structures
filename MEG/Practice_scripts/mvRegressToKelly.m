strcts = {'Th', 'CN', 'Pu', 'GP', 'Hpc', 'Amg', 'Acb'};
slctd_strc = [1,2,4];
%% multivariate regression
Y = squeeze(modulationIdx(:,4,:)); % to get the HLM values of 4 conditions separetely
X = LV(setxor(1:35,[23,28]),slctd_strc); %subjects 23 and 28 are excluded- the selected structures from the model are used in mvregress

[n,d] = size(Y);
Xmat  = [ones(n,1) X]; %add intercept 

Xcell = cell(1,n);
for i = 1:n
    Xcell{i} = [kron([Xmat(i,:)],eye(d))]; %to prepare the Xcell for mvregress
end
[beta,Sigma,E,CovB,logL] = mvregress(Xcell,Y);
% Partial Slopes - Standard Errors Bx
PartialStdErr = diag(sqrt(CovB)); % mvregress
%t-statistics (Beta / PartialStdErr)
tRatio = beta ./ PartialStdErr;
%p-values
pVals = 2*(1-tcdf(abs(tRatio),n-2));
% Summary Table of the Regression Model
RegressionSummary = [beta, PartialStdErr, tRatio, pVals];
