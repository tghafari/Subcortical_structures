% Statistical analysis of findings

%% z-test to compare correlation coefficients

rho_1 = .262; %bring whichever rho you are looking into
rho_2 = .337;
n_1 = 6;
n_2 = 6;
mu_z1 = sqrt(log((1+rho_1)/(1-rho_1)));
mu_z2 = sqrt(log((1+rho_2)/(1-rho_2)));
sigma = sqrt((1/(n_1-3))+(1/(n_2-3)));
Z = (mu_z1 - mu_z2)/sigma;

if Z > 1.96 || Z < -1.96
    sprintf('Z = %.3f',Z)
    disp('difference is significant')
else
    disp('difference is not significant')
end

%% plot confidence interval for pearson correlations
structures={'Thal','Caud','Puta','Pall','Hipp','Amyg','Accu'};
load('Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices\Beta_estimates_alpha_load')
%Preallocation
rho_values = zeros(3,4); rho_p_values = zeros(3,4);
rho_xpos = zeros(3,4); rho_xneg = zeros(3,4);
beta_values = zeros(3,4); beta_p_values = zeros(3,4);
beta_xpos = zeros(3,4); beta_xneg = zeros(3,4);

y_disp = 1:2:7;
ypos = zeros(1,4); yneg = zeros(1,4);

selected_strcts = [1,2,4];
for subStr = 1:length(selected_strcts)
    for ld = 1:4
%         rho_values(subStr,ld)    = pearRHLM_ld{selected_strcts(subStr),ld}(2,1);
%         rho_p_values(subStr,ld)  = pearPHLM_ld{selected_strcts(subStr),ld}(2,1);
%         rho_xpos(subStr,ld)      = RU_HLM_ld{selected_strcts(subStr),ld}(2,1);
%         rho_xneg(subStr,ld)      = RL_HLM_ld{selected_strcts(subStr),ld}(2,1);
        beta_values(subStr,ld)   = table2array(beta_est{ld}(subStr+1,1));
        beta_p_values(subStr,ld) = table2array(beta_est{ld}(subStr+1,4));
        beta_xpos(subStr,ld)     = table2array(beta_est{ld}(subStr+1,2));
        beta_xneg(subStr,ld)     = beta_xpos(subStr,ld);
    end
    
%     figure(1) %rhos
%     subplot(2,2,subStr)
%     errorbar(rho_values(subStr,:),y_disp,yneg,ypos,rho_xneg(subStr,:),rho_xpos(subStr,:),'o')
%     title(structures(selected_strcts(subStr)))
%     sgtitle('Pearson Correlation in comparison')
%     yticks([1 3 5 7]);
%     yticklabels({'rho_1','rho_2','rho_3','rho_4'})
    figure(2) %betas
    subplot(2,2,subStr)
    errorbar(beta_values(subStr,:),y_disp,yneg,ypos,beta_xneg(subStr,:),beta_xpos(subStr,:),'o')
    title(structures(selected_strcts(subStr)))
    sgtitle('Beta Estimates in comparison')
    yticks([1 3 5 7]);
    yticklabels({'Beta_1','Beta_2','Beta_3','Beta_4'})
    figure(3)
    hold on
    errorbar(1:4,beta_values(subStr,:),beta_xpos(subStr,:),'LineWidth',3)
    title('Beta Estimates in Comparison')
    xlim([0 5]); xticks(1:4); xticklabels({'Condition_1','Condition_2','Condition_3','Condition_4'});xtickangle(45);
    ylim([-6 2]); ylabel('Beta estimate')
    legend({'Th','CN','GP'})
    grid on
end

%% Plotting all betas in one plot
%Preallocation
beta_values = zeros(3,4); 
beta_xpos = zeros(3,4); beta_xneg = zeros(3,4);
y_disp = 1:2:5;
ypos = zeros(1,3); yneg = zeros(1,3);

selected_strcts = [1,2,4];
for subStr = 1:length(selected_strcts)
    for ld = 1:4
        beta_values(subStr,ld) = table2array(beta_est{ld}(subStr+1,1));
        beta_xpos(subStr,ld)     = table2array(beta_est{ld}(subStr+1,2));
        beta_xneg(subStr,ld)     = beta_xpos(subStr,ld);
    end
end
figure()
for ld= 1:4
    subplot(2,2,ld)
    errorbar(beta_values(:,ld),y_disp,yneg,ypos,beta_xneg(:,ld),beta_xpos(:,ld),'o')
    title(['condition: ',num2str(ld)])
    ylim([0 6]); yticks(1:2:5); yticklabels({'Th','CN','GP'});
    xlim([-5.5 3]); xlabel('Beta estimate')
end
sgtitle('Beta Estimates in Comparison')

%% comparing rhos/betas two by two
selected_strcts = [1,2,4];
n_1 = 33;
n_2 = 33;
for subStr = 1:length(selected_strcts)
    for ld_1 = 1:3
        rho_1 = abs(rho_values(subStr,ld_1)); 
        mu_z1 = sqrt(log((1+rho_1)/(1-rho_1)));
        beta_1 = beta_values(subStr,ld_1);
        for ld_2 = ld_1+1:4
        rho_2 = abs(rho_values(subStr,ld_2));
        mu_z2 = sqrt(log((1+rho_2)/(1-rho_2)));
        sigma = sqrt((1/(n_1-3))+(1/(n_2-3)));
        Z{subStr}(ld_1,ld_2) = (mu_z1 - mu_z2)/sigma;
        beta_2 = beta_values(subStr,ld_2);
        p_beta{subStr}(ld_1,ld_2) = ranksum(beta_1,beta_2);
        end
    end
    disp(find(Z{subStr}(:) < -1.96 | Z{subStr}(:) > 1.96))
end
