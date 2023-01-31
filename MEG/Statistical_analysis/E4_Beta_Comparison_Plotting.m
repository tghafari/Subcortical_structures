function [] = E4_Beta_Comparison_Plotting(slctd_strct,n_BS,plotting)

% plots bootstrapped beta estimates and their differences in conditions and structures
% slctd_strct = which structures to do comarison on ([1,2,4]; %Th, CN and GP: selected based
%                                                    on best model on all conditions combined)
% n_BS = number of bootstrapping iterations (n_BS = 10000;)
% plotting = plot beta distributions in cnd and str (1) OR beta differences
%            in cond and str (2)

saveFolderMat = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices\'; %Windows
load([saveFolderMat filesep 'Beta_pVal_BS']); %n_BS by structure by load condition
structs = {'Th', 'CN', 'Pu', 'GP', 'Hpc', 'Amg', 'Acb'};
%Preallocation
diff_conds = nan(n_BS+1,3,6);  %difference between conditions -> 3-4,2-4,2-3,1-4,1-3,1-2
diff_strcts = nan(n_BS+1,3,4); %difference between structures in each cond -> Th-CN, Th-GP, CN-GP


%% Statistical analysis according to: https://garstats.wordpress.com/2017/03/01/comp2dcorr/
alpha = .05; % CI = quantile(diff_betas,[alpha/2, 1-alpha/2]);
lower = round(alpha/2*n_BS); upper = n_BS-lower;

%sort beta estimates for structs and conditions
sortd_beta = sort(dist_beta);
CI_betas(:,:,:) = [sortd_beta(lower+1,:,:); sortd_beta(upper,:,:)]; %[upper;lower] for each strcut (col) in each cond (dim=3)

%conditions difference
diff_ld_cmbs = combnk(1:4,2); %subtract conditions two-by-two for all structs
for ld_cmbs = 1:6
diff_conds(:,:,ld_cmbs) = dist_beta(:,:,diff_ld_cmbs(ld_cmbs,1)) - dist_beta(:,:,diff_ld_cmbs(ld_cmbs,2));
end
sortd_diff_conds = sort(diff_conds);
CI_cnds(:,:,:) = [sortd_diff_conds(lower+1,:,:); sortd_diff_conds(upper,:,:)]; %[upper;lower] for each strcut (col) in each cond (dim=3)

%structures difference
diff_strct_cmbs = combnk(1:3,2); %subtract structrs 2by2 in each condition
for str_cmbs = 1:3
diff_strcts(:,str_cmbs,:) = dist_beta(:,diff_strct_cmbs(str_cmbs,1),:) - dist_beta(:,diff_strct_cmbs(str_cmbs,2),:);
end
sortd_diff_strcts = sort(diff_strcts);
CI_strcts(:,:,:) = [sortd_diff_strcts(lower+1,:,:); sortd_diff_strcts(upper,:,:)]; %[upper;lower] for each strcut (col) in each cond (dim=3)

%% illustrate bootstrapped betas
if plotting == 1
%Preallocation
pVal_B = nan(3,4); p_value_B_cnd_str = nan(3,4); %works for both beta in conds and in structs

%for conditions
cmp_cnd = colormap('jet');
for str = 1:3
    figure; hold on
    for B_cnd = 1:4
        %distribution of differences
        subplot(2,1,1); hold on
        h(B_cnd) = histogram(dist_beta(:,str,B_cnd),50);
        h(B_cnd).FaceColor = cmp_cnd(B_cnd*25,:); h(B_cnd).EdgeColor = 'w';
        xlim([floor(min(min(dist_beta(:,str,:))))  ceil(max(max(dist_beta(:,str,:))))])
        %confidence intervals
        subplot(2,1,2); hold on
        plot([dist_beta(end,str,B_cnd),dist_beta(end,str,B_cnd)],[B_cnd-.5,B_cnd+.5],'Color',cmp_cnd(B_cnd*25,:),'LineWidth',2)
        plot([CI_betas(1,str,B_cnd),CI_betas(2,str,B_cnd)],[B_cnd,B_cnd],'Color',cmp_cnd(B_cnd*25,:),'LineWidth',1)
        v = axis; plot([0 0],[v(3) v(4)],'k--','LineWidth',1)
        xlim([floor(min(min(dist_beta(:,str,:))))  ceil(max(max(dist_beta(:,str,:))))])
        ylim([0 5])
        %calculate p-values
        pVal_B(str,B_cnd) = mean(sortd_beta(:,str,B_cnd)<0);
        p_value_B_cnd_str(str,B_cnd) = 2*min(pVal_B(str,B_cnd),1-pVal_B(str,B_cnd));
    end
    sgtitle(['Structure: ', structs{slctd_strct(str)}]);
    legend(h,{'condition 1','condition 2','condition 3','condition 4'})
end
beta_pval_tbl = array2table(p_value_B_cnd_str,'VariableNames',{'Cond 1','Cond 2','Cond 3','Cond 4'},...
    'RowNames',{'Th','CN','GP'});

%for structures
cmp_str = colormap('colorcube');
for ld = 1:4
    figure; hold on
    for B_Str = 1:3
        %distribution of differences
        subplot(2,1,1); hold on
        h(B_Str) = histogram(dist_beta(:,B_Str,ld),50);
        h(B_Str).FaceColor = cmp_str(B_Str*8,:); h(B_Str).EdgeColor = 'w';
        xlim([floor(min(min(dist_beta(:,:,ld))))  ceil(max(max(dist_beta(:,:,ld))))])
        %confidence intervals
        subplot(2,1,2); hold on
        plot([dist_beta(end,B_Str,ld),dist_beta(end,B_Str,ld)],[B_Str-.5,B_Str+.5],'Color',cmp_str(B_Str*8,:),'LineWidth',2)
        plot([CI_betas(1,B_Str,ld),CI_betas(2,B_Str,ld)],[B_Str,B_Str],'Color',cmp_str(B_Str*8,:),'LineWidth',1)
        v = axis; plot([0 0],[v(3) v(4)],'k--','LineWidth',1)
        xlim([floor(min(min(dist_beta(:,:,ld))))  ceil(max(max(dist_beta(:,:,ld))))])
        ylim([0 4])
        %calculate p-values (it is the same as in conds)
    end
    sgtitle(['Conditions: ', num2str(ld)]);
    legend(h,{'Th','CN','GP'})
end

%% illustrate bootstrap difference distributions
elseif plotting == 2
%Preallocation
pVal_cnd = nan(3,6); p_value_cnd = nan(3,6);
pVal_str = nan(3,4); p_value_str = nan(3,4);

%for difference in conditions
cmp_cnd = colormap('jet');
for str = 1:2
    figure; hold on
    for df_cnd = 1:6
        %distribution of differences
        subplot(2,1,1); hold on
        h(df_cnd) = histogram(diff_conds(:,str,df_cnd),50);
        h(df_cnd).FaceColor = cmp_cnd(df_cnd*25,:); h(df_cnd).EdgeColor = 'w';
        xlim([floor(min(min(diff_conds(:,str,:))))  ceil(max(max(diff_conds(:,str,:))))])
        %confidence intervals
        subplot(2,1,2); hold on
        plot([diff_conds(end,str,df_cnd),diff_conds(end,str,df_cnd)],[df_cnd-.5,df_cnd+.5],'Color',cmp_cnd(df_cnd*25,:),'LineWidth',2)
        plot([CI_cnds(1,str,df_cnd),CI_cnds(2,str,df_cnd)],[df_cnd,df_cnd],'Color',cmp_cnd(df_cnd*25,:),'LineWidth',1)
        v = axis; plot([0 0],[v(3) v(4)],'k--','LineWidth',1)
        xlim([floor(min(min(diff_conds(:,str,:))))  ceil(max(max(diff_conds(:,str,:))))])
        ylim([0 7])
        %calculate p-values
        pVal_cnd(str,df_cnd) = mean(sortd_diff_conds(:,str,df_cnd)<0);
        p_value_cnd(str,df_cnd) = 2*min(pVal_cnd(str,df_cnd),1-pVal_cnd(str,df_cnd));
    end
    sgtitle(['Structure: ', structs{slctd_strct(str)}]);
    legend(h,{'3 - 4','2 - 4','2 - 3','1 - 4','1 - 3','1 - 2'})
end
cnds_pval_tbl = array2table(p_value_cnd,'VariableNames',{'3 - 4','2 - 4','2 - 3','1 - 4','1 - 3','1 - 2'},...
    'RowNames',{'Th','CN','GP'});

%for difference in structures
cmp_str = colormap('colorcube');
for ld = 1:2
    figure; hold on
    for df_S = 1:3
        %distribution of differences
        subplot(2,1,1); hold on
        h(df_S) = histogram(diff_strcts(:,df_S,ld),50);
        h(df_S).FaceColor = cmp_str(df_S*8,:); h(df_S).EdgeColor = 'w';
        xlim([floor(min(min(diff_strcts(:,:,ld))))  ceil(max(max(diff_strcts(:,:,ld))))])
        %confidence intervals
        subplot(2,1,2); hold on
        plot([diff_strcts(end,df_S,ld),diff_strcts(end,df_S,ld)],[df_S-.5,df_S+.5],'Color',cmp_str(df_S*8,:),'LineWidth',2)
        plot([CI_strcts(1,ld),CI_strcts(2,ld)],[df_S,df_S],'Color',cmp_str(df_S*8,:),'LineWidth',1)
        v = axis; plot([0 0],[v(3) v(4)],'k--','LineWidth',1)
        xlim([floor(min(min(diff_strcts(:,:,ld))))  ceil(max(max(diff_strcts(:,:,ld))))])
        ylim([0 4])
        %calculate p-values
        pVal_str(df_S,ld) = mean(sortd_diff_strcts(:,ld)<0);
        p_value_str(df_S,ld) = 2*min(pVal_str(df_S,ld),1-pVal_str(df_S,ld));
    end
    sgtitle(['Conditions: ', num2str(ld+2)]);
%     legend(h,{'Th - CN','Th - GP','CN - GP'})
end
strs_pval_tbl = array2table(p_value_str,'VariableNames',{'Cond 3','Cond 4'},...
    'RowNames',{'Th - CN'});
end
% save([saveFolderMat filesep 'Beta_p_values'],'beta_pval_tbl','cnds_pval_tbl','strs_pval_tbl');
