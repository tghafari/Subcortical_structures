function []=D2_ALI_HLMPlots(badSubs,signal,grp_indiv,append)
% correlation plot + average
%badSubs: the code of bad subjects ([23,28] for now)
%which signal you want to run the code on? 1-> alpha 2-> RFT
%signal = alpha or RFT
%grp_indiv: do you want to run on grp ROI -> 1 or individual ROI -> 2
%append: 1-> if append the load conditions 2- if load conditions separate
%3-> load conditions separated and subtracted

%system on which we are running the script
% saveFolderMat = '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/'; %Portal
% saveFolderMat = '/Volumes/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/'; %Mac
saveFolderMat = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\'; %Windows
cd(saveFolderMat)

numSubj = setxor(1:35,badSubs);
if append == 1
    modulationIdx = zeros(length(numSubj),6); %MI_R; MI_L; MI_diff (LI); MI_add(HLM); MI group; HLM group
elseif append == 2
    modulationIdx = zeros(length(numSubj),4,4);
end

for subj=1:length(numSubj)
    if strmatch(signal,'alpha')
        %rank and median split MI_R-MI_L (MI_R is right sensors, MI_L is left
        %sensors)
        if append == 1
            load([saveFolderMat 'indiv_level/Modulation_Index/Sub' num2str(numSubj(subj)) filesep 'MI_Data_dt_sym_salient']) %'MI_Data_dt'
            modulationIdx(subj,1) = MI_R;
            modulationIdx(subj,2) = MI_L;
            modulationIdx(subj,3) = MI_R-MI_L;    %ALI
            modulationIdx(subj,4) = MI_R+MI_L;    %HLM
        elseif append == 2
            load([saveFolderMat 'indiv_level/Modulation_Index/Sub' num2str(numSubj(subj)) filesep 'MI_Data_dt_sym_load']) %'MI_Data_dt'
            for ld = 1:4
                modulationIdx(subj,1,ld) = MI_R(ld);
                modulationIdx(subj,2,ld) = MI_L(ld);
                modulationIdx(subj,3,ld) = MI_R(ld)-MI_L(ld);    %ALI
                modulationIdx(subj,4,ld) = MI_R(ld)+MI_L(ld);    %HLM
            end
        end
    else
        if grp_indiv == 1
            load([saveFolderMat 'indiv_level/Modulation_Index/Sub' num2str(numSubj(subj)) filesep 'RFT_dt_MI_Data']) %'MI_Data_dt'
            modulationIdx(subj,1) = MI_avgTime_RFT_left_right_sens;
            modulationIdx(subj,2) = MI_avgTime_RFT_right_left_sens;
            modulationIdx(subj,3) = MI_avgTime_RFT_left_right_sens - MI_avgTime_RFT_right_left_sens;    %ALI
            modulationIdx(subj,4) = MI_avgTime_RFT_left_right_sens + MI_avgTime_RFT_right_left_sens;    %HLM
        elseif grp_indiv == 2
            if append == 1
                load([saveFolderMat 'indiv_level/Modulation_Index/Sub' num2str(numSubj(subj)) filesep 'RFT_dt_MI_Data_indiv_ROI']) %'MI_Data_dt'
                modulationIdx(subj,1) = MI_avgTime_RFT_left_right_sens;
                modulationIdx(subj,2) = MI_avgTime_RFT_right_left_sens;
                modulationIdx(subj,3) = MI_avgTime_RFT_left_right_sens - MI_avgTime_RFT_right_left_sens;    %ALI
                modulationIdx(subj,4) = MI_avgTime_RFT_left_right_sens + MI_avgTime_RFT_right_left_sens;    %HLM
            else
                load([saveFolderMat 'indiv_level/Modulation_Index/Sub' num2str(numSubj(subj)) filesep 'RFT_dt_MI_Data_indiv_ROI_load']) %'MI_Data_dt'
                for ld = 1:4
                    modulationIdx(subj,1,ld) = MI_avgTime_RFT_left_right_sens(ld);
                    modulationIdx(subj,2,ld) = MI_avgTime_RFT_right_left_sens(ld);
                    modulationIdx(subj,3,ld) = MI_avgTime_RFT_left_right_sens(ld) - MI_avgTime_RFT_right_left_sens(ld);    %ALI
                    modulationIdx(subj,4,ld) = MI_avgTime_RFT_left_right_sens(ld) + MI_avgTime_RFT_right_left_sens(ld);    %HLM
                end
            end
        end
    end
end
if append == 1
    ALI_median = median(modulationIdx(:,3));
    HLM_median = median(modulationIdx(:,4));
    modulationIdx(modulationIdx(:,3)>ALI_median,5)  = 1; %higher ALI group
    modulationIdx(modulationIdx(:,3)<=ALI_median,5) = 2; %lower ALI group
    modulationIdx(modulationIdx(:,4)>HLM_median,6)  = 1; %higher HLM group
    modulationIdx(modulationIdx(:,4)<=HLM_median,6) = 2; %lower HLM group
elseif append == 2
    ALI_median(1:4) = median(modulationIdx(:,3,:));
    HLM_median(1:4) = median(modulationIdx(:,4,:));
else
    HLM_diff_load(:,1) = modulationIdx(:,4,1) -  modulationIdx(:,4,2); %low load - ligh load -> noisy distractor
    HLM_diff_load(:,2) = modulationIdx(:,4,3) -  modulationIdx(:,4,4); %low load - ligh load -> salient distractor
    HLM_diff_sal(:,1) = modulationIdx(:,4,1) -  modulationIdx(:,4,3); %noisy - salient -> low load
    HLM_diff_sal(:,2) = modulationIdx(:,4,2) -  modulationIdx(:,4,4); %noisy - salient ->high load
end

%% save
if strmatch(signal,'alpha')
    if append == 1
        save([saveFolderMat 'group_level/Lateralization_indices' filesep 'MI_all_dt_sym_salient'],'modulationIdx')
    else
        save([saveFolderMat 'group_level/Lateralization_indices' filesep 'MI_all_dt_sym_load'],'modulationIdx')
    end
else
    if grp_indiv == 1
        save([saveFolderMat 'group_level/Lateralization_indices' filesep 'RFT_MI_all_dt'],'modulationIdx')
    elseif grp_indiv == 2
        if append == 1
            save([saveFolderMat 'group_level/Lateralization_indices' filesep 'RFT_MI_indiv_dt'],'modulationIdx')
        elseif append == 2
            save([saveFolderMat 'group_level/Lateralization_indices' filesep 'RFT_MI_indiv_dt_load'],'modulationIdx')
        end
    end
end

%% load volume data
load([saveFolderMat 'group_level' filesep 'Lateralization_indices' filesep 'LV_all.mat'])
structures={'Thal','Caud','Puta','Pall','Hipp','Amyg','Accu'};
LV(badSubs,:)=[];
LV_median = median(LV(:,:));

for ii=1:7
    LVGroup(LV(:,ii)>LV_median(ii),ii)  = 1;
    LVGroup(LV(:,ii)<=LV_median(ii),ii) = 2;
end

%% correlation between each sub structure volume and MI_diff for all subjects
if append == 2
figure(ld)
    for subStr=1:7
        subplot(2,4,subStr)
        hold on
        %LV HLM correlation
        scatter(modulationIdx(:,4,ld),LV(:,subStr),'filled','b')
        xlabel('HLM'); ylabel('Lateralization Volume Index')
        title(structures{subStr})
    end
    sgtitle(['Load Condition: ' num2str(ld)],'FontWeight','bold')
elseif append == 1
    figure()
    for subStr=1:7
        subplot(2,4,subStr)
        hold on
        %       %LV ALI correlation
        %     scatter(modulationIdx(:,3),LV(:,subStr),'filled','r')
        %     xlabel('RFT-LI'); ylabel('Lateralization Volume')
        %       %group on HLM-LV ALI correlation
        %     scatter(modulationIdx(modulationIdx(1:33,6)==1,3)...
        %         ,LV(modulationIdx(1:33,6)==1,subStr),'b','filled')
        %     scatter(modulationIdx(modulationIdx(1:33,6)==2,3)...
        %         ,LV(modulationIdx(1:33,6)==2,subStr),'r','filled')
        %     scatter(mean(modulationIdx(modulationIdx(1:33,6)==1,3)),mean(LV(modulationIdx(1:33,6)==1,3))...
        %         ,120,'b','d','filled')
        %     scatter(mean(modulationIdx(modulationIdx(1:33,6)==2,3)),mean(LV(modulationIdx(1:33,6)==2,3))...
        %         ,120,'r','d','filled')
        %     xlabel('RFT-LI'); ylabel('Lateralization Volume')
        %LV HLM correlation
        scatter(modulationIdx(:,4),LV(:,subStr),'filled','b')
        xlabel('HLM'); ylabel('Lateralization Volume Index')
        %       %group on ALI-LV HLM correlation
        %     scatter(modulationIdx(modulationIdx(1:33,5)==1,4)...
        %         ,LV(modulationIdx(1:33,5)==1,subStr),'b','filled')
        %     scatter(modulationIdx(modulationIdx(1:33,5)==2,4)...
        %         ,LV(modulationIdx(1:33,5)==2,subStr),'r','filled')
        %     scatter(mean(modulationIdx(modulationIdx(1:33,5)==1,4)),mean(LV(modulationIdx(1:33,5)==1,3))...
        %         ,120,'b','d','filled')
        %     scatter(mean(modulationIdx(modulationIdx(1:33,5)==2,4)),mean(LV(modulationIdx(1:33,5)==2,3))...
        %         ,120,'r','d','filled')
        %     xlabel('RFT-HLM'); ylabel('Lateralization Volume Index')        
        title(structures{subStr})
    end 
else
    for subStr = 1:7
        figure(1); subplot(2,4,subStr)
        hold on
        %LV HLM correlation
        scatter(HLM_diff_load(:,2),LV(:,subStr),'filled','b')
        xlabel('HLM (salient distr, LL-HL)'); ylabel('Lateralization Volume Index')
        title(structures{subStr})
        figure(2); subplot(2,4,subStr)
        hold on
        %LV HLM correlation
        scatter(HLM_diff_sal(:,2),LV(:,subStr),'filled','b')
        xlabel('HLM (High load, sal - noisy)'); ylabel('Lateralization Volume Index')
        title(structures{subStr})
    end
end
% legend('high HLM group','low HLM group','high HLM average','low HLM average')
% legend('high ALI group','low ALI group','high ALI average','low ALI average')
% legend('high LV group','low LV group')

%% All subs LIs
%each subjects HLM and ALI
figure()
if append == 2
    for ld = 1:4
        subplot(2,2,ld)
        bar(1:length(numSubj), [modulationIdx(:,3,ld),modulationIdx(:,4,ld)],1) %LI and HLM
        title(['Load Condition ' num2str(ld)])
        xlabel('subject number'); ylabel('ALI & HLM')
        grid on
        legend('ALI', 'HLM')
    end
else
    bar(1:length(numSubj), [modulationIdx(:,3),modulationIdx(:,4)],1) %LI and HLM
    xlabel('subject number'); ylabel('ALI & HLM')
    grid on
    legend('ALI', 'HLM')
end
%histogram of distribution of HLMs
clmp=colormap('jet');
figure()
if append == 2
    for ld = 1:4
        subplot(2,2,ld)
        h=histfit(modulationIdx(:,4,ld),5,'normal');
        h(1).FaceColor=clmp(30*ld,:);
        p(ld)=signrank(modulationIdx(:,4,ld));
        title(['Load Condition ' num2str(ld)])
        txt=sprintf('p= %.4f',p(ld));
        text(-.015,18,txt)
        ylim([0 20])
        xlabel('HLM'); ylabel('Nr. of Subjects')
    end
elseif append == 1
    histfit(modulationIdx(:,4),5,'normal');
    [~, pValue, W]=swtest(modulationIdx(:,4),0.05); %Shapiro-Wilk test
    txt=sprintf('p= %.4f',pValue);
    text(-.015,18,txt)
    ylim([0 20])
    xlabel('HLM'); ylabel('Nr. of Subjects')
else
    histfit(HLM_diff_load(:,2),5,'normal');
    p=signrank(HLM_diff_load(:,2)); %wilcoxon rank sum (compare with median = 0)
    txt=sprintf('p= %.4f',p);
    text(-.015,18,txt)
    ylim([0 20])
    xlabel('HLM_LL - HLM_HL in salient distractor'); ylabel('Nr. of Subjects')    
end

%% Statistical analysis
if append == 1
    % h=zeros(7,1);   h2=zeros(7,1);
    % p2_HLM=zeros(7,1);   p_HLM=zeros(7,1);
    % p2_LI=zeros(7,1);   p_LI=zeros(7,1);
    rhoHLM=zeros(7,1); pHLM=zeros(7,1);
    rhoALI=zeros(7,1); pALI=zeros(7,1);
    pearRHLM=cell(7,4); pearPHLM=cell(7,4);
    RL_HLM=cell(7,4); RU_HLM=cell(7,4);
    pearRALI=zeros(7,1); pearPALI=zeros(7,1);
elseif append == 2
    rhoHLM=zeros(7,1,4); pHLM=zeros(7,1,4);
    rhoALI=zeros(7,1,4); pALI=zeros(7,1,4);
    pearRHLM_ld=cell(7,4); pearPHLM_ld=cell(7,4);
    RL_HLM_ld=cell(7,4); RU_HLM_ld=cell(7,4);
    pearRALI_ld=zeros(7,1); pearPALI_ld=zeros(7,1);
else %difference between HLM in salience/load
    pearRHLM_diff_ld=cell(7,1); pearPHLM_diff_ld=cell(7,1);
    RL_HLM_diff_ld=cell(7,1); RU_HLM_diff_ld=cell(7,1);
    pearRHLM_diff_sal=cell(7,1); pearPHLM_diff_sal=cell(7,1);
    RL_HLM_diff_sal=cell(7,1); RU_HLM_diff_sal=cell(7,1);
end

structures={'Thal','Caud','Puta','Pall','Hipp','Amyg','Accu'};

for subStr=1:7
    if append == 1
        %     % test LV based on HLM
        %     [h(subStr),p_HLM(subStr)]   = ttest2(LV(modulationIdx(1:33,6)==1,subStr),LV(modulationIdx(1:33,6)==2,subStr));
        %     [p2_HLM(subStr),h2(subStr)] = ranksum(LV(modulationIdx(1:33,6)==1,subStr),LV(modulationIdx(1:33,6)==2,subStr));
        %
        %     %test LV based on ALI
        %     [h(subStr),p_LI(subStr)]   = ttest2(LV(modulationIdx(1:33,5)==1,subStr),LV(modulationIdx(1:33,5)==2,subStr));
        %     [p2_LI(subStr),h2(subStr)] = ranksum(LV(modulationIdx(1:33,5)==1,subStr),LV(modulationIdx(1:33,5)==2,subStr));
        
        %Correlation between LV an HLM/ALI
        [rhoHLM(subStr),pHLM(subStr)] = corr(LV(:,subStr),modulationIdx(:,4),'Type','Spearman');
        [rhoALI(subStr),pALI(subStr)] = corr(LV(:,subStr),modulationIdx(:,3),'Type','Spearman');
       
        [pearRHLM{subStr},pearPHLM{subStr},RL_HLM{subStr},RU_HLM{subStr}] = corrcoef(LV(:,subStr),modulationIdx(:,4));
        [pearRALI(subStr),pearPALI(subStr)] = corr(LV(:,subStr),modulationIdx(:,3),'Type','Pearson');
               
    elseif append == 2
        %Correlation between LV an HLM/ALI
        for ld = 1:4
            [rhoHLM(subStr,ld),pHLM(subStr,ld)] = corr(LV(:,subStr),modulationIdx(:,4,ld),'Type','Spearman');
            [rhoALI(subStr,ld),pALI(subStr,ld)] = corr(LV(:,subStr),modulationIdx(:,3,ld),'Type','Spearman');
            
            [pearRHLM_ld{subStr,ld},pearPHLM_ld{subStr,ld},RL_HLM_ld{subStr,ld},RU_HLM_ld{subStr,ld}] = corrcoef(LV(:,subStr),modulationIdx(:,4,ld));
            [pearRALI_ld(subStr,ld),pearPALI_ld(subStr,ld)] = corr(LV(:,subStr),modulationIdx(:,3,ld),'Type','Pearson');
        end
    else
        %Correlation between LV and HLM difference in salience and load
        [pearRHLM_diff_ld{subStr},pearPHLM_diff_ld{subStr},RL_HLM_diff_ld{subStr},RU_HLM_diff_ld{subStr}] = corrcoef(LV(:,subStr),HLM_diff_load(:,2));        
        [pearRHLM_diff_sal{subStr},pearPHLM_diff_sal{subStr},RL_HLM_diff_sal{subStr},RU_HLM_diff_sal{subStr}] = corrcoef(LV(:,subStr),HLM_diff_sal(:,2));        
    end
end

%% average MI over all subs on each side's ROI
MI_overTime_all_L = zeros(101,length(numSubj));
MI_overTime_all_R = zeros(101,length(numSubj));

for subj = numSubj
    load([saveFolderMat '\indiv_level\Modulation_Index\Sub' num2str(subj) '\All_Data_dt_sym_salient'])
    MI_overTime_all_L(:,subj) = MI_overTime_inROI_L;
    MI_overTime_all_R(:,subj) = MI_overTime_inROI_R;
end

save([saveFolderMat '\group_level\Lateralization_indices\MI_eachSub_eachSide_dt_sym_salient'],'MI_overTime_all_R','MI_overTime_all_L')
% plot average MI + std
MI_overTime_all_L(:,badSubs)=[];
MI_overTime_all_R(:,badSubs)=[];

figure()
subplot(1,2,1)
errorbar(mean(MI_overTime_all_L,2),std(MI_overTime_all_L,[],2)/sqrt(length(numSubj)),'Color',[0 0 .7])
ylim([-.15 .15]); xlim([0 86]); xticks(0:25:86); xticklabels([-0.850,-0.600,-0.350,-0.100])
subplot(1,2,2)
errorbar(mean(MI_overTime_all_R,2),std(MI_overTime_all_R,[],2)/sqrt(length(numSubj)),'Color',[.2 .4 0])
ylim([-.15 .15]); xlim([0 86]); xticks(0:25:86); xticklabels([-0.850,-0.600,-0.350,-0.100])

