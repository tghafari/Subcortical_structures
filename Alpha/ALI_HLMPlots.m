clc; clear; close all

%% dot plot + average
% saveFolderMat = '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/'; %Portal
saveFolderMat = '/Volumes/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/'; %Mac

cd(saveFolderMat)
numSubj = setxor(1:35,[28,23]);
modulationIdx = zeros(length(numSubj),5); %MI_R; MI_L; MI_diff (LI); MI_add(HLM); MI group; HLM group

%rank and median split MI_R-MI_L
for subj=1:length(numSubj)        
    load([saveFolderMat 'indiv_level/Modulation_Index/Sub' num2str(numSubj(subj)) filesep 'MI_Data_dt'])
    modulationIdx(subj,1) = MI_R;
    modulationIdx(subj,2) = MI_L;
    modulationIdx(subj,3) = MI_R-MI_L;    %ALI
    modulationIdx(subj,4) = MI_R+MI_L;    %HLM
end

ALI_median = median(modulationIdx(:,3));
HLM_median = median(modulationIdx(:,4));
modulationIdx(modulationIdx(:,3)>ALI_median,5)  = 1; %higher ALI group
modulationIdx(modulationIdx(:,3)<=ALI_median,5) = 2; %lower ALI group
modulationIdx(modulationIdx(:,4)>HLM_median,6)  = 1; %higher HLM group
modulationIdx(modulationIdx(:,4)<=HLM_median,6) = 2; %lower HLM group
% modulationIdx(end+1,:)=mean(modulationIdx(:,:));

save([saveFolderMat 'group_level/Lateralization_indices' filesep 'MI_all_dt'],'modulationIdx')

%% correlation between each sub structure volume and MI_diff for all subjects
load([saveFolderMat 'group_level/Lateralization_indices/LV_all.mat'])
structures={'Thal','Caud','Puta','Pall','Hipp','Amyg','Accu'};
LV([23,28],:)=[];
LV_median = median(LV(:,:));

for ii=1:7
LVGroup(LV(:,ii)>LV_median(ii),ii)  = 1;
LVGroup(LV(:,ii)<=LV_median(ii),ii) = 2;
end

figure()
for subStr=1:7
    subplot(2,4,subStr)
    hold on
%       %LV ALI correlation    
%     scatter(modulationIdx(:,3),LV(:,subStr),'filled','r')

%     %group on MI-LV ALI correlation
%     scatter(modulationIdx(modulationIdx(:,5)==1,3),LV(modulationIdx(:,5)==1,subStr),'filled')
%     scatter(modulationIdx(modulationIdx(:,5)==2,3),LV(modulationIdx(:,5)==2,subStr),'filled')
%     %group on LV-LV ALI correlation
%     scatter(modulationIdx(LVGroup(1:30,subStr)==1,3),LV(LVGroup(1:30,subStr)==1,subStr),'filled')
%     scatter(modulationIdx(LVGroup(1:30,subStr)==2,3),LV(LVGroup(1:30,subStr)==2,subStr),'filled')
%     %group on HLM-LV ALI correlation
%     scatter(modulationIdx(modulationIdx(1:33,6)==1,3)...
%         ,LV(modulationIdx(1:33,6)==1,subStr),'b','filled')
%     scatter(modulationIdx(modulationIdx(1:33,6)==2,3)...
%         ,LV(modulationIdx(1:33,6)==2,subStr),'r','filled')
%     scatter(mean(modulationIdx(modulationIdx(1:33,6)==1,3)),mean(LV(modulationIdx(1:33,6)==1,3))...
%         ,120,'b','d','filled')
%     scatter(mean(modulationIdx(modulationIdx(1:33,6)==2,3)),mean(LV(modulationIdx(1:33,6)==2,3))...
%         ,120,'r','d','filled')

%     LV HLM correlation    
    scatter(modulationIdx(:,4),LV(:,subStr),'filled','b')
%     %group on ALI-LV HLM correlation
%     scatter(modulationIdx(modulationIdx(1:33,5)==1,4)...
%         ,LV(modulationIdx(1:33,5)==1,subStr),'b','filled')
%     scatter(modulationIdx(modulationIdx(1:33,5)==2,4)...
%         ,LV(modulationIdx(1:33,5)==2,subStr),'r','filled')
%     scatter(mean(modulationIdx(modulationIdx(1:33,5)==1,4)),mean(LV(modulationIdx(1:33,5)==1,3))...
%         ,120,'b','d','filled')
%     scatter(mean(modulationIdx(modulationIdx(1:33,5)==2,4)),mean(LV(modulationIdx(1:33,5)==2,3))...
%         ,120,'r','d','filled')

           
%     xlabel('Alpha Lateralization Index'); ylabel('Lateralization Volume')
    xlabel('Hemispheric Lateralization Modulation'); ylabel('Lateralization Volume Index')
    title(structures{subStr})
end
% legend('high ALI group','low ALI group','high ALI average','low ALI average')
% legend('high LV group','low LV group')
% legend('high HLM group','low HLM group','high HLM average','low HLM average')

%% Statistical analysis
h=zeros(7,1);   h2=zeros(7,1);
p=zeros(7,1);   p2=zeros(7,1);
rhoHLM=zeros(7,1); pHLM=zeros(7,1);
rhoALI=zeros(7,1); pALI=zeros(7,1);
pearRHLM=zeros(7,1); pearPHLM=zeros(7,1);
pearRALI=zeros(7,1); pearPALI=zeros(7,1);

structures={'Thal','Caud','Puta','Pall','Hipp','Amyg','Accu'};

for subStr=1:7
% %test LV based on HLM    
% [h(subStr),p(subStr)]   = ttest(LV(modulationIdx(1:33,6)==1,subStr),LV(modulationIdx(1:32,6)==2,subStr));
% [p2(subStr),h2(subStr)] = ranksum(LV(modulationIdx(1:33,6)==1,subStr),LV(modulationIdx(1:33,6)==2,subStr));
%test LV based on ALI
[h(subStr),p(subStr)]   = ttest(LV(modulationIdx(1:33,5)==1,subStr),LV(modulationIdx(1:32,5)==2,subStr));
[p2(subStr),h2(subStr)] = ranksum(LV(modulationIdx(1:33,5)==1,subStr),LV(modulationIdx(1:33,5)==2,subStr));

%Correlation between LV an HLM/ALI
[rhoHLM(subStr),pHLM(subStr)] = corr(LV(:,subStr),modulationIdx(:,4),'Type','Spearman');
[rhoALI(subStr),pALI(subStr)] = corr(LV(:,subStr),modulationIdx(:,3),'Type','Spearman');

[pearRHLM(subStr),pearPHLM(subStr)] = corr(LV(:,subStr),modulationIdx(:,4),'Type','Pearson');
[pearRALI(subStr),pearPALI(subStr)] = corr(LV(:,subStr),modulationIdx(:,3),'Type','Pearson');

end

