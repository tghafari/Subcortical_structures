clc; clear; close all

%% dot plot + average
saveFolder = 'Z:\Load\Results\FieldTrip Plots\matFiles';
cd(saveFolder)
numSubj = setxor(1:30,28);
modulationIdx = zeros(length(numSubj),5); %MI_R; MI_L; MI_diff (LI); MI_add(HLM); MI group

%rank and median split MI_R-MI_L
for subj=1:length(numSubj)        
    load([saveFolder filesep 'Sub' num2str(numSubj(subj)) filesep 'MI_Data'])
    modulationIdx(subj,1) = MI_R;
    modulationIdx(subj,2) = MI_L;
    modulationIdx(subj,3) = MI_R-MI_L;    %? it is Right - Left in the paper, right?
    modulationIdx(subj,4) = MI_R+MI_L;
end

MI_median = median(modulationIdx(:,3));
modulationIdx(modulationIdx(:,3)>MI_median,5)  = 1; %higher MI group
modulationIdx(modulationIdx(:,3)<=MI_median,5) = 2; %lower MI group
modulationIdx(end+1,:)=mean(modulationIdx(:,:));

save([saveFolder filesep 'AllSubFiles' filesep 'MI_all'],'modulationIdx')

%% correlation between each sub structure volume and MI_diff for all subjects

structures={'Thal','Caud','Puta','Pall','Hipp','Amyg','Accu'};
LV([23,28],:)=[];
LV_median = median(LV(:,:));

for ii=1:7
LVGroup(LV(:,ii)>LV_median(ii),ii)  = 1;
LVGroup(LV(:,ii)<=LV_median(ii),ii) = 2;
end

figure(1)
for subStr=1:7
    subplot(2,4,subStr)
    hold on
%     %group on MI
%     scatter(modulationIdx(modulationIdx(:,5)==1,3),LV(modulationIdx(:,5)==1,subStr),'filled')
%     scatter(modulationIdx(modulationIdx(:,5)==2,3),LV(modulationIdx(:,5)==2,subStr),'filled')
    %group on LV
    scatter(modulationIdx(LVGroup(1:30,subStr)==1,3),LV(LVGroup(1:30,subStr)==1,subStr),'filled')
    scatter(modulationIdx(LVGroup(1:30,subStr)==2,3),LV(LVGroup(1:30,subStr)==2,subStr),'filled')

    xlabel('Alpha Lateralization Index'); ylabel('Lateralization Volume')
    title(structures{subStr})
end
% legend('high ALI group','low ALI group')
legend('high LV group','low LV group')