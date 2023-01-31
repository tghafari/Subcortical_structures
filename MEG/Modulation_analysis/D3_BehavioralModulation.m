function [] = D4_BehavioralModulation(badSubs,append)
% Calculates behavioral modulation indices
%badSubs = code of bad subjects ([23,28])
%append = 1-> load conditions combined 2-> load conditions separately 

data_folder = 'Z:\Load\MEG Data\proc_data\'; %Windows
save_folder = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices';
cd(save_folder)
%% Calculate behavioral index
badSubs = [23,28];
numSubs = setxor(1:35,badSubs);
if append == 1
    BAP_RT_acc = nan(35,3);
else
    BAP_RT_acc = nan(35,4,3);
end

for subj = numSubs
    if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
    disp(['loading ' sub])
    load([data_folder sub filesep sub '_all_clean_dt.mat']);fprintf('Done\n')
    
    if append == 1
        behaviorMat = data{1,1}.exp_cfg.behavior;
        % separate behavioral measure from eachother
        beh_RT_right = nanmean(behaviorMat(behaviorMat(:,2)==2,7));
        beh_RT_left = nanmean(behaviorMat(behaviorMat(:,2)==1,7));
        beh_acc_right = nanmean(behaviorMat(behaviorMat(:,2)==2,10));
        beh_acc_left = nanmean(behaviorMat(behaviorMat(:,2)==1,10));
        beh_IES_right = beh_RT_right ./ beh_acc_right;
        beh_IES_left = beh_RT_left ./ beh_acc_left;
        
        % collect all subs' Behavioral Assymetries in Performance
        BAP_RT_acc(subj,1) = (beh_RT_right - beh_RT_left) / (beh_RT_right + beh_RT_left);
        BAP_RT_acc(subj,2) = (beh_acc_right - beh_acc_left) / (beh_acc_right + beh_acc_left);
        BAP_RT_acc(subj,3) = (beh_IES_right - beh_IES_left) / (beh_IES_right + beh_IES_left);
      
        BAP_RT_acc(badSubs,:)=[];
    elseif append == 2 %load conditions separately
        
        for ld = 1:4
            % separate behavioral measure from eachother
            beh_RT_right(ld) = nanmean(data{1,1}.right{ld}.behavior(:,7));
            beh_RT_left(ld) = nanmean(data{1,1}.left{ld}.behavior(:,7));
            beh_acc_right(ld) = nanmean(data{1,1}.right{ld}.behavior(:,10));
            beh_acc_left(ld) = nanmean(data{1,1}.left{ld}.behavior(:,10));
            beh_IES_right(ld) = beh_RT_right(ld) ./ beh_acc_right(ld);
            beh_IES_left(ld) = beh_RT_left(ld) ./ beh_acc_left(ld);
            
            % collect all subs' Behavioral Assymetries in Performance
            BAP_RT_acc(subj,ld,1) = (beh_RT_right(ld) - beh_RT_left(ld)) / (beh_RT_right(ld) + beh_RT_left(ld));
            BAP_RT_acc(subj,ld,2) = (beh_acc_right(ld) - beh_acc_left(ld)) / (beh_acc_right(ld) + beh_acc_left(ld));
            BAP_RT_acc(subj,ld,3) = (beh_IES_right(ld) - beh_IES_left(ld)) / (beh_IES_right(ld) + beh_IES_left(ld));
        end
        BAP_RT_acc(badSubs,:,:)=[];
    end
end
% if append == 1
% save([save_folder filesep 'BAP_RT_acc_allsubs'],'BAP_RT_acc')    
% else
% save([save_folder filesep 'BAP_RT_acc_allsubs_load'],'BAP_RT_acc')
% end

%% Load volume data
load([save_folder filesep 'LV_all.mat'])
structures={'Thal','Caud','Puta','Pall','Hipp','Amyg','Accu'};
LV(badSubs,:)=[];
LV_median = median(LV(:,:));

for ii=1:7
    LVGroup(LV(:,ii)>LV_median(ii),ii)  = 1;
    LVGroup(LV(:,ii)<=LV_median(ii),ii) = 2;
end

%% Correlation between each sub structure volume and behavioral measures for all subjects
if append == 1
    load([save_folder filesep 'BAP_RT_acc_allsubs.mat'])
else
    load([save_folder filesep 'BAP_RT_acc_allsubs_load.mat'])
end

for subStr=1:7
    if append == 1
        figure(1)
        subplot(2,4,subStr)
        hold on
        %LV RT correlation
        scatter(BAP_RT_acc(:,1),LV(:,subStr),'filled','r')
        xlabel('RT'); ylabel('Lateralization Volume')
        title(structures{subStr})
        
        figure(2)
        subplot(2,4,subStr)
        hold on
        %LV acc correlation
        scatter(BAP_RT_acc(:,2),LV(:,subStr),'filled','r')
        xlabel('acc'); ylabel('Lateralization Volume')
        title(structures{subStr})
        
    elseif append == 2
        for ld = 1:4
            figure(ld)
            subplot(2,4,subStr)
            hold on
            %LV RT correlation
            scatter(BAP_RT_acc(:,ld,1),LV(:,subStr),'filled','b')
            xlabel('RT'); ylabel('Lateralization Volume')
            title(structures{subStr}); sgtitle(['condition: ' num2str(ld)])
            
            figure(ld+4)
            subplot(2,4,subStr)
            hold on
            %LV acc correlation
            scatter(BAP_RT_acc(:,ld,2),LV(:,subStr),'filled','b')
            xlabel('acc'); ylabel('Lateralization Volume')
            title(structures{subStr}); sgtitle(['condition: ' num2str(ld)])
            
            figure(ld+8)
            subplot(2,4,subStr)
            hold on
            %LV IES correlation
            scatter(BAP_RT_acc(:,ld,3),LV(:,subStr),'filled','b')
            xlabel('IES'); ylabel('Lateralization Volume')
            title(structures{subStr}); sgtitle(['condition: ' num2str(ld)])
        end
    end
end

%% All subs BAPs
figure()
bar(1:length(numSubs), [BAP_RT_acc(:,1),BAP_RT_acc(:,2)],1) %RT and acc
xlabel('subject number'); ylabel('RT & accuracy')
grid on
legend('RT', 'accuracy')

%% Findings outliers
if append == 1
    upLim = mean(BAP_RT_acc(:,2)) + 2.5*std(BAP_RT_acc(:,2));
    downLim = mean(BAP_RT_acc(:,2)) - 2.5*std(BAP_RT_acc(:,2));
    for sub = 1:length(numSubs)
    if BAP_RT_acc(sub,2) >= upLim || BAP_RT_acc(sub,2) <= downLim
        disp(['Sub number ' num2str(numSubs(sub)) ' is an outlier in accuracy'])
    end
    end
else
    accSTDMat = BAP_RT_acc(:,:,2); accSTDMat = accSTDMat(:); %calculate std over all conditions
    upLim = mean(accSTDMat) + 3*std(accSTDMat);
    downLim = mean(accSTDMat) - 3*std(accSTDMat);
for ld = 1:4
    for sub = 1:length(numSubs)
    if BAP_RT_acc(sub,ld,2) >= upLim  || BAP_RT_acc(sub,ld,2) <= downLim
        disp(['Sub number ' num2str(numSubs(sub)) ' is an outlier in accuracy'])
    end
    end
end
end
%% Statistical analysis
if append == 1
    rhoRT = zeros(7,1);    pRT = zeros(7,1);
    rhoAcc = zeros(7,1);   pAcc = zeros(7,1);
    rhoIES = zeros(7,1);   pIES = zeros(7,1);
    pearRRT = cell(7,1);   pearPRT = cell(7,1);
    RL_RT = cell(7,1);     RU_RT = cell(7,1);
    pearRAcc = cell(7,1);  pearPAcc = cell(7,1);
    RL_acc = cell(7,1);    RU_acc = cell(7,1);
    pearRIES = cell(7,4) ;  pearPIES = cell(7,4);
    RL_IES = cell(7,4);    RU_IES = cell(7,4);
else
   rhoRT_ld = zeros(7,4);  pRT_ld = zeros(7,4);
   rhoAcc_ld = zeros(7,4); pAcc_ld = zeros(7,4);
   rhoIES_ld = zeros(7,4); pIES_ld = zeros(7,4);
   pearRRT_ld = cell(7,4); pearPRT_ld = cell(7,4);
   RL_RT_ld = cell(7,4);   RU_RT_ld = cell(7,4);
   pearRacc_ld = cell(7,4); pearPacc_ld = cell(7,4);
   RL_acc_ld = cell(7,4);  RU_acc_ld = cell(7,4);
   pearRIES_ld = cell(7,4) ; pearPIES_ld = cell(7,4);
   RL_IES_ld = cell(7,4);   RU_IES_ld = cell(7,4);
end
structures={'Thal','Caud','Puta','Pall','Hipp','Amyg','Accu'};

for subStr=1:7
    if append == 1
    %Correlation between LV and behavior
    %RT
    [rhoRT(subStr),pRT(subStr)] = corr(LV(:,subStr),BAP_RT_acc(:,1),'Type','Spearman');
    %Accuracy
    [rhoAcc(subStr),pAcc(subStr)] = corr(LV(:,subStr),BAP_RT_acc(:,2),'Type','Spearman');
    %IES
    [rhoIES(subStr),pIES(subStr)] = corr(LV(:,subStr),BAP_RT_acc(:,3),'Type','Spearman');
    %Corrleation with CI
    [pearRRT{subStr},pearPRT{subStr},RL_RT{subStr},RU_RT{subStr}] = corrcoef(LV(:,subStr),BAP_RT_acc(:,1));
    [pearRAcc{subStr},pearPAcc{subStr},RL_acc{subStr},RU_acc{subStr}] = corrcoef(LV(:,subStr),BAP_RT_acc(:,2));
    [pearRIES{subStr},pearPIES{subStr},RL_IES{subStr},RU_IES{subStr}] = corrcoef(LV(:,subStr),BAP_RT_acc(:,3));
    elseif append == 2
        %Correlation between LV and behavior in load conditions
        for ld = 1:4
            %RT
            [rhoRT_ld(subStr,ld),pRT_ld(subStr,ld)] = corr(LV(:,subStr),BAP_RT_acc(:,ld,1),'Type','Spearman');
            %Accuracy
            [rhoAcc_ld(subStr,ld),pAcc_ld(subStr,ld)] = corr(LV(:,subStr),BAP_RT_acc(:,ld,2),'Type','Spearman');
            %IES
            [rhoIES_ld(subStr,ld),pIES_ld(subStr,ld)] = corr(LV(:,subStr),BAP_RT_acc(:,ld,3),'Type','Spearman');
            %Corrleation with CI
            [pearRRT_ld{subStr,ld},pearPRT_ld{subStr,ld},RL_RT_ld{subStr,ld},RU_RT_ld{subStr,ld}] = corrcoef(LV(:,subStr),BAP_RT_acc(:,ld,1));
            [pearRacc_ld{subStr,ld},pearPacc_ld{subStr,ld},RL_acc_ld{subStr,ld},RU_acc_ld{subStr,ld}] = corrcoef(LV(:,subStr),BAP_RT_acc(:,ld,2));
            [pearRIES_ld{subStr,ld},pearPIES_ld{subStr,ld},RL_IES_ld{subStr,ld},RU_IES_ld{subStr,ld}] = corrcoef(LV(:,subStr),BAP_RT_acc(:,ld,2));

        end
    end
end



