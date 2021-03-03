%frequency analysis and TFRs
for subj = setxor(1:35,28)
%     subj=input('Enter Subject Number:\n');
    if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
    disp(['loading ' sub])
    load([data_folder sub filesep sub '_TFR_LF_dt_correct_only.mat']);fprintf('Done\n')
    %% Channel selection
    MEG_sens=strmatch('MEG',TFR.left.LF.ind{1}.label);
    sens_type=str2num(cellfun(@(x) x(end),TFR.left.LF.ind{1}.label(MEG_sens),'UniformOutput',1));
    
    planars=sort([MEG_sens(sens_type==2) ; MEG_sens(sens_type==3)]);
    mags=MEG_sens(sens_type==1);
    %% Append all configs and loads -- run on the output of Tjerk's script
    %append all configurations togeather
    TFR.attRight.LF = ft_appendfreq([],TFR.right.LF.ind{:}); %does this accurately append attR and attL data?
    TFR.attLeft.LF  = ft_appendfreq([],TFR.left.LF.ind{:});
    
    %% Calculate and plot MI over time
    cfg = [];
    cfg.latency     = [-.85 0]; %the time of interest
    cfg.frequency   = [8 13];   %frequency of interest
    cfg.avgoverfreq = 'yes';
    cfg.avgoverchan = 'yes';
    cfg.avgoverrpt  = 'yes';
    % [78,71,75,61]; %R [87,77,94,86]; %L is not correct  %run separately than R ROI
    % idxR(1:num); %Cecilia's ROI (both hemispheres)
    
    %plot for right ROI -- ROI_left is on right sensors
    cfg.channel     = {ROI_attL{:}};
    TFR_attRight_LF_chanSlctd_R = ft_selectdata(cfg,TFR.attRight.LF); %separately run for each hemisphere's ROI always with R-L
    TFR_attLeft_LF_chanSlctd_R  = ft_selectdata(cfg,TFR.attLeft.LF);
    
    pwrRight = squeeze(TFR_attRight_LF_chanSlctd_R.powspctrm); %separately run for each hemisphere's ROI always with R-L
    pwrLeft  = squeeze(TFR_attLeft_LF_chanSlctd_R.powspctrm); %? why is attLeft over right ROI is positive?
    
    MI_overTime_inROI_R = (pwrRight-pwrLeft)./(pwrRight+pwrLeft); %separately run for each hemisphere's ROI always with R-L
    
    figure(1); subplot(1,2,2);
    hold on;
    plot(TFR_attRight_LF_chanSlctd_R.time, MI_overTime_inROI_R, 'k', 'LineWidth', 2);
    
    xlabel('Time'); ylabel('Modulation Index'); title (['MI over R-ROIs in alpha - Sub' num2str(subj)])
    xlim(cfg.latency);
    line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;
    
    %plot for left ROI -- ROI-right is on left sensors
    cfg.channel     = {ROI_attR{:}};
    TFR_attRight_LF_chanSlctd_L = ft_selectdata(cfg,TFR.attRight.LF); %separately run for each hemisphere's ROI always with R-L
    TFR_attLeft_LF_chanSlctd_L  = ft_selectdata(cfg,TFR.attLeft.LF);
    
    pwrRight = squeeze(TFR_attRight_LF_chanSlctd_L.powspctrm);
    pwrLeft  = squeeze(TFR_attLeft_LF_chanSlctd_L.powspctrm);
    
    MI_overTime_inROI_L = (pwrRight-pwrLeft)./(pwrRight+pwrLeft);
    
    figure(1); subplot(1,2,1);
    hold on;
    plot(TFR_attRight_LF_chanSlctd_L.time, MI_overTime_inROI_L, 'k', 'LineWidth', 2);
    
    % legend('att_R', 'att_L', 'MI','Location','northwest');
    xlabel('Time'); ylabel('Modulation Index'); title (['MI over L-ROIs in alpha- Sub' num2str(subj)])
    xlim(cfg.latency);
    line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;
    
    %% Calculate and bar plot MI
    % cfg.avgovertime = 'yes';
    
    cfg.channel     = {ROI_attL{:}};
    TFR_attRight_LF_R = ft_selectdata(cfg,TFR.attRight.LF);
    TFR_attLeft_LF_R  = ft_selectdata(cfg,TFR.attLeft.LF);
    MI_R = (nanmean(TFR_attRight_LF_R.powspctrm)-nanmean(TFR_attLeft_LF_R.powspctrm))./(nanmean(TFR_attRight_LF_R.powspctrm)+nanmean(TFR_attLeft_LF_R.powspctrm)); %Run after L_ROI
    
    cfg.channel     = {ROI_attR{:}};
    TFR_attRight_LF_L = ft_selectdata(cfg,TFR.attRight.LF);
    TFR_attLeft_LF_L  = ft_selectdata(cfg,TFR.attLeft.LF);
    MI_L = (nanmean(TFR_attRight_LF_L.powspctrm)-nanmean(TFR_attLeft_LF_L.powspctrm))./(nanmean(TFR_attRight_LF_L.powspctrm)+nanmean(TFR_attLeft_LF_L.powspctrm)); %Run after R_ROI
    
    figure(2);
    bar([1,2],[MI_L,MI_R])
    set(gca,'XTicklabel',{'L__ROI','R__ROI'},'XTickLabelRotation',45); ylabel('MI(R-L)');
    %% Save necessary data
    clear TFR_trials TFR
%     mkdir([saveFolderMat filesep 'Sub' num2str(subj)])
    save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'MI_Data_dt'],'MI_L','MI_R','MI_overTime_inROI_L','MI_overTime_inROI_R');
    save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'All_Data_dt'],'TFR_attRight_LF_chanSlctd_R','TFR_attRight_LF_chanSlctd_L',...
        'TFR_attLeft_LF_chanSlctd_R','TFR_attLeft_LF_chanSlctd_L','TFR_attRight_LF_R','TFR_attRight_LF_L',...
        'TFR_attLeft_LF_R','TFR_attLeft_LF_L','MI_L','MI_R','MI_overTime_inROI_L','MI_overTime_inROI_R');
    saveas(figure(1),[saveFolderFig filesep num2str(subj) '_MIoverTime_dt.jpg']);
    saveas(figure(2),[saveFolderFig filesep num2str(subj) '_MIndex_dt.jpg']);
    
    %% clear for next sub
    clear;clc;clf;close all
    
    %% Load clean MEG file
    data_folder='/rds/projects/j/jenseno-avtemporal-attention/Load/MEG Data/proc_data/'; %Portal
    addpath /rds/projects/j/jenseno-avtemporal-attention/MATLAB/fieldtrip-20200320 %Portal
    load('/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/group_level/ROI_alpha/ROI_dt_attRigth.mat')
    load('/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/group_level/ROI_alpha/ROI_dt_attLeft.mat')
    saveFolderMat = '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/indiv_level/Modulation_Index/';
    saveFolderFig = '/rds/projects/j/jenseno-avtemporal-attention/Load/Results/FieldTripPlots/';
    
    % data_folder='Z:\Load\MEG Data\proc_data\'; %Windows
    % addpath Z:\MATLAB\fieldtrip-20200320 %Windows
    % load 'Z:\MATLAB\Perceptual Load\FieldTrip\Results\group_level/ROI_alpha/ROI_attRight'
    % load 'Z:\MATLAB\Perceptual Load\FieldTrip\Results\group_level/ROI_alpha/ROI_attLeftt'
    % saveFolder = 'Z:\MATLAB/Perceptual Load/FieldTrip/Results';
    
    ft_defaults
end