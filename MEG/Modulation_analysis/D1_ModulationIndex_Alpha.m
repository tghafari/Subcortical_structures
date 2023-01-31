function [] = D1_ModulationIndex_Alpha(system,subj,append,plot)
% Caculates modulation indices for each subject
% system: 1->portal 2->windows
% subj: number of subjects you want to run the function on
% append: do you want all conditions to be appended (->1) or load conditions separate (->2)
% plot: do you want to plot MI and ALI? 1-> yes 2-> no

%% Load clean MEG file
if nargin < 1
    system = input('Which system do you want to run on?(1:portal/2:windows');
    subj = input('Enter Subject Number:\n');
    append = input('Would you like to append all load conditions?(1:yes/2:no)');
end

switch system
    case 1
        data_folder='/rds/projects/j/jenseno-avtemporal-attention/Load/MEG Data/proc_data/'; %Portal
        addpath /rds/projects/j/jenseno-avtemporal-attention/MATLAB/fieldtrip-20210328 %Portal
        load '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/group_level/Alpha/ROI_alpha/ROI_dt_right_sym.mat'
        load '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/group_level/Alpha/ROI_alpha/ROI_dt_left_sym.mat'
        saveFolderMat = '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/indiv_level/Modulation_Index/';
%         saveFolderFig = '/rds/projects/j/jenseno-avtemporal-attention/Load/Results/FieldTripPlots/';
        ft_defaults
    case 2
            data_folder='Z:\Load\MEG Data\proc_data\'; %Windows
            addpath Z:\MATLAB\fieldtrip-20210328 %Windows
            load 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Alpha\ROI_alpha\ROI_dt_left_sym.mat'
            load 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Alpha\ROI_alpha\ROI_dt_right_sym.mat'
            saveFolderMat = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\indiv_level\Modulation_Index\';
            ft_defaults
end

%% frequency analysis and TFRs
if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
disp(['loading ' sub])
load([data_folder sub filesep sub '_TFR_LF_dt_correct_only.mat']);fprintf('Done\n')
%% Channel selection
MEG_sens=strmatch('MEG',TFR.left.LF.ind{1}.label);
sens_type=str2num(cellfun(@(x) x(end),TFR.left.LF.ind{1}.label(MEG_sens),'UniformOutput',1));

planars=sort([MEG_sens(sens_type==2) ; MEG_sens(sens_type==3)]);
mags=MEG_sens(sens_type==1);
%% Append all configs and loads 
if append == 1
    TFR_trials.attRight.LF = ft_appendfreq([],TFR_trials.right.LF{:,:}); %if you only want to analyse salient distractor use [3,4].
    TFR_trials.attLeft.LF  = ft_appendfreq([],TFR_trials.left.LF{:,:});
    
    %% Calculate and plot MI over time
    cfg = [];
    cfg.latency     = [-.85 .15]; %the time of interest
    cfg.frequency   = [8 13];     %frequency of interest
    cfg.avgoverfreq = 'yes';
    cfg.avgoverchan = 'yes';
    cfg.avgoverrpt  = 'yes';
    % [78,71,75,61]; %R [87,77,94,86]; %L is not correct  %run separately than R ROI
    % idxR(1:num); %Cecilia's ROI (both hemispheres)
    
    %plot for right ROI -- ROI attLeft is on right sensors
    cfg.channel     = {ROI_lbl_R{:}}; %{ROI_attL{:}};
    TFR_attRight_LF_chanSlctd_R = ft_selectdata(cfg,TFR_trials.attRight.LF); %separately run for each hemisphere's ROI always with R-L
    TFR_attLeft_LF_chanSlctd_R  = ft_selectdata(cfg,TFR_trials.attLeft.LF);
    
    pwrRight = squeeze(TFR_attRight_LF_chanSlctd_R.powspctrm); %separately run for each hemisphere's ROI always with R-L
    pwrLeft  = squeeze(TFR_attLeft_LF_chanSlctd_R.powspctrm); %? why is attLeft over right ROI is positive?
    
    MI_overTime_inROI_R = (pwrRight-pwrLeft)./(pwrRight+pwrLeft); %separately run for each hemisphere's ROI always with R-L

    %plot for left ROI -- ROI attRight is on left sensors
    cfg.channel     = {ROI_lbl_L{:}}; %{ROI_attR{:}};
    TFR_attRight_LF_chanSlctd_L = ft_selectdata(cfg,TFR_trials.attRight.LF); %separately run for each hemisphere's ROI always with R-L
    TFR_attLeft_LF_chanSlctd_L  = ft_selectdata(cfg,TFR_trials.attLeft.LF);
    
    pwrRight = squeeze(TFR_attRight_LF_chanSlctd_L.powspctrm);
    pwrLeft  = squeeze(TFR_attLeft_LF_chanSlctd_L.powspctrm);
    
    MI_overTime_inROI_L = (pwrRight-pwrLeft)./(pwrRight+pwrLeft);
    
    % Calculate and bar plot MI- averages over time by nanmean
    MI_R = nanmean(MI_overTime_inROI_R);
    MI_L = nanmean(MI_overTime_inROI_L);
 %% Plotting
    if plot == 1
    figure(1); 
    % on left sensors
    subplot(1,2,1);
    hold on;
    plot(TFR_attRight_LF_chanSlctd_L.time, MI_overTime_inROI_L, 'k', 'LineWidth', 2);
    
    % legend('att_R', 'att_L', 'MI','Location','northwest');
    xlabel('Time'); ylabel('Modulation Index'); title (['MI over L-ROIs in alpha- Sub' num2str(subj)])
    xlim(cfg.latency);
    line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;
    % on right sensors
    subplot(1,2,2);
    hold on;
    plot(TFR_attRight_LF_chanSlctd_R.time, MI_overTime_inROI_R, 'k', 'LineWidth', 2);
    
    xlabel('Time'); ylabel('Modulation Index'); title (['MI over R-ROIs in alpha - Sub' num2str(subj)])
    xlim(cfg.latency);
    line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;

    figure();
    bar([1,2],[MI_L,MI_R])
    xticklabels({'L__sens','R__sens'});
    title('MI in attention right - attention left');
    end
else 
    %% if conditions are separated- only append the frequency tagging configs

    for ld = 1:4
        TFR_trials.attRight.LF{ld} = ft_appendfreq([],TFR_trials.right.LF{[1 2],ld});
        TFR_trials.attLeft.LF{ld} = ft_appendfreq([],TFR_trials.left.LF{[1 2],ld});
    end
    
    cfg = [];
    cfg.latency     = [-.85 .15]; %the time of interest
    cfg.frequency   = [8 13];     %frequency of interest
    cfg.avgoverfreq = 'yes';
    cfg.avgoverchan = 'yes';
    cfg.avgoverrpt  = 'yes';
    
    for ld = 1:4 %load conditions
        %plot for right ROI -- ROI attLeft is on right sensors
        cfg.channel     = {ROI_lbl_R{:}}; %{ROI_attL{:}};
        TFR_attRight_LF_chanSlctd_R{ld} = ft_selectdata(cfg,TFR_trials.attRight.LF{ld}); %separately run for each hemisphere's ROI always with R-L
        TFR_attLeft_LF_chanSlctd_R{ld}  = ft_selectdata(cfg,TFR_trials.attLeft.LF{ld});
        
        pwrRight_R{ld} = squeeze(TFR_attRight_LF_chanSlctd_R{ld}.powspctrm); %separately run for each hemisphere's ROI always with R-L
        pwrLeft_R{ld}  = squeeze(TFR_attLeft_LF_chanSlctd_R{ld}.powspctrm); %? why is attLeft over right ROI is positive?
        MI_overTime_inROI_R{ld} = (pwrRight_R{ld}-pwrLeft_R{ld})./(pwrRight_R{ld}+pwrLeft_R{ld}); %separately run for each hemisphere's ROI always with R-L
        
        %plot for left ROI -- ROI attRight is on left sensors
        cfg.channel     = {ROI_lbl_L{:}}; %{ROI_attR{:}};
        TFR_attRight_LF_chanSlctd_L{ld} = ft_selectdata(cfg,TFR_trials.attRight.LF{ld}); %separately run for each hemisphere's ROI always with R-L
        TFR_attLeft_LF_chanSlctd_L{ld}  = ft_selectdata(cfg,TFR_trials.attLeft.LF{ld});
        
        pwrRight_L{ld} = squeeze(TFR_attRight_LF_chanSlctd_L{ld}.powspctrm);
        pwrLeft_L{ld}  = squeeze(TFR_attLeft_LF_chanSlctd_L{ld}.powspctrm);
        MI_overTime_inROI_L{ld} = (pwrRight_L{ld}-pwrLeft_L{ld})./(pwrRight_L{ld}+pwrLeft_L{ld});
        
        % Calculate and bar plot MI - averages over time by nanmean        
        MI_R(ld) = nanmean(MI_overTime_inROI_R{ld}); 
        MI_L(ld) = nanmean(MI_overTime_inROI_L{ld});
     end
end
%% Save necessary data
save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'MI_Data_dt_sym_salient'],'MI_L','MI_R','MI_overTime_inROI_L','MI_overTime_inROI_R');
save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'All_Data_dt_sym_salient'],'TFR_attRight_LF_chanSlctd_R','TFR_attRight_LF_chanSlctd_L',...
    'TFR_attLeft_LF_chanSlctd_R','TFR_attLeft_LF_chanSlctd_L','MI_L','MI_R','MI_overTime_inROI_L','MI_overTime_inROI_R');
if plot == 1
saveas(figure(1),[saveFolderFig filesep num2str(subj) '_MIoverTime_dt.jpg']);
saveas(figure(2),[saveFolderFig filesep num2str(subj) '_MIndex_dt.jpg']);
end
end