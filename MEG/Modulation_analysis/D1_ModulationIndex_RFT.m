function []=D1_ModulationIndex_RFT(subj,system,ev_ind,grp_indiv,append,plot)
%[]=B1_ModulationIndex_RFT(subj,system,ev_ind,grp_indiv,append,plot)
%Calculates modulation index for rapid frequency tagging power
%inputs:
% subj: code of subject you want to analyse
% system: code of system you want to run the program on (1-> portal,
% 2-> office pc)
% ev_ind: 1-> run it on evoked data, 2-> run it on induced data
% grp_indiv: 1-> group ROI, 2-> individual ROIs
% append: 1-> append all load conditions 2-> analyse load conditions
% separately
% plot: 1-> plot MI over time and LI, 2-> do not plot


if nargin < 6
    if nargin < 1
        subj = input('which subject:');
        system = input('On which system you want to run the program (1/2)?');
        ev_ind = input('evoked or induced (1/2)?');
        grp_indiv = input('group ROI or individual RO (1/2)I?');
        append = input('append load conditions (1/2)?');
        plot = input('Do you want to plot everything (1/2)?');
    elseif nargin < 2
        system = input('On which system you want to run the program?');
        ev_ind = input('evoked or induced?');
        grp_indiv = input('group ROI or individual ROI?');
        append = input('append load conditions (1/2)?');
        plot = input('Do you want to plot everything (1/2)?');
    elseif nargin < 3
        ev_ind = input('evoked or induced?');
        grp_indiv = input('group ROI or individual ROI?');
        append = input('append load conditions (1/2)?');
        plot = input('Do you want to plot everything (1/2)?');
    elseif nargin < 4
        grp_indiv = input('group ROI or individual ROI?');
        append = input('append load conditions (1/2)?');
        plot = input('Do you want to plot everything (1/2)?');
    elseif nargin < 5
        append = input('append load conditions (1/2)?');
        plot = input('Do you want to plot everything (1/2)?');
    else
        plot = input('Do you want to plot everything (1/2)?');
    end
end


%% Load clean MEG file
switch system
    case 1 %portal
        data_folder='/rds/projects/j/jenseno-avtemporal-attention/Load/MEG Data/proc_data/'; %Portal
        addpath /rds/projects/j/jenseno-avtemporal-attention/MATLAB/fieldtrip-20210328 %Portal
        roi_loadPath = '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/group_level/RFT/ROI_RFT/';
        saveFolderMat = '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/indiv_level/Modulation_Index/';
        saveFolderFig = '/rds/projects/j/jenseno-avtemporal-attention/Load/Results/FieldTripPlots/';
    case 2 %office PC
        data_folder='Z:\Load\MEG Data\proc_data\'; %Windows
        addpath Z:\MATLAB\fieldtrip-20210328 %Windows
        roi_loadPath = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\RFT\ROI_RFT\';
        saveFolderMat = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\indiv_level\Modulation_Index\';
        saveFolderFig = 'Z:\Load\Results\FieldTripPlots\';
end

ft_defaults

% Load appropriate roi
if grp_indiv == 1
    load([roi_loadPath 'ROI_L_sens'])
    load([roi_loadPath 'ROI_R_sens'])
elseif grp_indiv ==  2
    load([roi_loadPath 'ROI_L_sens_indiv'])
    load([roi_loadPath 'ROI_R_sens_indiv'])
    load([roi_loadPath 'selected_indiv_ROIs'])
end


%% frequency analysis and TFRs
if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end

disp(['loading ' sub])
if ev_ind == 1
    load([data_folder sub filesep sub '_TFR_RFT_correct_only.mat']);fprintf('Done\n')
else
    load([data_folder sub filesep sub '_TFR_RFT_dt_correct_only.mat']);fprintf('Done\n')
end

%% Channel selection
switch ev_ind
    case 1
        MEG_sens  = strmatch('MEG',TFR.left.RFT.ev.f1{1,1}.label);
        sens_type = str2num(cellfun(@(x) x(end),TFR.left.RFT.ev.f1{1,1}.label(MEG_sens),'UniformOutput',1));
    case 2
        MEG_sens  = strmatch('MEG',TFR.left.RFT.ind.f1{1,1}.label);
        sens_type = str2num(cellfun(@(x) x(end),TFR.left.RFT.ind.f1{1,1}.label(MEG_sens),'UniformOutput',1));
end

planars = sort([MEG_sens(sens_type==2) ; MEG_sens(sens_type==3)]);
mags    = MEG_sens(sens_type==1);

%% Rename/Separate files
if append == 2
    for ld = 1:4
        cue_right_63_flicker_right{ld}  = TFR.right.RFT.ind.f1{2,ld}; %-->ROIs on left
        dist_right_63_flicker_right{ld} = TFR.left.RFT.ind.f1{2,ld};  %-->ROIs on left
        cue_left_63_flicker_left{ld}  = TFR.left.RFT.ind.f1{1,ld};    %-->ROIs on right
        dist_left_63_flicker_left{ld} = TFR.right.RFT.ind.f1{1,ld};   %-->ROIs on right
        
        cue_right_70_flicker_right{ld}  = TFR.right.RFT.ind.f2{1,ld}; %-->ROIs on left
        dist_right_70_flicker_right{ld} = TFR.left.RFT.ind.f2{1,ld};  %-->ROIs on left
        cue_left_70_flicker_left{ld}  = TFR.left.RFT.ind.f2{2,ld};    %-->ROIs on right
        dist_left_70_flicker_left{ld} = TFR.right.RFT.ind.f2{2,ld};   %-->ROIs on right
    end
elseif append ==1
    
    switch ev_ind
        case 1 %using evoked data
            cue_right_63_flicker_right  = TFR.attRight.ev.f1{2}; %-->ROIs on left
            dist_right_63_flicker_right = TFR.attLeft.ev.f1{2};  %-->ROIs on left
            cue_left_63_flicker_left  = TFR.attLeft.ev.f1{1};    %-->ROIs on right
            dist_left_63_flicker_left = TFR.attRight.ev.f1{1};   %-->ROIs on right
            
            cue_right_70_flicker_right  = TFR.attRight.ev.f2{1}; %-->ROIs on left
            dist_right_70_flicker_right = TFR.attLeft.ev.f2{1};  %-->ROIs on left
            cue_left_70_flicker_left  = TFR.attLeft.ev.f2{2};    %-->ROIs on right
            dist_left_70_flicker_left = TFR.attRight.ev.f2{2};   %-->ROIs on right
        case 2 %using induced data
            cue_right_63_flicker_right  = TFR.attRight.ind.f1{2}; %-->ROIs on left
            dist_right_63_flicker_right = TFR.attLeft.ind.f1{2};  %-->ROIs on left
            cue_left_63_flicker_left  = TFR.attLeft.ind.f1{1};    %-->ROIs on right
            dist_left_63_flicker_left = TFR.attRight.ind.f1{1};   %-->ROIs on right
            
            cue_right_70_flicker_right  = TFR.attRight.ind.f2{1}; %-->ROIs on left
            dist_right_70_flicker_right = TFR.attLeft.ind.f2{1};  %-->ROIs on left
            cue_left_70_flicker_left  = TFR.attLeft.ind.f2{2};    %-->ROIs on right
            dist_left_70_flicker_left = TFR.attRight.ind.f2{2};   %-->ROIs on right
    end
end
%% Calculate and plot MI over time

%find last non-nan datapoint
if append == 2
    t_end = TFR.left.RFT.ind.f1{1, 1}.time(find(~isnan(squeeze(TFR.left.RFT.ind.f1{1, 1}.powspctrm(1,1,1,:))),1,'last'));
else
    if ev_ind == 1
        t_end = TFR.attLeft.ev.f1{1, 1}.time(find(~isnan(squeeze(TFR.attLeft.ev.f1{1, 1}.powspctrm(1,1,1,:))),1,'last'));
    else
        t_end = TFR.attLeft.ind.f1{1, 1}.time(find(~isnan(squeeze(TFR.attLeft.ind.f1{1, 1}.powspctrm(1,1,1,:))),1,'last'));
    end
end
cfg = [];
cfg.latency     = [-.85 .15]; %=>target locked %for cue locked=>[0 t_end]; %the time of interest
cfg.avgoverfreq = 'yes';
if append == 1; cfg.avgoverrpt  = 'yes'; end
cfg.avgoverchan = 'yes';
cfg.nanmean     = 'yes';

%% right RFT MI on left ROIs
disp('Calculating MIs...')
%%%%% on left sensors %%%%%
switch grp_indiv
    case 1 %on group ROIs
        cfg.channel = {ROI_labels_left_sens{:}};
    case 2 %on individual ROIs
        if chos_roi_L(subj,:)>5; chos_roi_L(subj,:)=chos_roi_L(subj,:)/100; end
        cfg.channel = {ROI_all_subs_left_sens_lbl{chos_roi_L(subj,:),subj}};
end
%%%%% right 63 right att- right dist %%%%%
cfg.frequency = [63 63];
if append == 2
    for ld = 1:4
        cue_right_63_right_overTime_L{ld}  = ft_selectdata(cfg,cue_right_63_flicker_right{ld});
        dist_right_63_right_overTime_L{ld} = ft_selectdata(cfg,dist_right_63_flicker_right{ld});
        
        cue_right_63_right_overTime_pow_L(:,ld)  = squeeze(cue_right_63_right_overTime_L{ld}.powspctrm); %separately run for each hemisphere's ROI always with R-L
        dist_right_63_right_overTime_pow_L(:,ld) = squeeze(dist_right_63_right_overTime_L{ld}.powspctrm);
        
        MI_overTime_63_right_left_sens(:,ld) = (cue_right_63_right_overTime_pow_L(:,ld) - dist_right_63_right_overTime_pow_L(:,ld))...
            ./(cue_right_63_right_overTime_pow_L(:,ld) + dist_right_63_right_overTime_pow_L(:,ld));
    end
else
    cue_right_63_right_overTime_L  = ft_selectdata(cfg,cue_right_63_flicker_right);
    dist_right_63_right_overTime_L = ft_selectdata(cfg,dist_right_63_flicker_right);
    
    cue_right_63_right_overTime_pow_L  = squeeze(cue_right_63_right_overTime_L.powspctrm); %separately run for each hemisphere's ROI always with R-L
    dist_right_63_right_overTime_pow_L = squeeze(dist_right_63_right_overTime_L.powspctrm);
    
    MI_overTime_63_right_left_sens = (cue_right_63_right_overTime_pow_L - dist_right_63_right_overTime_pow_L)...
        ./(cue_right_63_right_overTime_pow_L + dist_right_63_right_overTime_pow_L);
end
%%%%% right 70 right att- right dist %%%%%
cfg.frequency = [70 70];
if append == 2
    for ld = 1:4
        cue_right_70_right_overTime_L{ld}  = ft_selectdata(cfg,cue_right_70_flicker_right{ld});
        dist_right_70_right_overTime_L{ld} = ft_selectdata(cfg,dist_right_70_flicker_right{ld});
        
        cue_right_70_right_overTime_pow_L(:,ld)  = squeeze(cue_right_70_right_overTime_L{ld}.powspctrm);
        dist_right_70_right_overTime_pow_L(:,ld) = squeeze(dist_right_70_right_overTime_L{ld}.powspctrm);
        
        MI_overTime_70_right_left_sens(:,ld) = (cue_right_70_right_overTime_pow_L(:,ld) - dist_right_70_right_overTime_pow_L(:,ld))...
            ./(cue_right_70_right_overTime_pow_L(:,ld) + dist_right_70_right_overTime_pow_L(:,ld));
    end
else
    cue_right_70_right_overTime_L  = ft_selectdata(cfg,cue_right_70_flicker_right);
    dist_right_70_right_overTime_L = ft_selectdata(cfg,dist_right_70_flicker_right);
    
    cue_right_70_right_overTime_pow_L  = squeeze(cue_right_70_right_overTime_L.powspctrm);
    dist_right_70_right_overTime_pow_L = squeeze(dist_right_70_right_overTime_L.powspctrm);
    
    MI_overTime_70_right_left_sens = (cue_right_70_right_overTime_pow_L - dist_right_70_right_overTime_pow_L)...
        ./(cue_right_70_right_overTime_pow_L + dist_right_70_right_overTime_pow_L);
end
%calculate MI-RFT right (MI_63_right + MI_70_right) on left sensors
if append == 2
    for ld = 1:4
        MI_overTime_RFT_right_left_sens(:,ld) = mean([MI_overTime_63_right_left_sens(:,ld),MI_overTime_70_right_left_sens(:,ld)],2);
    end
else
    MI_overTime_RFT_right_left_sens = mean([MI_overTime_63_right_left_sens,MI_overTime_70_right_left_sens],2);
end
%plot MI over time right RFT on left sensors
if plot == 1
    figure(1);
    if append == 2
        for ld = 1:4
            subplot(2,2,ld);
            hold on;
            plot(dist_right_70_right_overTime_L{ld}.time, MI_overTime_RFT_right_left_sens(:,ld), 'k', 'LineWidth', 2);
            
            xlabel('Time'); ylabel('Modulation Index'); title (['MI L-ROIs R-RFT - Sub' num2str(subj) ' Cond:' num2str(ld)])
            xlim(cfg.latency);
            line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;
        end
    else
        hold on;
        plot(dist_right_70_right_overTime_L.time, MI_overTime_RFT_right_left_sens, 'k', 'LineWidth', 2);
        
        xlabel('Time'); ylabel('Modulation Index'); title (['MI over L-ROIs in right RFT - Sub' num2str(subj)])
        xlim(cfg.latency);
        line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;
    end
end
%% left RFT MI on right ROIs

%%%%% on right sensors %%%%%
switch grp_indiv
    case 1 %on group ROIs
        cfg.channel = {ROI_labels_right_sens{:}};
    case 2 %on individual ROIs
        if chos_roi_R(subj,:)>5; chos_roi_R(subj,:)=chos_roi_R(subj,:)/100; end
        cfg.channel = {ROI_all_subs_right_sens_lbl{chos_roi_R(subj,:),subj}};
end

%%%%% left 63 left att- left dist %%%%%
cfg.frequency = [63 63];
if append == 2
    for ld = 1:4
        cue_left_63_left_overTime_R{ld}  = ft_selectdata(cfg,cue_left_63_flicker_left{ld});
        dist_left_63_left_overTime_R{ld} = ft_selectdata(cfg,dist_left_63_flicker_left{ld});
        
        cue_left_63_left_overTime_pow_R(:,ld) = squeeze(cue_left_63_left_overTime_R{ld}.powspctrm); %separately run for each hemisphere's ROI always with R-L
        dist_left_63_left_overTime_pow_R(:,ld) = squeeze(dist_left_63_left_overTime_R{ld}.powspctrm);
        %for left we do dist - cue to make it similar to alpha
        MI_overTime_63_left_right_sens(:,ld) = (dist_left_63_left_overTime_pow_R(:,ld) - cue_left_63_left_overTime_pow_R(:,ld))...
            ./(dist_left_63_left_overTime_pow_R(:,ld) + cue_left_63_left_overTime_pow_R(:,ld));
    end
else
    cue_left_63_left_overTime_R  = ft_selectdata(cfg,cue_left_63_flicker_left);
    dist_left_63_left_overTime_R = ft_selectdata(cfg,dist_left_63_flicker_left);
    
    cue_left_63_left_overTime_pow_R  = squeeze(cue_left_63_left_overTime_R.powspctrm); %separately run for each hemisphere's ROI always with R-L
    dist_left_63_left_overTime_pow_R = squeeze(dist_left_63_left_overTime_R.powspctrm);
    %for left we do dist - cue to make it similar to alpha
    MI_overTime_63_left_right_sens = (dist_left_63_left_overTime_pow_R - cue_left_63_left_overTime_pow_R)...
        ./(dist_left_63_left_overTime_pow_R + cue_left_63_left_overTime_pow_R);
end
%%%%%  left 70 left att- left dist %%%%%
cfg.frequency = [70 70];
if append == 2
    for ld = 1:4
        cue_left_70_left_overTime_R{ld}  = ft_selectdata(cfg,cue_left_70_flicker_left{ld});
        dist_left_70_left_overTime_R{ld} = ft_selectdata(cfg,dist_left_70_flicker_left{ld});
        
        cue_left_70_left_overTime_pow_R(:,ld)  = squeeze(cue_left_70_left_overTime_R{ld}.powspctrm);
        dist_left_70_left_overTime_pow_R(:,ld) = squeeze(dist_left_70_left_overTime_R{ld}.powspctrm);
        
        MI_overTime_70_left_right_sens(:,ld) = (dist_left_70_left_overTime_pow_R(:,ld) - cue_left_70_left_overTime_pow_R(:,ld))...
            ./(dist_left_70_left_overTime_pow_R(:,ld) + cue_left_70_left_overTime_pow_R(:,ld));
    end
else
    cue_left_70_left_overTime_R  = ft_selectdata(cfg,cue_left_70_flicker_left);
    dist_left_70_left_overTime_R = ft_selectdata(cfg,dist_left_70_flicker_left);
    
    cue_left_70_left_overTime_pow_R  = squeeze(cue_left_70_left_overTime_R.powspctrm);
    dist_left_70_left_overTime_pow_R = squeeze(dist_left_70_left_overTime_R.powspctrm);
    
    MI_overTime_70_left_right_sens = (dist_left_70_left_overTime_pow_R - cue_left_70_left_overTime_pow_R)...
        ./(dist_left_70_left_overTime_pow_R + cue_left_70_left_overTime_pow_R);
end
%calculate MI-RFT left (MI_63_right + MI_70_right) on right sensors
if append == 2
    for ld = 1:4
        MI_overTime_RFT_left_right_sens(:,ld) = mean([MI_overTime_63_left_right_sens(:,ld),MI_overTime_70_left_right_sens(:,ld)],2);
    end
else
    MI_overTime_RFT_left_right_sens = mean([MI_overTime_63_left_right_sens,MI_overTime_70_left_right_sens],2);
end
%plot MI over time left RFT on right sensors
if plot == 1
    figure(2);
    if append == 2
        for ld = 1:4
            subplot(2,2,ld);
            hold on;
            plot(dist_left_70_left_overTime_R{ld}.time, MI_overTime_RFT_left_right_sens(:,ld), 'k', 'LineWidth', 2);
            
            xlabel('Time'); ylabel('Modulation Index'); title (['MI R-ROIs L-RFT - Sub' num2str(subj) ' Cond:' num2str(ld)])
            xlim(cfg.latency);
            line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;
        end
    else
        hold on;
        plot(dist_left_70_left_overTime_R.time, MI_overTime_RFT_left_right_sens, 'k', 'LineWidth', 2);
        
        xlabel('Time'); ylabel('Modulation Index'); title (['MI over R-ROIs in left RFT - Sub' num2str(subj)])
        xlim(cfg.latency);
        line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;
    end
end
%% Calculate and bar plot MI -- right RFT on left sensors left RFT on right sensors
if append == 1 %all this section for append == 1 can be changed to how append == 2 was calculated
    cfg.avgovertime = 'yes';
    
    %%%%% left sensors, right RFT %%%%
    switch grp_indiv
        case 1 %on group ROIs
            cfg.channel = {ROI_labels_left_sens{:}};
        case 2 %on individual ROIs
            if chos_roi_L(subj,:)>5; chos_roi_L(subj,:)=chos_roi_L(subj,:)/100; end
            cfg.channel = {ROI_all_subs_left_sens_lbl{chos_roi_L(subj,:),subj}};
    end
    
    %%%%% right 63 right att- right dist %%%%%
    cfg.frequency = [63 63];
    
    cue_right_63_right_avgTime_L  = ft_selectdata(cfg,cue_right_63_flicker_right);
    dist_right_63_right_avgTime_L = ft_selectdata(cfg,dist_right_63_flicker_right);
    
    cue_right_63_right_avgTime_pow_L  = squeeze(cue_right_63_right_avgTime_L.powspctrm); %separately run for each hemisphere's ROI always with R-L
    dist_right_63_right_avgTime_pow_L = squeeze(dist_right_63_right_avgTime_L.powspctrm);
    
    MI_avgTime_63_right_left_sens = (cue_right_63_right_avgTime_pow_L - dist_right_63_right_avgTime_pow_L)...
        ./(cue_right_63_right_avgTime_pow_L + dist_right_63_right_avgTime_pow_L);
    
    %%%%% right 70 right att- right dist %%%%%
    cfg.frequency = [70 70];
    
    cue_right_70_right_avgTime_L  = ft_selectdata(cfg,cue_right_70_flicker_right);
    dist_right_70_right_avgTime_L = ft_selectdata(cfg,dist_right_70_flicker_right);
    
    cue_right_70_right_avgTime_pow_L  = squeeze(cue_right_70_right_avgTime_L.powspctrm);
    dist_right_70_right_avgTime_pow_L = squeeze(dist_right_70_right_avgTime_L.powspctrm);
    
    MI_avgTime_70_right_left_sens = (cue_right_70_right_avgTime_pow_L - dist_right_70_right_avgTime_pow_L)...
        ./(cue_right_70_right_avgTime_pow_L + dist_right_70_right_avgTime_pow_L);
    
    % Calculate avg MI-RFT right (MI_63_right + MI_70_right) on left sensors
    MI_avgTime_RFT_right_left_sens = mean([MI_avgTime_63_right_left_sens,MI_avgTime_70_right_left_sens],2);
    
    %%%%% right sensors, left RFT %%%%%
    switch grp_indiv
        case 1 %on group ROIs
            cfg.channel = {ROI_labels_right_sens{:}};
        case 2 %on individual ROIs
            if chos_roi_R(subj,:)>5; chos_roi_R(subj,:)=chos_roi_R(subj,:)/100; end
            cfg.channel = {ROI_all_subs_right_sens_lbl{chos_roi_R(subj,:),subj}};
    end
    
    %%%%% left 63 left att- left dist %%%%%
    cfg.frequency = [63 63];
    
    cue_left_63_left_avgTime_R  = ft_selectdata(cfg,cue_left_63_flicker_left);
    dist_left_63_left_avgTime_R = ft_selectdata(cfg,dist_left_63_flicker_left);
    
    cue_left_63_left_avgTime_pow_R  = squeeze(cue_left_63_left_avgTime_R.powspctrm); %separately run for each hemisphere's ROI always with R-L
    dist_left_63_left_avgTime_pow_R = squeeze(dist_left_63_left_avgTime_R.powspctrm);
    
    MI_avgTime_63_left_right_sens = (dist_left_63_left_avgTime_pow_R - cue_left_63_left_avgTime_pow_R)...
        ./(dist_left_63_left_avgTime_pow_R + cue_left_63_left_avgTime_pow_R);
    
    %%%%%  left 70 left att- left dist %%%%%
    cfg.frequency = [70 70];
    
    cue_left_70_left_avgTime_R  = ft_selectdata(cfg,cue_left_70_flicker_left);
    dist_left_70_left_avgTime_R = ft_selectdata(cfg,dist_left_70_flicker_left);
    
    cue_left_70_left_avgTime_pow_R  = squeeze(cue_left_70_left_avgTime_R.powspctrm);
    dist_left_70_left_avgTime_pow_R = squeeze(dist_left_70_left_avgTime_R.powspctrm);
    
    MI_avgTime_70_left_right_sens = (dist_left_70_left_avgTime_pow_R - cue_left_70_left_avgTime_pow_R)...
        ./(dist_left_70_left_avgTime_pow_R + cue_left_70_left_avgTime_pow_R);
    
    %calculate MI-RFT left (MI_63_right + MI_70_right) on right sensors
    MI_avgTime_RFT_left_right_sens = mean([MI_avgTime_63_left_right_sens,MI_avgTime_70_left_right_sens],2);
    if plot == 1
        figure(3);
        bar([1,2],[MI_avgTime_RFT_right_left_sens,MI_avgTime_RFT_left_right_sens])
        set(gca,'XTicklabel',{'L__ROI','R__ROI'},'XTickLabelRotation',45); ylabel('MI');
    end
elseif append == 2 %this basically outputs the same results as selecting data with cfg.avgovertime = 'yes'
    MI_avgTime_63_right_left_sens = nanmean(MI_overTime_63_right_left_sens,1);
    MI_avgTime_70_right_left_sens = nanmean(MI_overTime_63_right_left_sens,1);
    MI_avgTime_RFT_right_left_sens = mean([MI_avgTime_63_right_left_sens;MI_avgTime_70_right_left_sens],1);
    MI_avgTime_63_left_right_sens = nanmean(MI_overTime_63_left_right_sens,1);
    MI_avgTime_70_left_right_sens = nanmean(MI_overTime_70_left_right_sens,1);
    MI_avgTime_RFT_left_right_sens = mean([MI_avgTime_63_left_right_sens;MI_avgTime_70_left_right_sens],1);
    if plot == 1
        figure(3)
        for ld = 1:4
            subplot(2,2,ld)
            bar([1,2],[MI_avgTime_RFT_right_left_sens(ld),MI_avgTime_RFT_left_right_sens(ld)])
            set(gca,'XTicklabel',{'L__ROI','R__ROI'},'XTickLabelRotation',45); ylabel('MI');
            title(['RFT_LI in Cond:' num2str(ld)])
        end
    end
end
%% Save necessary data
clear TFR_trials TFR
disp('saving...')
if append == 2
    save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'RFT_dt_MI_Data_indiv_ROI_load'],'MI_overTime_RFT_left_right_sens','MI_overTime_RFT_right_left_sens',...
        'MI_avgTime_RFT_left_right_sens','MI_avgTime_RFT_right_left_sens');
    save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'RFT_dt_All_Data_indiv_ROI_load'],'MI_avgTime_70_left_right_sens','MI_avgTime_70_right_left_sens',...
        'MI_avgTime_63_right_left_sens','MI_avgTime_63_left_right_sens','MI_overTime_70_left_right_sens','MI_overTime_70_right_left_sens',...
        'MI_overTime_63_right_left_sens','MI_overTime_63_left_right_sens');
    if plot == 1
        saveas(figure(1),[saveFolderFig filesep num2str(subj) '_MIoverTime_RFT_dt_indiv_ROI_load.jpg']);
        saveas(figure(2),[saveFolderFig filesep num2str(subj) '_MIndex_RFT_dt_indiv_ROI_load.jpg']);
    end
else
if grp_indiv == 2
    save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'RFT_dt_MI_Data_indiv_ROI'],'MI_overTime_RFT_left_right_sens','MI_overTime_RFT_right_left_sens',...
        'MI_avgTime_RFT_left_right_sens','MI_avgTime_RFT_right_left_sens');
    save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'RFT_dt_All_Data_indiv_ROI'],'MI_avgTime_70_left_right_sens','MI_avgTime_70_right_left_sens',...
        'MI_avgTime_63_right_left_sens','MI_avgTime_63_left_right_sens','MI_overTime_70_left_right_sens','MI_overTime_70_right_left_sens',...
        'MI_overTime_63_right_left_sens','MI_overTime_63_left_right_sens');
    if plot == 1
        saveas(figure(1),[saveFolderFig filesep num2str(subj) '_MIoverTime_RFT_dt_indiv_ROI.jpg']);
        saveas(figure(2),[saveFolderFig filesep num2str(subj) '_MIndex_RFT_dt_indiv_ROI.jpg']);
    end
else
    if ev_ind == 2
        save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'RFT_dt_MI_Data'],'MI_overTime_RFT_left_right_sens','MI_overTime_RFT_right_left_sens',...
            'MI_avgTime_RFT_left_right_sens','MI_avgTime_RFT_right_left_sens');
        save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'RFT_dt_All_Data'],'MI_avgTime_70_left_right_sens','MI_avgTime_70_right_left_sens',...
            'MI_avgTime_63_right_left_sens','MI_avgTime_63_left_right_sens','MI_overTime_70_left_right_sens','MI_overTime_70_right_left_sens',...
            'MI_overTime_63_right_left_sens','MI_overTime_63_left_right_sens');
        if plot == 1
            saveas(figure(1),[saveFolderFig filesep num2str(subj) '_MIoverTime_RFT_dt.jpg']);
            saveas(figure(2),[saveFolderFig filesep num2str(subj) '_MIndex_RFT_dt.jpg']);
        end
    elseif ev_ind == 1
        save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'RFT_MI_Data'],'MI_overTime_RFT_left_right_sens','MI_overTime_RFT_right_left_sens',...
            'MI_avgTime_RFT_left_right_sens','MI_avgTime_RFT_right_left_sens');
        save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'RFT_All_Data'],'MI_avgTime_70_left_right_sens','MI_avgTime_70_right_left_sens',...
            'MI_avgTime_63_right_left_sens','MI_avgTime_63_left_right_sens','MI_overTime_70_left_right_sens','MI_overTime_70_right_left_sens',...
            'MI_overTime_63_right_left_sens','MI_overTime_63_left_right_sens');
        if plot == 1
            saveas(figure(1),[saveFolderFig filesep num2str(subj) '_MIoverTime_RFT.jpg']);
            saveas(figure(2),[saveFolderFig filesep num2str(subj) '_MIndex_RFT.jpg']);
        end
    end
end
end
close all
