%load/select all subjs TFR
clear;clc;clf;close all

%% Load clean MEG file
data_folder='/rds/projects/j/jenseno-avtemporal-attention/Projects/Load/MEG Data/proc_data/'; %Portal
addpath /rds/projects/j/jenseno-avtemporal-attention/Programming/MATLAB/fieldtrip-20210328 %Portal
saveFolder = '/rds/projects/j/jenseno-avtemporal-attention/Programming/MATLAB/Perceptual_Load/FieldTrip/Results/group_level/Alpha';
load([saveFolder '/ROI_alpha/ROI_dt_right_sym.mat']);
load([saveFolder '/ROI_alpha/ROI_dt_left_sym.mat']);

% saveFolder = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Alpha'; %Windows
% data_folder='Z:\Load\MEG Data\proc_data\'; %Windows
% addpath Z:\MATLAB\fieldtrip-20210328 %Windows

ft_defaults

badSubs = [23,28];
numSub = setxor(1:35,badSubs);
TFR_attRight_alpha_all_subs_ROI_R = cell(length(numSub),1); %for all TFR
TFR_attLeft_alpha_all_subs_ROI_R  = cell(length(numSub),1);
TFR_attRight_alpha_all_subs_ROI_L = cell(length(numSub),1); %for all TFR
TFR_attLeft_alpha_all_subs_ROI_L  = cell(length(numSub),1);
TFR_dt_allsubs_LF_ROI_R_cntrst = cell(length(numSub),1);
TFR_dt_allsubs_LF_ROI_L_cntrst = cell(length(numSub),1);

for subj=numSub
    if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
    disp(['loading ' sub])
    load([data_folder sub filesep sub '_TFR_LF_dt_correct_only.mat', 'TFR']);fprintf('Done\n')
    
    %% Channel selection
    MEG_sens=strmatch('MEG',TFR.left.LF.ind{1}.label);
    sens_type=str2num(cellfun(@(x) x(end),TFR.left.LF.ind{1}.label(MEG_sens),'UniformOutput',1));
    
    planars=sort([MEG_sens(sens_type==2) ; MEG_sens(sens_type==3)]);
    mags=MEG_sens(sens_type==1);
    
    %% Append all configs and loads 
    TFR.attRight.LF = ft_appendfreq([],TFR.right.LF.ind{:}); 
    TFR.attLeft.LF  = ft_appendfreq([],TFR.left.LF.ind{:});
    
    %% Choose necessary data
    cfg = [];
    cfg.latency     = [-.85 0]; %the time of interest is 850ms before target onset
    cfg.frequency   = [2 30];   %frequency of interest
    cfg.avgoverrpt  = 'yes';
    cfg.nanmean     = 'yes';
    cfg.avgoverchan = 'yes';
    cfg.channel     = {ROI_lbl_R{:}}; %only choose planars for ROI selection
    
%     %select data on right ROI sensors -- for any other analyses
%     TFR_attRight_alpha_ROI_R = ft_selectdata(cfg,TFR.attRight.LF); %separately run for each hemisphere's ROI always with R-L
%     TFR_attLeft_alpha_ROI_R  = ft_selectdata(cfg,TFR.attLeft.LF);
    
    %select data on right ROI sensors -- for TFR plots
    TFR_attRight_LF_ROI_R = ft_selectdata(cfg,TFR.attRight.LF); %separately run for each hemisphere's ROI always with R-L
    TFR_attLeft_LF_ROI_R  = ft_selectdata(cfg,TFR.attLeft.LF);


%     TFR_attRight_alpha_all_subs_ROI_R{subj,1} = TFR_attRight_alpha_ROI_R;
%     TFR_attLeft_alpha_all_subs_ROI_R{subj,1} = TFR_attLeft_alpha_ROI_R;
    
%     TFR_dt_allsubs_alpha_ROI_R_cntrst{subj,1} = TFR_attLeft_alpha_ROI_R;
%     TFR_dt_allsubs_alpha_ROI_R_cntrst{subj,1}.powspctrm = TFR_attRight_alpha_ROI_R.powspctrm - TFR_attLeft_alpha_ROI_R.powspctrm;

    % for TFR plots
    TFR_dt_allsubs_LF_ROI_R_cntrst{subj,1} = TFR_attLeft_LF_ROI_R;
    TFR_dt_allsubs_LF_ROI_R_cntrst{subj,1}.powspctrm = TFR_attRight_LF_ROI_R.powspctrm - TFR_attLeft_LF_ROI_R.powspctrm;

    cfg.channel     = {ROI_lbl_L{:}}; %only choose planars for ROI selection
%     %select data on left ROI sensors -- for any other analyses
%     TFR_attRight_alpha_ROI_L = ft_selectdata(cfg,TFR.attRight.LF); %separately run for each hemisphere's ROI always with R-L
%     TFR_attLeft_alpha_ROI_L  = ft_selectdata(cfg,TFR.attLeft.LF);

    %select data on left ROI sensors -- for TFR plots
    TFR_attRight_LF_ROI_L = ft_selectdata(cfg,TFR.attRight.LF); %separately run for each hemisphere's ROI always with R-L
    TFR_attLeft_LF_ROI_L  = ft_selectdata(cfg,TFR.attLeft.LF);

%     TFR_attRight_alpha_all_subs_ROI_L{subj,1} = TFR_attRight_alpha_ROI_L;
%     TFR_attLeft_alpha_all_subs_ROI_L{subj,1} = TFR_attLeft_alpha_ROI_L;
    
%     TFR_dt_allsubs_alpha_ROI_L_cntrst{subj,1} = TFR_attLeft_alpha_ROI_L;
%     TFR_dt_allsubs_alpha_ROI_L_cntrst{subj,1}.powspctrm = TFR_attRight_alpha_ROI_L.powspctrm - TFR_attLeft_alpha_ROI_L.powspctrm;
    
    % for TFR plots
    TFR_dt_allsubs_LF_ROI_L_cntrst{subj,1} = TFR_attLeft_LF_ROI_L;
    TFR_dt_allsubs_LF_ROI_L_cntrst{subj,1}.powspctrm = TFR_attRight_LF_ROI_L.powspctrm - TFR_attLeft_LF_ROI_L.powspctrm;

end
disp('Saving contrast TFR...')
% save([saveFolder filesep 'TFR_dt_all_subs_LF_cntrst'],'TFR_dt_allsubs_LF_ROI_R_cntrst',...
%     'TFR_dt_allsubs_LF_ROI_L_cntrst','-v7.3')
disp('done')

%% Average over all subjects 
disp('Grand averaging...')  
cfg=[];
frq_grnd_avg_ROI_R = ft_freqgrandaverage(cfg,TFR_dt_allsubs_LF_ROI_R_cntrst{numSub,1});
frq_grnd_avg_ROI_L = ft_freqgrandaverage(cfg,TFR_dt_allsubs_LF_ROI_L_cntrst{numSub,1});

save([saveFolder filesep 'TFR_dt_all_subs_LF_cntrst_grndavg'],'frq_grnd_avg_ROI_R',...
    'frq_grnd_avg_ROI_L','-v7.3')
%% plot TFR 
cfg = [];
cfg.baseline = 'no';
% cfg.baselinetype = 'relative';
cfg.maskstyle = 'saturation';
cfg.xlim = 'maxmin';  %specification of time in s (here time before target onset)
cfg.ylim = 'maxmin'; %frequency specification, e.g. alpha
cfg.zlim = 'maxmin'; %scale specification; [0.4 2]; [-2.5e-27 2.5e-27];
cfg.layout = 'neuromag306all.lay';
cfg.avgoverchan = 'yes';

% cfg.channel = {ROI_lbl_R(:)};

figure
% subplot(1,2,2)
% title('Right ROI')
ft_singleplotTFR(cfg,frq_grnd_avg_ROI_R);
saveas(gcf,'Right_ROI_TFR.jpg')

figure
% cfg.channel = {ROI_lbl_R(:)};
% subplot(1,2,1)
% title('Left ROI')
ft_singleplotTFR(cfg,frq_grnd_avg_ROI_L);
saveas(gcf,'Left_ROI_TFR.jpg')




