%% ROI analysis
clear;clc;clf;close all

%% Load clean MEG file
% data_folder='/rds/projects/j/jenseno-avtemporal-attention/Load/MEG Data/proc_data/'; %Portal
% addpath /rds/projects/j/jenseno-avtemporal-attention/MATLAB/fieldtrip-20200320 %Portal
% load '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual Load/FieldTrip/TjerksCodes/ROI_alpha_4_correct_only'
% ft_defaults

data_folder='Z:\Load\MEG Data\proc_data\'; %Windows
addpath Z:\MATLAB\fieldtrip-20200320 %Windows
ft_defaults

saveFolder = 'Z:\Load\Results\Group\Alpha';

badSubs=28;
numSub = setxor(1:35,badSubs);
TFR_attRight_alpha_all_subs = cell(length(numSub),1);
TFR_attLeft_alpha_all_subs  = cell(length(numSub),1);


for subj=numSub
    if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
    disp(['loading ' sub])
    load([data_folder sub filesep sub '_TFR_LF_correct_only.mat']);fprintf('Done\n')
    
    %% Channel selection
    MEG_sens=strmatch('MEG',TFR.left.LF.ind{1}.label);
    sens_type=str2num(cellfun(@(x) x(end),TFR.left.LF.ind{1}.label(MEG_sens),'UniformOutput',1));
    
    planars=sort([MEG_sens(sens_type==2) ; MEG_sens(sens_type==3)]);
    mags=MEG_sens(sens_type==1);
    
    %% Append all configs and loads -- run on the output of Tjerk's script
    TFR.attRight.LF = ft_appendfreq([],TFR.right.LF.ind{:}); %does this accurately append attR and attL data?
    TFR.attLeft.LF  = ft_appendfreq([],TFR.left.LF.ind{:});
    
    %% Choose necessary data
    cfg = [];
    cfg.latency     = [.150 1]; %the time of interest
    cfg.frequency   = [8 13];   %frequency of interest
    cfg.avgoverfreq = 'yes';
    cfg.avgoverrpt  = 'yes';
    cfg.avgovertime = 'yes';
    cfg.nanmean     = 'yes';
    cfg.channel     = planars; %only choose planars for ROI selection
    %select data
    TFR_attRight_alpha_planars = ft_selectdata(cfg,TFR.attRight.LF); %separately run for each hemisphere's ROI always with R-L
    TFR_attLeft_alpha_planars  = ft_selectdata(cfg,TFR.attLeft.LF);

    TFR_attRight_alpha_all_subs{subj,1} = TFR_attRight_alpha_planars;
    TFR_attLeft_alpha_all_subs{subj,1}  = TFR_attLeft_alpha_planars;
end
save([saveFolder filesep 'TFR_all_subs_alpha'],'TFR_attRight_alpha_all_subs','TFR_attLeft_alpha_all_subs')

%% Average over all subjects and contrast R vs. L

cfg = [];
frq_grnd_avg_attRight = ft_freqgrandaverage(cfg,TFR_attRight_alpha_all_subs{:,1});
frq_grnd_avg_attLeft  = ft_freqgrandaverage(cfg,TFR_attLeft_alpha_all_subs{:,1});


R_L_frq_contrast = squeeze(frq_grnd_avg_attRight.powspctrm)-squeeze(frq_grnd_avg_attLeft.powspctrm);
[BRL,idxRL] = sortrows(RLContrast,'descend'); % a vector of data ordered from highest to lowest difference in R-L






