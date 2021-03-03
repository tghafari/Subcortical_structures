function [] = ROISelection_RFT(badSubs,ROI_num)

% ROI analysis -> Calculates ROI for RFT based on grand average
%inputs: 
%   badSubs: subjects you want to reject
%   ROI_num: number of sensors to constitute the ROI for each side & config

%% Load clean MEG file
data_folder='/rds/projects/j/jenseno-avtemporal-attention/Load/MEG Data/proc_data/'; %Portal
addpath /rds/projects/j/jenseno-avtemporal-attention/MATLAB/fieldtrip-20200320 %Portal
% load '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual Load/FieldTrip/TjerksCodes/ROI_alpha_4_correct_only'

% data_folder='Z:\Load\MEG Data\proc_data\'; %Windows
% addpath Z:\MATLAB\fieldtrip-20200320 %Windows

ft_defaults

saveFolder = '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/group_level';
% saveFolder = 'Z:\Load\Results\MATLAB_groupLevels\Alpha';

if nargin < 2
    if nargin < 1
        ROI_num = input('How many sensors in each ROI?');
        badSubs = input('Enter the code of bad subjects:');        
    else
        ROI_num = input('How many sensors in each ROI?');
    end
end

numSub = setxor(1:35,badSubs);

%Preallocate
RFT.c1.group.cue_left_63_left    = cell(length(numSub),1);
RFT.c1.group.dist_left_63_left   = cell(length(numSub),1);
RFT.c2.group.cue_right_63_right  = cell(length(numSub),1);
RFT.c2.group.dist_right_63_right = cell(length(numSub),1);

RFT.c1.group.cue_right_70_right  = cell(length(numSub),1);
RFT.c1.group.dist_right_70_right = cell(length(numSub),1);
RFT.c2.group.cue_left_70_left    = cell(length(numSub),1);
RFT.c2.group.dist_left_70_left   = cell(length(numSub),1);

for subj=numSub
    
    if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
    disp(['loading ' sub])
    load([data_folder sub filesep sub '_TFR_RFT_tl_correct_only.mat']);fprintf('Done\n')
    
    %% Channel selection
    
    MEG_sens  = strmatch('MEG',TFR.left.RFT.ev.f1{1,1}.label);
    sens_type = str2num(cellfun(@(x) x(end),TFR.left.RFT.ev.f1{1,1}.label(MEG_sens),'UniformOutput',1));
    
    planars = sort([MEG_sens(sens_type==2) ; MEG_sens(sens_type==3)]);
    mags    = MEG_sens(sens_type==1);
    
    %% Rename/Separate files
    
    cue_left_63_flicker_left  = TFR.attLeft.ev.f1{1};    %-->ROIs on right
    dist_left_63_flicker_left = TFR.attRight.ev.f1{1};   %-->ROIs on right
    cue_right_63_flicker_right  = TFR.attRight.ev.f1{2}; %-->ROIs on left
    dist_right_63_flicker_right = TFR.attLeft.ev.f1{2};  %-->ROIs on left
    
    cue_right_70_flicker_right  = TFR.attRight.ev.f2{1}; %-->ROIs on left
    dist_right_70_flicker_right = TFR.attLeft.ev.f2{1};  %-->ROIs on left
    cue_left_70_flicker_left  = TFR.attLeft.ev.f2{2};    %-->ROIs on right
    dist_left_70_flicker_left = TFR.attRight.ev.f2{2};   %-->ROIs on right
    
    %% Choose necessary data
    
    disp('selecting data')
    %find last non-nan datapoint
    t_end = TFR.attLeft.ev.f1{1, 1}.time(find(~isnan(squeeze(TFR.attLeft.ev.f1{1, 1}.powspctrm(1,1,1,:))),1,'last'));
    %is the squeeze correct? on which dimenstion should it not be nan?
    
    cfg = [];
    cfg.latency     = [0 t_end]; %the time of interest
    cfg.avgoverfreq = 'yes';
    cfg.avgoverrpt  = 'yes';
    cfg.avgovertime = 'yes';
    cfg.nanmean     = 'yes';
    cfg.channel     = planars; %only choose planars for ROI selection
    
    %select data separately for 63 and 70
    cfg.frequency = [63 63];   %frequency of config 1
    RFT.c1.indiv.cue_left_63_left    = ft_selectdata(cfg,cue_left_63_flicker_left);
    RFT.c1.indiv.dist_left_63_left   = ft_selectdata(cfg,dist_left_63_flicker_left);
    RFT.c2.indiv.cue_right_63_right  = ft_selectdata(cfg,cue_right_63_flicker_right);
    RFT.c2.indiv.dist_right_63_right = ft_selectdata(cfg,dist_right_63_flicker_right);
    
    cfg.frequency = [70 70];   %frequency of config 2
    RFT.c1.indiv.cue_right_70_right  = ft_selectdata(cfg,cue_right_70_flicker_right);
    RFT.c1.indiv.dist_right_70_right = ft_selectdata(cfg,dist_right_70_flicker_right);
    RFT.c2.indiv.cue_left_70_left    = ft_selectdata(cfg,cue_left_70_flicker_left);
    RFT.c2.indiv.dist_left_70_left   = ft_selectdata(cfg,dist_left_70_flicker_left);
    
    %save subjects data in the main cell
    RFT.c1.group.cue_left_63_left{subj,1}    = RFT.c1.indiv.cue_left_63_left;
    RFT.c1.group.dist_left_63_left{subj,1}   = RFT.c1.indiv.dist_left_63_left;
    RFT.c2.group.cue_right_63_right{subj,1}  = RFT.c2.indiv.cue_right_63_right;
    RFT.c2.group.dist_right_63_right{subj,1} = RFT.c2.indiv.dist_right_63_right;
    
    RFT.c1.group.cue_right_70_right{subj,1}  = RFT.c1.indiv.cue_right_70_right;
    RFT.c1.group.dist_right_70_right{subj,1} = RFT.c1.indiv.dist_right_70_right;
    RFT.c2.group.cue_left_70_left{subj,1}    = RFT.c2.indiv.cue_left_70_left;
    RFT.c2.group.dist_left_70_left{subj,1}   = RFT.c2.indiv.dist_left_70_left;
    
end
% disp('Saving BU')
% save([saveFolder filesep 'RFT/TFR_all_subs_RFT_BU'],'RFT','-v7.3')

%% Average over all subjects and contrast R vs. L

disp('grand averaging over all configs')
cfg = [];
RFT.c1.group.frq_grnd_avg.cue_left_63_left    = ft_freqgrandaverage(cfg,RFT.c1.group.cue_left_63_left{numSub,1});
RFT.c1.group.frq_grnd_avg.dist_left_63_left   = ft_freqgrandaverage(cfg,RFT.c1.group.dist_left_63_left{numSub,1});
RFT.c2.group.frq_grnd_avg.cue_right_63_right  = ft_freqgrandaverage(cfg,RFT.c2.group.cue_right_63_right{numSub,1});
RFT.c2.group.frq_grnd_avg.dist_right_63_right = ft_freqgrandaverage(cfg,RFT.c2.group.dist_right_63_right{numSub,1});

RFT.c1.group.frq_grnd_avg.cue_right_70_right  = ft_freqgrandaverage(cfg,RFT.c1.group.cue_right_70_right{numSub,1});
RFT.c1.group.frq_grnd_avg.dist_right_70_right = ft_freqgrandaverage(cfg,RFT.c1.group.dist_right_70_right{numSub,1});
RFT.c2.group.frq_grnd_avg.cue_left_70_left    = ft_freqgrandaverage(cfg,RFT.c2.group.cue_left_70_left{numSub,1});
RFT.c2.group.frq_grnd_avg.dist_left_70_left   = ft_freqgrandaverage(cfg,RFT.c2.group.dist_left_70_left{numSub,1});

disp('Calculating ROI')
frqCntrst.flick63.left  = squeeze(RFT.c1.group.frq_grnd_avg.cue_left_63_left.powspctrm) ... 
                        - squeeze(RFT.c1.group.frq_grnd_avg.dist_left_63_left.powspctrm);   % sens on right for 63Hz
frqCntrst.flick63.right = squeeze(RFT.c2.group.frq_grnd_avg.cue_right_63_right.powspctrm) ... 
                        - squeeze(RFT.c2.group.frq_grnd_avg.dist_right_63_right.powspctrm); % sens on left for 63Hz

frqCntrst.flick70.right = squeeze(RFT.c1.group.frq_grnd_avg.cue_right_70_right.powspctrm) ... 
                        - squeeze(RFT.c1.group.frq_grnd_avg.dist_right_70_right.powspctrm); % sens on left for 70Hz
frqCntrst.flick70.left  = squeeze(RFT.c2.group.frq_grnd_avg.cue_left_70_left.powspctrm) ... 
                        - squeeze(RFT.c2.group.frq_grnd_avg.dist_left_70_left.powspctrm);   % sens on right for 70Hz

[sens_flick63_l,frqCntrst.flick63.left.ordIdx] = sortrows(frqCntrst.flick63.left,'ascend');   % a vector of frq contrasts ordered from highest to lowest difference when 63 flickers on left
[sens_flick63_r,frqCntrst.flick63.right.ordIdx] = sortrows(frqCntrst.flick63.right,'ascend'); % a vector of frq contrasts ordered from highest to lowest difference when 63 flickers on right

[sens_flick70_r,frqCntrst.flick70.right.ordIdx] = sortrows(frqCntrst.flick70.right,'ascend'); % a vector of frq contrasts ordered from highest to lowest difference when 70 flickers on right
[sens_flick70_l,frqCntrst.flick70.left.ordIdx] = sortrows(frqCntrst.flick70.left,'ascend');   % a vector of frq contrasts ordered from highest to lowest difference when 70 flickers on left

ROI.flick63.left  = RFT.c1.group.frq_grnd_avg.cue_left_63_left.label(frqCntrst.flick63.left.ordIdx(1:ROI_num)); 
ROI.flick63.right = RFT.c2.group.frq_grnd_avg.cue_right_63_right.label(frqCntrst.flick63.right.ordIdx(1:ROI_num));

ROI.flick70.right = RFT.c1.group.frq_grnd_avg.cue_right_70_right.label(frqCntrst.flick70.right.ordIdx(1:ROI_num)); 
ROI.flick70.left  = RFT.c2.group.frq_grnd_avg.cue_left_70_left.label(frqCntrst.flick70.left.ordIdx(1:ROI_num));

disp('saving ROI')
save([saveFolder filesep 'RFT/ROI_RFT/ROI_63'],'ROI.flick63.left','ROI.flick63.right')
save([saveFolder filesep 'RFT/ROI_RFT/ROI_70'],'ROI.flick70.right','ROI.flick70.left')

end