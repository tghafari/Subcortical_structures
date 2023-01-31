%% ROI analysis
clear;clc;clf;close all

%% Load clean MEG file
% data_folder='/rds/projects/j/jenseno-avtemporal-attention/Load/MEG Data/proc_data/'; %Portal
% addpath /rds/projects/j/jenseno-avtemporal-attention/MATLAB/fieldtrip-20210328 %Portal
% saveFolder = '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual_Load/FieldTrip/Results/group_level/Alpha';
 
saveFolder = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Alpha'; %Windows
data_folder='Z:\Load\MEG Data\proc_data\'; %Windows
addpath Z:\MATLAB\fieldtrip-20210328 %Windows

ft_defaults

badSubs=[23,28];
numSub = setxor(1:35,badSubs);
% TFR_attRight_alpha_all_subs = cell(length(numSub),1); %for all TFR
% TFR_attLeft_alpha_all_subs  = cell(length(numSub),1);
TFR_attRight_alpha_all_subs_pow = nan(102,length(numSub)); %for power only
TFR_attLeft_alpha_all_subs_pow = nan(102,length(numSub));

for subj=numSub
    if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
    disp(['loading ' sub])
    load([data_folder sub filesep sub '_TFR_LF_dt_correct_only.mat']);fprintf('Done\n')
    
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
    cfg.latency     = [-.85 0]; %the time of interest
    cfg.frequency   = [8 13];   %frequency of interest
    cfg.avgoverfreq = 'yes';
    cfg.avgoverrpt  = 'yes';
    cfg.avgovertime = 'yes';
    cfg.nanmean     = 'yes';
    cfg.channel     = planars; %only choose planars for ROI selection
    %select data
    TFR_attRight_alpha_planars = ft_selectdata(cfg,TFR.attRight.LF); %separately run for each hemisphere's ROI always with R-L
    TFR_attLeft_alpha_planars  = ft_selectdata(cfg,TFR.attLeft.LF);

%     % Choosing all the TFR data to select ROI
%     TFR_attRight_alpha_all_subs{subj,1} = TFR_attRight_alpha_planars;
%     TFR_attLeft_alpha_all_subs{subj,1}  = TFR_attLeft_alpha_planars;
     
    % Only choosing powerspectrum for selecting ROIs
    TFR_attRight_alpha_all_subs_pow(:,subj) = TFR_attRight_alpha_planars.powspctrm;
    TFR_attLeft_alpha_all_subs_pow(:,subj)  = TFR_attLeft_alpha_planars.powspctrm;

end
disp('Saving BU...')
% save([saveFolder filesep 'TFR_dt_all_subs_alpha_pow'],'TFR_attRight_alpha_all_subs','TFR_attLeft_alpha_all_subs')

%% Average over all subjects and contrast R vs. L only using power spectrums
disp('Selecting ROI')
load([saveFolder filesep 'comb_planar_labels_alpha_sym']) %loads labels and a matrix containing indices of sensors in the order of accordance in left and right
% with normalization
num_ROI = 5; %number of sensors
[~,id_ord]=sort(ord_sens); %get the index of each row (use as reference to sort all channels)

RLPowContrast_nl = (TFR_attRight_alpha_all_subs_pow - TFR_attLeft_alpha_all_subs_pow) ./ (TFR_attRight_alpha_all_subs_pow + TFR_attLeft_alpha_all_subs_pow);
RLPowContrast_nl_avg = nanmean(RLPowContrast_nl,2);
%if you don't need symmetric ROIs:
%[BRL,idxRL] = sortrows(RLPowContrast_nl_avg,'ascend'); %for R-L and vice
%versa for L-R
% ROI_attR = alpha_labels(idxRL(1:num_ROI)); %5 sensors will make ROI
% ROI_attL = alpha_labels(idxLR(1:num_ROI));
% inxRL = idxRL(1:num_ROI);
% inxLR = idxLR(1:num_ROI);


% for symmetric ROI
RLPowContrast_nl_avg_ord = RLPowContrast_nl_avg(id_ord,:); %sort powspctrms in a way that first left is in accordance with first right and so on
RLPowContrast_left_nl_avg = RLPowContrast_nl_avg_ord(1:51,:); %first half are left sensors
RLPowContrast_right_nl_avg = RLPowContrast_nl_avg_ord(52:102,:); %second half are corresponding right sensors
RLPowContrast_right_left_sensors_diff = RLPowContrast_right_nl_avg - RLPowContrast_left_nl_avg; %subtract power in left sensors from right sensors to get symmetric ROIs

[BRL,idxRL] = sortrows(RLPowContrast_right_left_sensors_diff,'descend'); 
% a vector of data ordered from highest to lowest negative difference in R-L

%Choose the first 5 idxRL from the ordered list of labels in right and left
% load([saveFolder filesep 'labels_in_accordance']) %these labels are not
% combined (thus the *2)
% ROI_lbl_R = {labels_accord{idxRL(1:5)*2-1,3}};
% ROI_lbl_L = {labels_accord{idxRL(1:5)*2-1,1}};

%easier to do it manually
ROI_lbl_R = {'MEG2312+2313','MEG2322+2323','MEG2032+2033','MEG2432+2433','MEG2442+2443'};
ROI_lbl_L = {'MEG1912+1913','MEG1942+1943','MEG2042+2043','MEG1642+1643','MEG1632+1633'};
idxR = [87,88,77,93,94]; %label numbers in alpha_labels cell
idxL = [71,74,78,62,61];

disp('saving ROI')
save([saveFolder filesep 'ROI_alpha/ROI_dt_right_sym'],'ROI_lbl_R','idxR')
save([saveFolder filesep 'ROI_alpha/ROI_dt_left_sym'],'ROI_lbl_L','idxL')
save([saveFolder filesep 'alpha_avg_contrasts_nl'],'RLPowContrast_nl_avg','LRPowContrast_nl_avg');

%% Average over all subjects and contrast R vs. L using the whole TFR data -- takes much longer
% disp('Calculating ROI')
% cfg = [];
% frq_grnd_avg_attRight = ft_freqgrandaverage(cfg,TFR_attRight_alpha_all_subs{numSub,1});
% frq_grnd_avg_attLeft  = ft_freqgrandaverage(cfg,TFR_attLeft_alpha_all_subs{numSub,1});
% 
% RLPowContrast = squeeze(frq_grnd_avg_attRight.powspctrm)-squeeze(frq_grnd_avg_attLeft.powspctrm)./squeeze(frq_grnd_avg_attRight.powspctrm)+squeeze(frq_grnd_avg_attLeft.powspctrm);
% LRPowContrast = squeeze(frq_grnd_avg_attLeft.powspctrm)-squeeze(frq_grnd_avg_attRight.powspctrm)./squeeze(frq_grnd_avg_attLeft.powspctrm)+squeeze(frq_grnd_avg_attRight.powspctrm);
% 
% [BRL,idxRL] = sortrows(RLPowContrast,'ascend'); % a vector of data ordered from highest to lowest negative difference in R-L
% [BLR,idxLR] = sortrows(LRPowContrast,'ascned'); % a vector of data ordered from highest to lowest negative difference in L-R
% 
% ROI_attR = frq_grnd_avg_attRight.label(idxRL(1:5)); %5 sensors will make ROI
% ROI_attL = frq_grnd_avg_attLeft.label(idxLR(1:5));
% %note that the indices above are basically the 5 first and 5 last indices
% of one matrix as they contain the highest difference in R-L and L-R

% disp('saving ROI')
% save([saveFolder filesep 'ROI_alpha/ROI_dt_attRigth'],'ROI_attR','idxRL')
% save([saveFolder filesep 'ROI_alpha/ROI_dt_attLeft'],'ROI_attL','idxLR')


