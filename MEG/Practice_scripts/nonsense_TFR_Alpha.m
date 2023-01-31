%frequency analysis and TFRs
clear;clc;
%% Load clean MEG file
% data_folder='/rds/projects/j/jenseno-avtemporal-attention/Load/MEG Data/proc_data/'; %Portal
% addpath /rds/projects/j/jenseno-avtemporal-attention/MATLAB/fieldtrip-20200320 %Portal
% load '/rds/projects/j/jenseno-avtemporal-attention/MATLAB/Perceptual Load/FieldTrip/TjerksCodes/ROI_alpha_4_correct_only'
% ft_defaults

data_folder='Z:\Load\MEG Data\proc_data\'; %Windows
addpath Z:\MATLAB\fieldtrip-20200320 %Windows
load 'Z:\MATLAB\Perceptual Load\FieldTrip\TjerksCodes\ROI_alpha_4_correct_only'
ft_defaults

subj=input('Enter Subject Number:\n');
if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
disp(['loading ' sub])
load([data_folder sub filesep sub '_TFR_LF_correct_only.mat']);fprintf('Done\n')

%% Channel selection
MEG_sens=strmatch('MEG',TFR.left.LF.ind{1}.label);
sens_type=str2num(cellfun(@(x) x(end),TFR.left.LF.ind{1}.label(MEG_sens),'UniformOutput',1)); 

planars=sort([MEG_sens(sens_type==2) ; MEG_sens(sens_type==3)]); 
mags=MEG_sens(sens_type==1);

%Posterior sensor selection - FROM TJERK
load Z:\MATLAB\fieldtrip-20200320\template\layout\neuromag306cmb.lay
Cluster_nrs=[152,16:25,264];

post_roi=[];
post_roi_label='';
for i=1:length(Cluster_nrs)
    inds=strmatch(['MEG' int2str(Cluster_nrs(i))],lay.label);
    post_roi=[post_roi ; inds];
    post_roi_label=[post_roi_label ; lay.label(inds)];
end

%% Append all configs and loads -- run on the output of Tjerk's script
%append all configurations togeather
TFR.attRight.LF = ft_appendfreq([],TFR.right.LF.ind{:}); 
TFR.attLeft.LF  = ft_appendfreq([],TFR.left.LF.ind{:});

%==========================================================================%
c1_left=cell(1,4); c1_right=cell(1,4); c2_left=cell(1,4); c2_right=cell(1,4);
cnt=1;

for l = 1:4
    c1_left{cnt}=TFR.left.LF.ind{1,l};
    c1_right{cnt}=TFR.right.LF.ind{1,l};
    c2_left{cnt}=TFR.left.LF.ind{2,l};
    c2_right{cnt}=TFR.right.LF.ind{2,l};
    cnt=cnt+1;
end

%attention (cue) direction only--ignore frequency tagging and load 
TFR.attLeft.LF=ft_appendfreq([],c1_left{:},c2_left{:});
TFR.attRight.LF=ft_appendfreq([],c1_right{:},c2_right{:});
%==========================================================================%

%% Select ROI- based on highest right-left contrast in alpha, the first 10 channels 

num = 20; %number of channels

cfg=[];
cfg.latency     = [.15 1]; %the time of interest
cfg.frequency   = [8 13]; %frequency of interest
% cfg.avgovertime = 'yes';
cfg.avgoverfreq = 'yes';
cfg.avgoverrpt  = 'yes';
cfg.channel     = planars;% make a pre selection of the channles to pick from (let's say, take the best contrast considering all the occipital and parietal channels, and this is done as in the first lines of the code by using ft_channelselection)

TFR_attRight_LF_timSlctd = ft_selectdata(cfg,TFR.attRight.LF);
TFR_attLeft_LF_timSlctd  = ft_selectdata(cfg,TFR.attLeft.LF); %this is the power of all sensors for alpha frequency averaged over timwin of interest in att left condition

RLContrast = nanmean(TFR_attRight_LF_timSlctd.powspctrm,3)-nanmean(TFR_attLeft_LF_timSlctd.powspctrm,3); %shouldn't be diffrent than normalized over the sum of powspctrm
LRContrast = nanmean(TFR_attLeft_LF_timSlctd.powspctrm,3)-nanmean(TFR_attRight_LF_timSlctd.powspctrm,3);

[BRL,idxRL] = sortrows(RLContrast,'descend'); % a vector of data ordered from highest to lowest difference in R-L
[BLR,idxLR] = sortrows(LRContrast,'descend'); % a vector of data ordered from highest to lowest difference in R-L

ROIRL = TFR_attRight_LF_timSlctd.label(idxRL(1:num));
ROILR = TFR_attLeft_LF_timSlctd.label(idxLR(1:num));

%% plot MI - average over ROI and timwin 

cfg = [];
cfg.latency     = [.15 1]; %the time of interest
cfg.frequency   = [8 13];   %frequency of interest
cfg.avgoverfreq = 'yes';
% cfg.avgoverchan = 'yes';
cfg.avgoverrpt  = 'yes'; 
% cfg.channel     = {ROI_right{:}}; 
% cfg.channel     = {ROI_left{:}};
% [78,71,75,61]; %R [87,77,94,86]; %L is not correct  %run separately than R ROI
% idxR(1:num); %Cecilia's ROI (both hemispheres)

TFR_attRight_LF_chanSlctd = ft_selectdata(cfg,TFR.attRight.LF); %separately run for each hemisphere's ROI always with R-L 
TFR_attLeft_LF_chanSlctd  = ft_selectdata(cfg,TFR.attLeft.LF);

pwrRight = squeeze(TFR_attRight_LF_chanSlctd.powspctrm); %separately run for each hemisphere's ROI always with R-L 
pwrLeft  = squeeze(TFR_attLeft_LF_chanSlctd.powspctrm); %? why is attLeft over right ROI is positive?

MI_overTime_inROI = (pwrRight-pwrLeft)./(pwrRight+pwrLeft); %separately run for each hemisphere's ROI always with R-L 


figure(1); subplot(1,2,1);
hold on;
% plot(TFR_attRight_LF_chanSlctd.time,pwrRight,'b');
% plot(TFR_attLeft_LF_chanSlctd.time, pwrLeft,'r');
plot(TFR_attRight_LF_chanSlctd.time, MI_overTime_inROI, 'k', 'LineWidth', 2);
% legend('att_R', 'att_L', 'MI','Location','northwest');
xlabel('Time'); ylabel('Modulation Index'); title ('MI over L-ROIs in alpha-Sub13')
xlim([-.3 .5]); 
line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;

%plot the single modulation index
cfg.avgovertime = 'yes';
TFR_attRight_LF = ft_selectdata(cfg,TFR.attRight.LF);
TFR_attLeft_LF  = ft_selectdata(cfg,TFR.attLeft.LF);

MI_L = (TFR_attRight_LF.powspctrm-TFR_attLeft_LF.powspctrm)./(TFR_attRight_LF.powspctrm+TFR_attLeft_LF.powspctrm); %Run after R_ROI
MI_R = (TFR_attRight_LF.powspctrm-TFR_attLeft_LF.powspctrm)./(TFR_attRight_LF.powspctrm+TFR_attLeft_LF.powspctrm); %Run after L_ROI

figure(); 
bar([1,2],[MI_L,MI_R])
set(gca,'XTicklabel',{'L__ROI','R__ROI'},'XTickLabelRotation',45); ylabel('MI(R-L)');

figure()
ALI = MI_R-MI_L;
bar(ALI)
ylabel('ALI(attR-attL)')

%% Calculate ALI
cfg=[];
cfg.latency     = [-.3 .5]; %the time of interest
cfg.frequency   = [8 13]; %frequency of interest
cfg.avgoverchan = 'yes';
cfg.avgoverfreq = 'yes';

% cfg.channel     = channels_Left;% make a pre selection of the channles to pick from (let's say, take the best contrast considering all the occipital and parietal channels, and this is done as in the first lines of the code by using ft_channelselection)
TFR.attLeft.LF.chanSlctd  = ft_selectdata(cfg,TFR.attLeft.LF);

% cfg.channel     = channels_Right; %same as before but right hemisphere
TFR.attRight.LF.chanSlctd = ft_selectdata(cfg,TFR.attRight.LF);

powerLeft   = squeeze(TFR.attLeft.LF.chanSlctd.powspctrm);   %power in left hemisphere at different timepoints
powerRight  = squeeze(TFR.attRight.LF.chanSlctd.powspctrm); 
powerDiffRL = (powerRight-powerLeft);%./(powerRight+powerLeft); %this is how ALI is computed basically

figure(1); clf; 
hold on;
plot(TFR.attLeft.LF.chanSlctd.time, mean(powerLeft), 'r');
plot(TFR.attRight.LF.chanSlctd.time, mean(powerRight), 'b');
plot(TFR.attRight.LF.chanSlctd.time, mean(powerDiffRL), 'k', 'LineWidth', 2);
legend('L sens', 'R sens', 'Difference','Location','northwest');
xlabel('Time'); ylabel('Lateralization Index'); title ('Alpha Lateralization Index selected sens')
xlim([-.3 .5]); line(ylim,[0,0],'Color', 'b','LineWidth',1.5); box on;


%% Plot TFRs 

cfg=[];
cfg.baseline=[-2 -.4];
cfg.baselinetype='relative';
cfg.xlim=[-2 1];
cfg.zlim=[.4 2]; %for plotTFR
% cfg.ylim=[-3e-12 2e-12]; %for plotER
% cfg.channel='MEG2032+2033'; %for singleplot
% cfg.showlabels='yes'; %for singleplot
cfg.colorbar='yes';
cfg.layout='neuromag306cmb.lay'; %for multiplot

figure(1);ft_singleplotTFR(cfg,TFR.attLeft.LF); %low frequency -TFR
title('attention left')
figure(2);ft_singleplotTFR(cfg,TFR.attRight.LF); %low frequency -TFR
title('attention right')

figure(1); ft_multiplotTFR(cfg,TFR.attLeft.LF);
title('attention left')
figure(2); ft_multiplotTFR(cfg,TFR.attRight.LF);
title('attention right')

cfg_topo=[];
cfg_topo.baseline=[-.5 -.1];
cfg_topo.baselinetype='relative';
cfg_topo.xlim=[.5 1];
cfg_topo.ylim=[-.1 1.2];
cfg_topo.zlim=[.9 1.4];
cfg_topo.showlabels='no'; 
cfg_topo.layout='neuromag306cmb.lay'; 
figure(3);ft_topoplotTFR(cfg_topo,TFR.attLeft.LF); %low frequency -TFR
figure(4);ft_topoplotTFR(cfg_topo,TFR.attRight.LF); %low frequency -TFR


% figure(1);ft_multiplotTFR(cfg,PSD.low_attLeft.comb); %low frequency -PSD
% figure(2);ft_multiplotTFR(cfg,PSD.low_attRight.comb); %low frequency -PSD
% figure(3);ft_multiplotTFR(cfg,PSD.low_attLeft.tl.comb); %low frequency tl -PSD
% figure(4);ft_multiplotTFR(cfg,PSD.low_attRight.tl.comb); %low frequency tl -PSD
% 
% 


