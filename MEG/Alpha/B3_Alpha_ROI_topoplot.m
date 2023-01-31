%% All necessary topoplots
% Load necessary data
addpath Z:\MATLAB\fieldtrip-20210328 %Windows
saveFolder = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\';

ft_defaults

%% Alpha

load([saveFolder 'Alpha' filesep 'TFR_dt_all_subs_alpha_pow'])
load([saveFolder 'Alpha' filesep 'comb_planar_labels_alpha'])
load([saveFolder 'Alpha' filesep 'alpha_avg_contrasts_nl'])
load([saveFolder 'Alpha' filesep 'ROI_alpha' filesep 'ROI_dt_right_sym'])
load([saveFolder 'Alpha' filesep 'ROI_alpha' filesep 'ROI_dt_left_sym'])

%% Topoplot difference (cue-distr) average power spectrum of all subjects

RL_alpha_pow_contrast.avg = RLPowContrast_nl_avg; 
RL_alpha_pow_contrast.label = alpha_labels;
RL_alpha_pow_contrast.time = 1;

LR_alpha_pow_contrast.avg = LRPowContrast_nl_avg; 
LR_alpha_pow_contrast.label = alpha_labels;
LR_alpha_pow_contrast.time = 1;


% Define topoplot congifs
cfg = [];
cfg.layout = 'neuromag306cmb'; %cmb for combined planar gradients
cfg.zlim = 'maxabs'; % maxabs gives you even distribution of positive and negative scale (best for contrasts)
cfg.gridscale = 200; %not necessary, but makes it prettier
cfg.style = 'straight'; %also optional, but this leaves out the contour lines, an aesthetic choice
cfg.comment = 'no'; %to suppress the text usually plotted with the topo
cfg.markercolor = [0.3 0.3 0.3]; %Plots the sensors positions in light grey instead of black (looks better with ROI plotted)

% Topoplot the contrast 
figure(1);
ft_topoplotER(cfg,RL_alpha_pow_contrast)
title('att right - att left')
colorbar

figure(2);
ft_topoplotER(cfg,LR_alpha_pow_contrast)
title('att left - att right')
colorbar

%% Plot ROIs on top of them

% Prepare layouts for plotting
cfg = [];
cfg.layout = 'neuromag306cmb';
layout = ft_prepare_layout(cfg);

figure(1)
hold on
ft_plot_layout(layout,'chanindx',idxR(1:5),'label','no','box','no','pointsymbol','.','pointcolor',[1 1 1],'pointsize',4) %plot solid inside
ft_plot_layout(layout,'chanindx',idxR(1:5),'label','no','box','no','pointsymbol','o','pointcolor',[1 1 1],'pointsize',6) %plots outside ring

figure(1)
hold on
ft_plot_layout(layout,'chanindx',idxL(1:5),'label','no','box','no','pointsymbol','.','pointcolor',[1 1 1],'pointsize',4)
ft_plot_layout(layout,'chanindx',idxL(1:5),'label','no','box','no','pointsymbol','o','pointcolor',[1 1 1],'pointsize',6)


%% TFR plots
%the script is correct it just doens't show anything yet for RFT

cfg = [];
cfg.baseline   = 'no'; 
cfg.xlim       = [-0.85 -.15];
cfg.zlim       = [-3 3] ;	
cfg.ylim       = [50 90];
cfg.showlabels = 'no';	
cfg.channel    = ROI_labels_left_sens{2};
ft_singleplotTFR(cfg,TFR.attRight.ind.f1{2});


