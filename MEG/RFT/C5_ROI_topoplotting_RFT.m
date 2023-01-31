%% All necessary topoplots
% Load necessary data
addpath Z:\MATLAB\fieldtrip-20210328 %Windows
saveFolder = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\';

ft_defaults

%% RFT

% Load
load([saveFolder 'RFT' filesep 'AVG_all_subs_nl_pow_RFT.mat'])
load([saveFolder 'RFT' filesep 'comb_planar_labels.mat'])
load([saveFolder 'RFT' filesep 'normalized_cue_dist_diff_RFT.mat'])
load([saveFolder 'RFT' filesep 'ROI_RFT/ROI_R_sens'])
load([saveFolder 'RFT' filesep 'ROI_RFT/ROI_L_sens'])

%% Topoplot average power spectrum of all subejects

% Create the following struct with averaged powerspectrums and labels
    % note that I use ft_topoplotER here, not TFR (because there are no
    % frequencies here)

cue_left_f1f2_left.label = labels;
cue_left_f1f2_left.avg  = cue_left_f1f2_left_avg;
cue_left_f1f2_left.time = 1;
dist_left_f1f2_left.label = labels;
dist_left_f1f2_left.avg  = dist_left_f1f2_left_avg;
dist_left_f1f2_left.time = 1;

cue_right_f1f2_right.label = labels;
cue_right_f1f2_right.avg  = cue_right_f1f2_right_avg;
cue_right_f1f2_right.time = 1;
dist_right_f1f2_right.label = labels;
dist_right_f1f2_right.avg  = dist_right_f1f2_right_avg;
dist_right_f1f2_right.time = 1;

% Define topoplot congifs
cfg = [];
cfg.layout = 'neuromag306cmb'; %cmb for combined planar gradients
cfg.zlim = 'maxabs'; %scale; maxabs gives you even distribution of positive and negative scale (best for contrasts)
cfg.gridscale = 200; %not necessary, but makes it prettier
cfg.style = 'straight'; %also optional, but this leaves out the contour lines, an aesthetic choice
cfg.comment = 'no'; %to suppress the text usually plotted with the topo
cfg.markercolor = [0.3 0.3 0.3]; %Plots the sensors positions in light grey instead of black (looks better with ROI plotted)

figure(1);
ft_topoplotER(cfg,cue_left_f1f2_left)
title('cue left RFT left')
colorbar

figure(2);
ft_topoplotER(cfg,dist_left_f1f2_left)
title('dist left RFT left')
colorbar

figure(3);
ft_topoplotER(cfg,cue_right_f1f2_right)
title('cue right RFT right')
colorbar

figure(4);
ft_topoplotER(cfg,dist_right_f1f2_right)
title('dist right RFT right')
colorbar

%% Topoplot difference (cue-distr) average power spectrum of all subjects

diff_f1f2_left.label = labels;
diff_f1f2_left.avg  = RFT_c1c2_group_powCntrst_left;
diff_f1f2_left.time = 1;
diff_f1f2_right.label = labels;
diff_f1f2_right.avg  = RFT_c1c2_group_powCntrst_right;
diff_f1f2_right.time = 1;

% Define topoplot congifs
cfg = [];
cfg.layout = 'neuromag306cmb'; %cmb for combined planar gradients
cfg.zlim = [0.2 2]; % maxabs gives you even distribution of positive and negative scale (best for contrasts)
cfg.gridscale = 200; %not necessary, but makes it prettier
cfg.style = 'straight'; %also optional, but this leaves out the contour lines, an aesthetic choice
cfg.comment = 'no'; %to suppress the text usually plotted with the topo
cfg.markercolor = [0.3 0.3 0.3]; %Plots the sensors positions in light grey instead of black (looks better with ROI plotted)

figure(5); % where attention is on left (lit-up sensors on right)
ft_topoplotER(cfg,diff_f1f2_left)
title('difference in left')
colorbar

figure(6); % where attention is on right (lit-up sensors on left)
ft_topoplotER(cfg,diff_f1f2_right)
title('difference in right')
colorbar

%% Plot ROIs on top of them
    % Now to plot the sensor selection (ROI) on top of this you can use
    %ft_plot_lay with ROI indices
    
% Prepare layouts for plotting
cfg = [];
cfg.layout = 'neuromag306cmb';
layout = ft_prepare_layout(cfg);

figure(5)
hold on
ft_plot_layout(layout,'chanindx',ROI_inds_right_sens,'label','no','box','no','pointsymbol','.','pointcolor',[0 0 0],'pointsize',4) %plot solid inside
ft_plot_layout(layout,'chanindx',ROI_inds_right_sens,'label','no','box','no','pointsymbol','o','pointcolor',[0 0 0],'pointsize',6) %plots outside ring

figure(6)
hold on
ft_plot_layout(layout,'chanindx',ROI_inds_left_sens,'label','no','box','no','pointsymbol','.','pointcolor',[0 0 0],'pointsize',4)
ft_plot_layout(layout,'chanindx',ROI_inds_left_sens,'label','no','box','no','pointsymbol','o','pointcolor',[0 0 0],'pointsize',6)
