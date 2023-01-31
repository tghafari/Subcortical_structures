% function [] = RFT_ROI_topoplot_indiv(numSubs,badSubs,ROI_num,five_ROI)
% topoplot individual ROIs and save
% numSubs = 35; %number of all subjects
% badSubs = 28; %code of bad subjects
% ROI_num = 5;  %how many ROIs to plot
% five_ROI = 2; %1-> plot all ROI_num ROI, 2-> plot only the best ROI


% Load necessary data
addpath Z:\MATLAB\fieldtrip-20210328 %Windows
saveFolder = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\';
saveFolderFig = 'Z:\Load\Results\FieldTripPlots\';

ft_defaults

load([saveFolder 'RFT' filesep 'cue_dist_diff_RFT_subs_individually'])
load([saveFolder 'RFT' filesep 'comb_planar_labels.mat'])
load([saveFolder 'RFT' filesep 'ROI_RFT/ROI_R_sens_indiv'])
load([saveFolder 'RFT' filesep 'ROI_RFT/ROI_L_sens_indiv'])
load([saveFolder 'RFT' filesep 'ROI_RFT/selected_indiv_ROIs'])



for sub = setxor(1:numSubs,badSubs)
    diff_RFT_left.label = labels;
    diff_RFT_left.avg  = RFT_left_cue_left_dist_left_cntrst(:,sub);
    diff_RFT_left.time = 1;
    diff_RFT_right.label = labels;
    diff_RFT_right.avg  = RFT_right_cue_right_dist_right_cntrst(:,sub);
    diff_RFT_right.time = 1;
    
    % Define topoplot congifs
    cfg = [];
    cfg.layout = 'neuromag306cmb'; %cmb for combined planar gradients
    cfg.zlim = 'maxabs'; % maxabs gives you even distribution of positive and negative scale (best for contrasts)
    cfg.gridscale = 200; %not necessary, but makes it prettier
    cfg.style = 'straight'; %also optional, but this leaves out the contour lines, an aesthetic choice
    cfg.comment = 'no'; %to suppress the text usually plotted with the topo
    cfg.markercolor = [0.3 0.3 0.3]; %Plots the sensors positions in light grey instead of black (looks better with ROI plotted)
    
    figure(1); % where attention is on left (lit-up sensors on right)
    ft_topoplotER(cfg,diff_RFT_left)
    title('attention left, RFT left')
    
    figure(2); % where attention is on right (lit-up sensors on left)
    ft_topoplotER(cfg,diff_RFT_right)
    title('attention right, RFT right')
    
    %% Plot ROIs on top of them
    if five_ROI == 1
        % Now to plot the sensor selection (ROI) on top of this you can use
        %ft_plot_lay with ROI indices
        
        % Prepare layouts for plotting
        cfg = [];
        cfg.layout = 'neuromag306cmb';
        layout = ft_prepare_layout(cfg);
        clmp = colormap('jet');
        
        for roi = 1:ROI_num
            figure(1)
            hold on
            ft_plot_layout(layout,'chanindx',ROI_all_subs_right_sens_idx(roi,sub),'label','no','box','no','pointsymbol','.','pointcolor',clmp(roi*50,:),'pointsize',4) %plot solid inside
            ft_plot_layout(layout,'chanindx',ROI_all_subs_right_sens_idx(roi,sub),'label','no','box','no','pointsymbol','o','pointcolor',clmp(roi*50,:),'pointsize',6) %plots outside ring
            % colorbar
            
            figure(2)
            hold on
            ft_plot_layout(layout,'chanindx',ROI_all_subs_left_sens_idx(roi,sub),'label','no','box','no','pointsymbol','.','pointcolor',clmp(roi*50,:),'pointsize',4)
            ft_plot_layout(layout,'chanindx',ROI_all_subs_left_sens_idx(roi,sub),'label','no','box','no','pointsymbol','o','pointcolor',clmp(roi*50,:),'pointsize',6)
            % colorbar
        end
        saveas(figure(1),[saveFolderFig 'ROI' filesep num2str(sub) '_left_RFT_ROI.jpg']);
        saveas(figure(2),[saveFolderFig 'ROI' filesep num2str(sub) '_right_RFT_ROI.jpg']);
        close all
        
    % Plot only the selected ROIs on top of them
    elseif five_ROI == 2
        
        % Prepare layouts for plotting
        cfg = [];
        cfg.layout = 'neuromag306cmb';
        layout = ft_prepare_layout(cfg);
        
        figure(1)
        hold on
        if chos_roi_R(sub,:)>5; chos_roi_R(sub,:)=chos_roi_R(sub,:)/100; end
        ft_plot_layout(layout,'chanindx',ROI_all_subs_right_sens_idx(chos_roi_R(sub,:),sub),'label','no','box','no','pointsymbol','.','pointcolor',[1 0 0],'pointsize',4) %plot solid inside
        ft_plot_layout(layout,'chanindx',ROI_all_subs_right_sens_idx(chos_roi_R(sub,:),sub),'label','no','box','no','pointsymbol','o','pointcolor',[1 0 0],'pointsize',6) %plots outside ring
        % colorbar
        
        figure(2)
        hold on
        if chos_roi_L(sub,:)>5; chos_roi_L(sub,:)=chos_roi_L(sub,:)/100; end
        ft_plot_layout(layout,'chanindx',ROI_all_subs_left_sens_idx(chos_roi_L(sub,:),sub),'label','no','box','no','pointsymbol','.','pointcolor',[1 0 0],'pointsize',4)
        ft_plot_layout(layout,'chanindx',ROI_all_subs_left_sens_idx(chos_roi_L(sub,:),sub),'label','no','box','no','pointsymbol','o','pointcolor',[1 0 0],'pointsize',6)
        % colorbar
        
        saveas(figure(1),[saveFolderFig 'ROI' filesep num2str(sub) '_slctd_left_RFT_ROI.jpg']);
        saveas(figure(2),[saveFolderFig 'ROI' filesep num2str(sub) '_slctd_right_RFT_ROI.jpg']);
        close all
    end
end
