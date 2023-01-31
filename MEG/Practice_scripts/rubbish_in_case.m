%put too much effort into these, I'm saving them in case they were needed.
%But basically left RFt on left sensors (and vice verca) doesn't make sense
%as the other hemifield is being flickered with another frequency and isn't
%a good control.

%% calculate MI-RFT left (MI_63_right + MI_70_right) on left sensors
   %%%%% on right sensors %%%%%
    cfg.channel   = {ROI_labels_right_sens{:}}; 
       
    %%%%% right 63 right att- right dist %%%%%
    cfg.frequency = [63 63];   
    
    cue_right_63_right_overTime_R  = ft_selectdata(cfg,cue_right_63_flicker_right);
    dist_right_63_right_overTime_R = ft_selectdata(cfg,dist_right_63_flicker_right);
 
    cue_right_63_right_overTime_pow_R  = squeeze(cue_right_63_right_overTime_R.powspctrm); %separately run for each hemisphere's ROI always with R-L
    dist_right_63_right_overTime_pow_R = squeeze(dist_right_63_right_overTime_R.powspctrm); 
    
    MI_overTime_63_right_right_sens = (cue_right_63_right_overTime_pow_R - dist_right_63_right_overTime_pow_R)...
                                    ./(cue_right_63_right_overTime_pow_R + dist_right_63_right_overTime_pow_R); 
                              
    %%%%% right 70 right att- right dist %%%%%                  
    cfg.frequency = [70 70];  
   
    cue_right_70_right_overTime_R  = ft_selectdata(cfg,cue_right_70_flicker_right);               
    dist_right_70_right_overTime_R = ft_selectdata(cfg,dist_right_70_flicker_right);
    
    cue_right_70_right_overTime_pow_R  = squeeze(cue_right_70_right_overTime_R.powspctrm);          
    dist_right_70_right_overTime_pow_R = squeeze(dist_right_70_right_overTime_R.powspctrm);
    
    MI_overTime_70_right_right_sens = (cue_right_70_right_overTime_pow_R - dist_right_70_right_overTime_pow_R)...
                                    ./(cue_right_70_right_overTime_pow_R + dist_right_70_right_overTime_pow_R); 
    
%% calculate MI-RFT right (MI_63_right + MI_70_right) on right sensors
    MI_overTime_RFT_right_right_sens = MI_overTime_63_right_right_sens + MI_overTime_70_right_right_sens;
    
    %plot MI over time right RFT on right sensors
    figure(1); subplot(1,2,2);
    hold on;
    plot(dist_right_70_right_overTime_R.time, MI_overTime_RFT_right_right_sens, 'k', 'LineWidth', 2);
    
    xlabel('Time'); ylabel('Modulation Index'); title (['MI over R-ROIs in right RFT - Sub' num2str(subj)])
    xlim(cfg.latency);
    line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;

    %%%%% on left sensors %%%%%
    cfg.channel   = {ROI_labels_left_sens{:}}; 

    %%%%% left 63 left att- left dist %%%%%
    cfg.frequency = [63 63];   
    
    cue_left_63_left_overTime_L  = ft_selectdata(cfg,cue_left_63_flicker_left);
    dist_left_63_left_overTime_L = ft_selectdata(cfg,dist_left_63_flicker_left);
 
    cue_left_63_left_overTime_pow_L  = squeeze(cue_left_63_left_overTime_L.powspctrm); %separately run for each hemisphere's ROI always with R-L
    dist_left_63_left_overTime_pow_L = squeeze(dist_left_63_left_overTime_L.powspctrm); 
    
    MI_overTime_63_left_left_sens = (cue_left_63_left_overTime_pow_L - dist_left_63_left_overTime_pow_L)...
                                  ./(cue_left_63_left_overTime_pow_L + dist_left_63_left_overTime_pow_L); 
                              
    %%%%% left 70 left att- left dist %%%%%                  
    cfg.frequency = [70 70];  
   
    cue_left_70_left_overTime_L  = ft_selectdata(cfg,cue_left_70_flicker_left);               
    dist_left_70_left_overTime_L = ft_selectdata(cfg,dist_left_70_flicker_left);
    
    cue_left_70_left_overTime_pow_L  = squeeze(cue_left_70_left_overTime_L.powspctrm);          
    dist_left_70_left_overTime_pow_L = squeeze(dist_left_70_left_overTime_L.powspctrm);
    
    MI_overTime_70_left_left_sens = (cue_left_70_left_overTime_pow_L - dist_left_70_left_overTime_pow_L)...
                                   ./(cue_left_70_left_overTime_pow_L + dist_left_70_left_overTime_pow_L); 
    
    %calculate MI-RFT left (MI_63_right + MI_70_right) on left sensors
    MI_overTime_RFT_left_left_sens = MI_overTime_63_left_left_sens + MI_overTime_70_left_left_sens;
    
    %plot MI over time left RFT on left sensors
    figure(2); subplot(1,2,1);
    hold on;
    plot(dist_left_70_left_overTime_L.time, MI_overTime_RFT_left_left_sens, 'k', 'LineWidth', 2);
    
    xlabel('Time'); ylabel('Modulation Index'); title (['MI over L-ROIs in left RFT - Sub' num2str(subj)])
    xlim(cfg.latency);
    line(xlim,[0,0],'Color', 'b','LineWidth',1.5); box on;
