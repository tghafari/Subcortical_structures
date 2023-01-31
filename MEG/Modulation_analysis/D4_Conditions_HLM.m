function [] = D5_Conditions_HLM(badSubs)

data_folder = 'Z:\Load\MEG Data\proc_data\'; %Windows
save_folder = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices';
load 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Alpha\ROI_alpha\ROI_dt_attRight_nl'
load 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Alpha\ROI_alpha\ROI_dt_attLeft_nl'
cd(save_folder)

% badSubs = [23,28];

%% Calculate behavioral index
for subj = setxor(1:35,badSubs)
    if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
    disp(['loading ' sub])
    load([data_folder sub filesep sub '_TFR_LF_dt_correct_only.mat']);fprintf('Done\n')
    %% Channel selection
    MEG_sens = strmatch('MEG',TFR.left.LF.ind{1}.label);
    sens_type = str2num(cellfun(@(x) x(end),TFR_trials.left.LF{1}.label(MEG_sens),'UniformOutput',1));
    
    planars=sort([MEG_sens(sens_type==2) ; MEG_sens(sens_type==3)]);
    mags=MEG_sens(sens_type==1);
    
    %% Calculate and plot MI over time
    cfg = [];
    cfg.latency     = [-.85 0]; %the time of interest
    cfg.frequency   = [8 13];   %frequency of interest
    cfg.avgoverfreq = 'yes';
    cfg.avgoverchan = 'yes';
    cfg.avgoverrpt  = 'yes';
    
    for cnf=1:2
        for lds=1:4
            %plot for right ROI -- ROI attLeft is on right sensors
            cfg.channel     = {ROI_attL{:}};
            TFR_attRight_LF_chanSlctd_R{cnf,ld} = ft_selectdata(cfg,TFR_trials.right.LF{cnf,ld}); %separately run for each hemisphere's ROI always with R-L
            TFR_attLeft_LF_chanSlctd_R{cnf,ld}  = ft_selectdata(cfg,TFR_trials.left.LF{cnf,ld});
            
            pwrRight{cnf,ld} = squeeze(TFR_attRight_LF_chanSlctd_R{cnf,ld}.powspctrm); %separately run for each hemisphere's ROI always with R-L
            pwrLeft{cnf,ld}  = squeeze(TFR_attLeft_LF_chanSlctd_R{cnf,ld}.powspctrm); %? why is attLeft over right ROI is positive?
            
            MI_overTime_inROI_R{cnf,ld}(:,subj) = (pwrRight-pwrLeft{cnf,ld})./(pwrRight+pwrLeft{cnf,ld}); %separately run for each hemisphere's ROI always with R-L
            
            %plot for left ROI -- ROI attRight is on left sensors
            cfg.channel     = {ROI_attR{:}};
            TFR_attRight_LF_chanSlctd_L{cnf,ld} = ft_selectdata(cfg,TFR_trials.right.LF{cnf,ld}); %separately run for each hemisphere's ROI always with R-L
            TFR_attLeft_LF_chanSlctd_L{cnf,ld}  = ft_selectdata(cfg,TFR_trials.left.LF{cnf,ld});
            
            pwrRight{cnf,ld} = squeeze(TFR_attRight_LF_chanSlctd_L{cnf,ld}.powspctrm);
            pwrLeft{cnf,ld}  = squeeze(TFR_attLeft_LF_chanSlctd_L{cnf,ld}.powspctrm);
            
            MI_overTime_inROI_L{cnf,ld}(:,subj) = (pwrRight-pwrLeft{cnf,ld})./(pwrRight+pwrLeft{cnf,ld});
            
            %% Calculate and bar plot MI
            cfg.avgovertime = 'yes';
            
            cfg.channel     = {ROI_attL{:}};
            TFR_attRight_LF_R{cnf,ld} = ft_selectdata(cfg,TFR_trials.right.LF{cnf,ld});
            TFR_attLeft_LF_R{cnf,ld}  = ft_selectdata(cfg,TFR_trials.left.LF{cnf,ld});
            MI_R{cnf,ld}(:,subj) = (nanmean(TFR_attRight_LF_R{cnf,ld}.powspctrm)-nanmean(TFR_attLeft_LF_R{cnf,ld}.powspctrm))...
                ./(nanmean(TFR_attRight_LF_R{cnf,ld}.powspctrm)+nanmean(TFR_attLeft_LF_R{cnf,ld}.powspctrm)); %Run after L_ROI
            
            cfg.channel     = {ROI_attR{:}};
            TFR_attRight_LF_L{cnf,ld} = ft_selectdata(cfg,TFR_trials.right.LF{cnf,ld});
            TFR_attLeft_LF_L{cnf,ld}  = ft_selectdata(cfg,TFR_trials.left.LF{cnf,ld});
            MI_L{cnf,ld}(:,subj) = (nanmean(TFR_attRight_LF_L{cnf,ld}.powspctrm)-nanmean(TFR_attLeft_LF_L{cnf,ld}.powspctrm))...
                ./(nanmean(TFR_attRight_LF_L{cnf,ld}.powspctrm)+nanmean(TFR_attLeft_LF_L{cnf,ld}.powspctrm)); %Run after R_ROI
        end
    end
    
    %% Save necessary data
    clear TFR_trials TFR
    %     mkdir([saveFolderMat filesep 'Sub' num2str(subj)])
    save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'MI_Data_dt'],'MI_L','MI_R','MI_overTime_inROI_L','MI_overTime_inROI_R');
    save([saveFolderMat filesep 'Sub' num2str(subj) filesep 'All_Data_dt'],'TFR_attRight_LF_chanSlctd_R','TFR_attRight_LF_chanSlctd_L',...
        'TFR_attLeft_LF_chanSlctd_R','TFR_attLeft_LF_chanSlctd_L','TFR_attRight_LF_R','TFR_attRight_LF_L',...
        'TFR_attLeft_LF_R','TFR_attLeft_LF_L','MI_L','MI_R','MI_overTime_inROI_L','MI_overTime_inROI_R');
    saveas(figure(1),[saveFolderFig filesep num2str(subj) '_MIoverTime_dt.jpg']);
    saveas(figure(2),[saveFolderFig filesep num2str(subj) '_MIndex_dt.jpg']);
    
end
end
