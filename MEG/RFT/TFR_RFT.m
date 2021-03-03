function [] = RFT_TFR(sub,av_type)
%RFT TFR analysis

%% Load clean MEG file
proc_folder='/rds/projects/j/jenseno-avtemporal-attention/Load/MEG Data/proc_data/'; %Portal
addpath /rds/projects/j/jenseno-avtemporal-attention/MATLAB/fieldtrip-20200320 %Portal

% proc_folder='Z:\Load\MEG Data\proc_data\'; %Windows
% addpath Z:\MATLAB\fieldtrip-20200320 %Windows

ft_defaults

correct_only = 1;
if nargin<2
choice_made = 0;
fprintf('Select averaging type: \n')
fprintf('[1] RFT-induced\n')
fprintf('[2] RFT-evoked\n')  %better choice
fprintf('[3] HF-cue-locked\n')
fprintf('[4] HF-target-locked\n') %better choice

while ~choice_made
    choice=input('Enter number: ');
    if choice>5 || choice<1
        disp('Invalid choice, please try again');
    else
        av_type = choice;
        switch choice
            case 1
                disp('induced selected')
            case 2
                disp('evoked selected')
            case 3
                disp('cue-locked selected')
            case 4
                disp('target-locked selected')
        end
        choice_made = 1;
    end
end
end
if nargin <1
choice_made = 0;
sub_folders=dir([proc_folder filesep 'S*']);
fprintf('Select subject to analyse: \n')
switch av_type
    case 1
        cnt=1;
        for s=1:size(sub_folders,1)
            if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_all_clean.mat'])>0
                fprintf(['[' int2str(cnt) '] ' sub_folders(s).name '\n'])
                sel(cnt)=s;
                cnt=cnt+1;
            end
        end
    case 2
        cnt=1;
        for s=1:size(sub_folders,1)
            if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_tl.mat'])>0
                fprintf(['[' int2str(cnt) '] ' sub_folders(s).name '\n'])
                sel(cnt)=s;
                cnt=cnt+1;
            end
        end
        
    case 3
        cnt=1;
        for s=1:size(sub_folders,1)
            if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_all_clean.mat'])>0
                fprintf(['[' int2str(cnt) '] ' sub_folders(s).name '\n'])
                sel(cnt)=s;
                cnt=cnt+1;
            end
        end
        
    case 4
        cnt=1;
        for s=1:size(sub_folders,1)
            if exist([proc_folder sub_folders(s).name filesep sub_folders(s).name '_all_clean_dt.mat'])>0
                fprintf(['[' int2str(cnt) '] ' sub_folders(s).name '\n'])
                sel(cnt)=s;
                cnt=cnt+1;
            end
        end
        
end

while ~choice_made
    choice=input('Enter number: ');
    if choice>(cnt-1) || choice<1
        disp('Invalid choice, please try again');
    else
        sub=sub_folders(sel(choice)).name;
        disp(['Subject ' sub ' selected'])
        choice_made=1;
    end
end
end
switch av_type
    case 1
        fprintf(['Loading ' proc_folder sub filesep sub '_all_clean.mat...'])
        load([proc_folder sub filesep sub '_all_clean.mat']);disp('Done')
    case 2
        fprintf(['Loading ' proc_folder sub filesep sub '_tl.mat...'])
        load([proc_folder sub filesep sub '_tl.mat']);disp('Done')
    case 3
        fprintf(['Loading ' proc_folder sub filesep sub '_all_clean.mat...'])
        load([proc_folder sub filesep sub '_all_clean.mat']);disp('Done')
    case 4
        fprintf(['Loading ' proc_folder sub filesep sub '_all_clean_dt.mat...'])
        load([proc_folder sub filesep sub '_all_clean_dt.mat']);disp('Done')
end

%% RFT timecourse

if av_type < 3
    %[1] high frequencies - high frequency resolution
    cfg              = [];
    cfg.foi          = 63;
    cfg.pad          = 4;%'nextpow2';
    cfg.toi          = -2.3:0.01:1;
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.keeptrials   = 'yes';
    cfg.t_ftimwin    = repmat(0.5,1,length(cfg.foi));%10./cfg.foi;%
    
    
    for c=1:2 % 2 RFT configurations
        for l=1:4 % 4 load conditions
            
            if av_type == 1 % Induced
                disp('Computing high frequency TFR for RFT-induced')
                
                cfg.channel      = strmatch('MEG',data{1}.left{1}.meg.label);
                
                %left
                if correct_only
                    cfg.trials = find(data{c}.left{l}.behavior(:,10)>0);
                end
                
                cfg.foi=63;
                disp(['Calculating 63Hz RFT trials for configuration ' int2str(c) ', attention left, load condition ' int2str(l) ' (induced)'])
                TFR.left.RFT.ind.f1{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg);
                TFR.left.RFT.ind.f1{c,l}.grad = data{c}.grad;
                
                disp(['Calculating 70Hz RFT trials for configuration ' int2str(c) ', attention left, load condition ' int2str(l)  ' (induced)'])
                cfg.foi = 70;
                TFR.left.RFT.ind.f2{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg);
                TFR.left.RFT.ind.f2{c,l}.grad = data{c}.grad;
                
                %right
                if correct_only
                    cfg.trials = find(data{c}.right{l}.behavior(:,10)>0);
                end
                
                cfg.foi=63;
                disp(['Calculating 63Hz RFT trials for configuration ' int2str(c) ', attention right, load condition ' int2str(l) ' (induced)'])
                TFR.right.RFT.ind.f1{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg);
                TFR.right.RFT.ind.f1{c,l}.grad = data{c}.grad;
                
                cfg.foi = 70;
                disp(['Calculating 70Hz RFT trials for configuration ' int2str(c) ', attention right, load condition ' int2str(l) ' (induced)'])
                TFR.right.RFT.ind.f2{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg);
                TFR.right.RFT.ind.f2{c,l}.grad = data{c}.grad;
                
            elseif av_type == 2 % Evoked
                disp('Computing high frequency TFR for RFT-evoked')
                
                cfg.channel = strmatch('MEG',tl{1}.left{1}.label);
                
                %all trials for evoked (since there is only one) -> what does
                %this mean?
                cfg.trials='all';
                
                cfg.foi=63;
                disp(['Calculating 63Hz RFT trials for configuration ' int2str(c) ', attention left, load condition ' int2str(l) ' (evoked)'])
                if correct_only
                    TFR.left.RFT.ev.f1{c,l} = ft_freqanalysis(cfg,tl_correct{c}.left{l});
                    TFR.left.RFT.ev.f1{c,l}.grad = tl_correct{c}.left{l}.grad;
                else
                    TFR.left.RFT.ev.f1{c,l} = ft_freqanalysis(cfg,tl{c}.left{l}); %#ok<UNRCH>
                    TFR.left.RFT.ev.f1{c,l}.grad = tl{c}.left{l}.grad;
                end
                
                cfg.foi=70;
                disp(['Calculating 70Hz RFT trials for configuration ' int2str(c) ', attention left, load condition ' int2str(l) ' (evoked)'])
                if correct_only
                    TFR.left.RFT.ev.f2{c,l} = ft_freqanalysis(cfg,tl_correct{c}.left{l});
                    TFR.left.RFT.ev.f2{c,l}.grad = tl_correct{c}.left{l}.grad;
                else
                    TFR.left.RFT.ev.f2{c,l} = ft_freqanalysis(cfg,tl{c}.left{l}); %#ok<UNRCH>
                    TFR.left.RFT.ev.f2{c,l}.grad = tl{c}.left{l}.grad;
                end
                
                cfg.foi=63;
                disp(['Calculating 63Hz RFT trials for configuration ' int2str(c) ', attention right, load condition ' int2str(l) ' (evoked)'])
                if correct_only
                    TFR.right.RFT.ev.f1{c,l} = ft_freqanalysis(cfg,tl_correct{c}.right{l});
                    TFR.right.RFT.ev.f1{c,l}.grad = tl_correct{c}.right{l}.grad;
                else
                    TFR.right.RFT.ev.f1{c,l} = ft_freqanalysis(cfg,tl{c}.right{l}); %#ok<UNRCH>
                    TFR.right.RFT.ev.f1{c,l}.grad = tl{c}.right{l}.grad;
                end
                
                cfg.foi=70;
                disp(['Calculating 70Hz RFT trials for configuration ' int2str(c) ', attention right, load condition ' int2str(l) ' (evoked)'])
                if correct_only
                    TFR.right.RFT.ev.f2{c,l} = ft_freqanalysis(cfg,tl_correct{c}.right{l});
                    TFR.right.RFT.ev.f2{c,l}.grad = tl_correct{c}.right{l}.grad;
                else
                    TFR.right.RFT.ev.f2{c,l} = ft_freqanalysis(cfg,tl{c}.right{l}); %#ok<UNRCH>
                    TFR.right.RFT.ev.f2{c,l}.grad = tl{c}.right{l}.grad;
                end
                
            end
        end
    end
    
    %combine planar data
    disp('Creating combined planar data')
    for c=1:2 %both configurations
        for l=1:4 %4 load conditions
            if av_type == 1
                TFR.left.RFT.ind.f1{c,l}  = ft_combineplanar([],TFR.left.RFT.ind.f1{c,l});
                TFR.right.RFT.ind.f1{c,l} = ft_combineplanar([],TFR.right.RFT.ind.f1{c,l});
                
                TFR.left.RFT.ind.f2{c,l}  = ft_combineplanar([],TFR.left.RFT.ind.f2{c,l});
                TFR.right.RFT.ind.f2{c,l} = ft_combineplanar([],TFR.right.RFT.ind.f2{c,l});
                
            elseif av_type == 2
                TFR.left.RFT.ev.f1{c,l}  = ft_combineplanar([],TFR.left.RFT.ev.f1{c,l});
                TFR.right.RFT.ev.f1{c,l} = ft_combineplanar([],TFR.right.RFT.ev.f1{c,l});
                
                TFR.left.RFT.ev.f2{c,l}  = ft_combineplanar([],TFR.left.RFT.ev.f2{c,l});
                TFR.right.RFT.ev.f2{c,l} = ft_combineplanar([],TFR.right.RFT.ev.f2{c,l});
            end
        end
    end
    
    
    % Append allloads -- configs should remain separated
    for c=1:2
        if av_type == 1
            TFR.attLeft.ind.f1{c}  = ft_appendfreq([],TFR.left.RFT.ind.f1{c,:});
            TFR.attRight.ind.f1{c} = ft_appendfreq([],TFR.right.RFT.ind.f1{c,:});
            TFR.attLeft.ind.f2{c}  = ft_appendfreq([],TFR.left.RFT.ind.f2{c,:});
            TFR.attRight.ind.f2{c} = ft_appendfreq([],TFR.right.RFT.ind.f2{c,:});
        elseif av_type == 2
            TFR.attLeft.ev.f1{c}  = ft_appendfreq([],TFR.left.RFT.ev.f1{c,:});
            TFR.attRight.ev.f1{c} = ft_appendfreq([],TFR.right.RFT.ev.f1{c,:});
            TFR.attLeft.ev.f2{c}  = ft_appendfreq([],TFR.left.RFT.ev.f2{c,:});
            TFR.attRight.ev.f2{c} = ft_appendfreq([],TFR.right.RFT.ev.f2{c,:});
        end
    end
end

%% HF TFR
if av_type > 2
    switch av_type
        case 3 %Cue-locked
            disp('Computing high frequency TFR- cue locked')
            %[1] high frequencies - high frequency resolution
            cfg              = [];
            cfg.foi          = 50:0.5:100;
            cfg.pad          = 4;%'nextpow2';
            cfg.toi          = -2.3:0.01:1;
            cfg.output       = 'pow';
            cfg.method       = 'mtmconvol';
            cfg.taper        = 'hanning';
            cfg.t_ftimwin    = repmat(0.5,1,length(cfg.foi));%10./cfg.foi;%
            cfg.channel      = strmatch('MEG',data{1}.left{1}.meg.label);
            
            for c=1:2 %both configurations
                for l=1:4 %4 load conditions
                    
                    disp(['Calculating TFR for configuration ' int2str(c) ', attention left, load condition ' int2str(l) ' (induced)'])
                    if correct_only
                        cfg.trials = find(data{c}.left{l}.behavior(:,10)>0);
                        beh.left{c,l}=data{c}.left{l}.behavior(cfg.trials,:);
                    else
                        beh.left{c,l}=data{c}.left{l}.behavior;
                    end
                    TFR.left.HF.ind{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg);
                    TFR.left.HF.ind{c,l}.grad = data{c}.grad;
                    
                    if correct_only
                        cfg.trials = find(data{c}.right{l}.behavior(:,10)>0);
                        beh.right{c,l}=data{c}.right{l}.behavior(cfg.trials,:);
                    else
                        beh.right{c,l}=data{c}.right{l}.behavior;
                    end
                    
                    disp(['Calculating TFR for configuration ' int2str(c) ', attention right, load condition ' int2str(l) ' (induced)'])
                    TFR.right.HF.ind{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg);
                    TFR.right.HF.ind{c,l}.grad = data{c}.grad;
                    
                    %all trials for evoked (since there is only one)
                    cfg.trials='all';
                    
                    disp(['Calculating TFR for configuration ' int2str(c) ', attention left, load condition ' int2str(l) ' (evoked)'])
                    if correct_only
                        TFR.left.HF.ev{c,l} = ft_freqanalysis(cfg,tl_correct{c}.left{l});
                    else
                        TFR.left.HF.ev{c,l} = ft_freqanalysis(cfg,tl{c}.left{l});
                    end
                    TFR.left.HF.ev{c,l}.grad = data{c}.grad;
                    
                    disp(['Calculating TFR for configuration ' int2str(c) ', attention right, load condition ' int2str(l) ' (evoked)'])
                    if correct_only
                        TFR.right.HF.ev{c,l} = ft_freqanalysis(cfg,tl_correct{c}.right{l});
                    else
                        TFR.right.HF.ev{c,l} = ft_freqanalysis(cfg,tl{c}.right{l});
                    end
                    TFR.right.HF.ev{c,l}.grad = data{c}.grad;
                end
            end
            
        case 4 %Discrimination target-locked
            disp('Computing high frequency TFR- target locked')
            %[1] high frequencies - high frequency resolution
            cfg              = [];
            cfg.foi          = 50:0.5:100;
            cfg.pad          = 4;%'nextpow2';
            cfg.toi          = -2.3:0.01:1;
            cfg.output       = 'pow';
            cfg.method       = 'mtmconvol';
            cfg.taper        = 'hanning';
            cfg.t_ftimwin    = repmat(0.5,1,length(cfg.foi));
            cfg.channel      = strmatch('MEG',data{1}.left{1}.meg_dt.label);
            
            for c=1:2 %both configurations
                for l=1:4 %4 load conditions
                    
                    disp(['Calculating TFR for configuration ' int2str(c) ', attention left, load condition ' int2str(l) ' (induced)'])
                    if correct_only
                        cfg.trials = find(data{c}.left{l}.behavior(:,10)>0);
                        beh.left{c,l}=data{c}.left{l}.behavior(cfg.trials,:);
                    else
                        beh.left{c,l}=data{c}.left{l}.behavior;
                    end
                    TFR.left.HF.ind{c,l} = ft_freqanalysis(cfg,data{c}.left{l}.meg_dt);
                    TFR.left.HF.ind{c,l}.grad = data{c}.grad;
                    
                    if correct_only
                        cfg.trials = find(data{c}.right{l}.behavior(:,10)>0);
                        beh.right{c,l}=data{c}.right{l}.behavior(cfg.trials,:);
                    else
                        beh.right{c,l}=data{c}.right{l}.behavior;
                    end
                    
                    disp(['Calculating TFR for configuration ' int2str(c) ', attention right, load condition ' int2str(l) ' (induced)'])
                    TFR.right.HF.ind{c,l} = ft_freqanalysis(cfg,data{c}.right{l}.meg_dt);
                    TFR.right.HF.ind{c,l}.grad = data{c}.grad;
                end
            end
    end
    
    %combine planar data
    disp('Creating combined planar data')
    for c=1:2 %both configurations
        for l=1:4 %4 load conditions
            TFR.left.HF.ind{c,l}=ft_combineplanar([],TFR.left.HF.ind{c,l});
            TFR.right.HF.ind{c,l}=ft_combineplanar([],TFR.right.HF.ind{c,l});
            
            if av_type == 3
                TFR.left.HF.ev{c,l}=ft_combineplanar([],TFR.left.HF.ev{c,l});
                TFR.right.HF.ev{c,l}=ft_combineplanar([],TFR.right.HF.ev{c,l});
            end
        end
    end
    
    %append all configurations togeather
        TFR.attLeft.HF.ind  = ft_appendfreq([],TFR.left.HF.ind{:});
        TFR.attRight.HF.ind = ft_appendfreq([],TFR.right.HF.ind{:});
    if av_type ==3
        TFR.attLeft.HF.ev  = ft_appendfreq([],TFR.left.HF.ev{:});
        TFR.attRight.HF.ev = ft_appendfreq([],TFR.right.HF.ev{:});
    end
end
%% Save data
fprintf(['Saving data for ' sub '...\n'])
switch av_type
    case 1
        if correct_only
            save([proc_folder sub filesep sub '_TFR_RFT_correct_only.mat'],'TFR','beh','-v7.3');
        else
            save([proc_folder sub filesep sub '_TFR_RFT.mat'],'TFR','beh','-v7.3');
        end
    case 2
        if correct_only
            save([proc_folder sub filesep sub '_TFR_RFT_tl_correct_only.mat'],'TFR','-v7.3');
        else
            save([proc_folder sub filesep sub '_TFR_RFT_tl.mat'],'TFR','-v7.3');
        end
    case 3
        if correct_only
            save([proc_folder sub filesep sub '_TFR_HF_correct_only.mat'],'TFR','beh','-v7.3');
        else
            save([proc_folder sub filesep sub '_TFR_HF.mat'],'TFR','beh','-v7.3');
        end
    case 4
        if correct_only
            save([proc_folder sub filesep sub '_TFR_HF_correct_only_dt.mat'],'TFR','beh','-v7.3');
        else
            save([proc_folder sub filesep sub '_TFR_HF_dt.mat'],'TFR','beh','-v7.3');
        end
        
end
fprintf('..Done\n')

end