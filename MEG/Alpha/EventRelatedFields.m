%calculates event related fields
%1.select trials of specific condition
%% Load clean MEG file
%folders
data_folder='/rds/projects/j/jenseno-avtemporal-attention/Load/proc_data/';
%/rds/projects/2017/jenseno-01/Tjerk/Load2/proc_data/';

%set path
addpath /rds/projects/j/jenseno-avtemporal-attention/MATLAB/fieldtrip-20200320
ft_defaults
subj=input('Enter Subject Number:\n');
if numel(num2str(subj))==1; sub=['S0' num2str(subj)]; else; sub=['S' num2str(subj)]; end
disp(['loading ' sub])
load([data_folder sub filesep sub '_all_clean.mat']);fprintf('Done\n')

%% Create timelocks per condition
disp('Creating timelocks')
% for c=1:2 %two rapid frequency configurations
%     for l=1:4 %4 load conditions
%         tl{c}.left{l}=ft_timelockanalysis([],data{c}.left{l}.meg);
%         tl{c}.right{l}=ft_timelockanalysis([],data{c}.right{l}.meg);
%         
%         %now with only correct trials        
%         cfg.trials = find(data{c}.left{l}.behavior(:,10)>0);                
%         tl_correct{c}.left{l}=ft_timelockanalysis(cfg,data{c}.left{l}.meg);
%         cfg.trials = find(data{c}.right{l}.behavior(:,10)>0);                
%         tl_correct{c}.right{l}=ft_timelockanalysis(cfg,data{c}.right{l}.meg);        
%     end
% end
% 
% %ignore load condition and attention direction
% cnt=1;
% for l=1:4
%     config{1}{cnt}=ft_appenddata([],data{1}.left{l}.meg,data{1}.right{l}.meg);
%     config{2}{cnt}=ft_appenddata([],data{2}.left{l}.meg,data{2}.right{l}.meg);
%     cnt=cnt+1;
% end

% %append
% c1=ft_appenddata([],config{1}{:});
% c2=ft_appenddata([],config{2}{:});
% all=ft_appenddata([],c1,c2);


cnt=1;
for l=1:4
    c1_left{cnt}=data{1}.left{l}.meg;
    c1_right{cnt}=data{1}.right{l}.meg;
    c2_left{cnt}=data{2}.left{l}.meg;
    c2_right{cnt}=data{2}.right{l}.meg;
    cnt=cnt+1;
end

%attention (cue) direction only--ignore frequency tagging
attLeft=ft_appenddata([],c1_left{:},c2_left{:});
attRight=ft_appenddata([],c1_right{:},c2_right{:});

cfg_sel=[];
cfg_sel.latency=[0.5 1];
attLeft=ft_selectdata(cfg_sel,attLeft); 
attRight=ft_selectdata(cfg_sel,attRight); 
tl_attLeft=ft_timelockanalysis([],attLeft);
tl_attRight=ft_timelockanalysis([],attRight);

%% Combine gradiometers and plot 

cfg=[];
cfg.method='sum';
cmb_tl_attLeft=ft_combineplanar(cfg,tl_attLeft);
cmb_tl_attRight=ft_combineplanar(cfg,tl_attRight);

%Plot ERFs
figure(1);
cfg=[];
cfg.layout='neuromag306cmb.lay';
cfg.baseline=[-.5 -.1]; %check the bseline with manuscript
cfg.xlim=[0.5 1];
cfg.ylim=[-3e-12 2e-12];
ft_multiplotER(cfg,cmb_tl_attLeft,cmb_tl_attRight); 

%Plot Topoplots
% cfg=[];
% cfg.layout='neuromag306cmb.lay';
% cfg.baseline=[-.5 -.1]; %check the bseline with manuscript
% cfg.xlim=[.5 .8];
cfg.zlim=[-3e-12 2e-12];
figure(2);ft_topoplotER(cfg,cmb_tl_attLeft);
figure(3);ft_topoplotER(cfg,cmb_tl_attRight); 




