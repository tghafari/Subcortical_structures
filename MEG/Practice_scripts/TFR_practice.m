%% Prepare data--no need, you can simply run Tjerk's B4 and B3 scripts
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
tl_attLeft = ft_timelockanalysis([],attLeft);
tl_attRight = ft_timelockanalysis([],attRight);

%Define frequency analysis cgf-lower freq

cfg=[];
cfg.output='pow';
cfg.channel='MEGGRAD';
cfg.taper='hanning';
cfg.method='mtmconvol'; %'mtmfft';
cfg.foi=1:.5:30;
cfg.t_ftimwin=ones(length(cfg.foi),1)*.5;
cfg.toi=.5:.05:1;
cfg.keeptrials='no';
cfg.pad='nextpow2';

%two orthogonal gradometers
PSD.low_attLeft.ind.grad=ft_freqanalysis(cfg,attLeft);
PSD.low_attRight.ind.grad=ft_freqanalysis(cfg,attRight);
PSD.low_attLeft.tl.grad =ft_freqanalysis(cfg,tl_attLeft);
PSD.low_attRight.tl.grad=ft_freqanalysis(cfg,tl_attRight);

%combine gradometers
cfg_cmb=[];cfg_cmb.method='sum';
PSD.low_attLeft.comb=ft_combineplanar(cfg_cmb,PSD.low_attLeft.ind.grad);
PSD.low_attRight.comb=ft_combineplanar(cfg_cmb,PSD.low_attRight.ind.grad);
PSD.low_attLeft.tl.comb=ft_combineplanar(cfg_cmb,PSD.low_attLeft.tl.grad);
PSD.low_attRight.tl.comb=ft_combineplanar(cfg_cmb,PSD.low_attRight.tl.grad);

% %magnetometers

% cfg.channel=mags;
% PSD.attLeft.ind.mag =ft_freqanalysis(cfg,attLeft);
% PSD.attRight.ind.mag =ft_freqanalysis(cfg,attRight);
% PSD.attLeft.ev.mag =ft_freqanalysis(cfg,tl_attLeft);
% PSD.attRight.ev.mag =ft_freqanalysis(cfg,tl_attRight);