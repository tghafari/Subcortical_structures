%Incomplete
%ft_timelockanalysis (event-related field vs time) of ft_freqanalysis
%(oscillatory activity power/coherence vs frequency)

figure;
cfg = [];
cfg.channel=planars;
% cfg.layout       = 'neuromag306all.lay';
cfg.channel='MEG2032+2033';
cfg.baseline = [-2.3 -1.3]; 
cfg.xlim = [-1.3 -0.3]; 
cfg.ylim = [-3e-12 3e-11]; 
ft_singleplotER(cfg, tl_attLeft);

figure;
cfg = [];
cfg.baseline     = [-2.3 -1.3]; 
cfg.baselinetype = 'relative';
cfg.xlim         = [-2.3 .5];
cfg.zlim         = [0.4 2] ;	        
cfg.layout       = 'neuromag306all.lay';
% 'circular';
ft_multiplotER(cfg,tl_attRight); 
% tl_attLeft);


%ft_freqanalysis
figure;
cfg = [];
cfg.baseline     = [-0.5 -0.3]; 
cfg.baselinetype = 'relative';
cfg.xlim=[-0.5 1];
cfg.zlim         = [0.4 2] ;	        
cfg.layout = 'neuromag306all.lay';
ft_multiplotTFR(cfg, tfr_low_left_c);


