dirlist  = dir('Z:\MATLAB\fieldtrip-20200320\template\layout\*');
filename = {dirlist(strncmp('neuromag306',{dirlist.name},11)).name};
for i=1:length(filename)
cfg = [];
cfg.layout = filename{i};
layout = ft_prepare_layout(cfg);

figure
ft_plot_layout(layout);
h = title(filename{i});
set(h, 'Interpreter', 'none');

[p, f, x] = fileparts(filename{i});
% print([lower(f) '.png'], '-dpng');
end