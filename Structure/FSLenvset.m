setenv( 'FSLDIR', '/usr/local/fsl');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

% set FSL environment
setenv('FSLOUTPUTTYPE','NIFTI_GZ'); % this to tell what the output type would be
% That?s it and you can launch any FSL command as follows,
cmd = '/usr/local/fsl/bin/feat design.fsf';
system(cmd);