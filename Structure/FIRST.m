% set FSL environment -- check without the last 4 lines
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE','NIFTI_GZ'); % this to tell what the output type would be
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;


addpath '/Volumes/jenseno-avtemporal-attention/Load/MRI_data';
filePath='/Volumes/jenseno-avtemporal-attention/Load/MRI_data';
 
% filePath='/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/';
% addpath '/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/';

for nSub=16
Subject=['S' num2str(nSub)];
cmd=['/usr/local/fsl/bin/run_first_all -i ' filePath '/Raw_Data/' Subject '_mri.nii -o ' Subject]; %error here
mkdir([filePath '/Processed Data/'],[Subject '-FIRST'])
cd([filePath '/Processed Data/' Subject '-FIRST'])
disp(['segmenting ' Subject])
system(cmd)
fprintf('%s done\n',Subject)
end

% cmdCheck=['/usr/local/fsl/bin/first_roi_slicesdir ' filePath '/' Subject
% '_mri.nii ' Subject '-L_Amyg_first.nii.gz']; I dunno how it works yet.
