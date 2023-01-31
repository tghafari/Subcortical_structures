clc;clear

% set FSL environment -- check without the last 4 lines
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE','NIFTI_GZ'); % this to tell what the output type would be
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

addpath '/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/';
filePath='/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/';

% addpath '/Volumes/jenseno-avtemporal-attention/Load/MRI_data/';
% filePath='/Volumes/jenseno-avtemporal-attention/Load/MRI_data/';
 
% filePath='/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/';
% addpath '/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/';

% %Donders MRI to check
% addpath '/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/Donders_Data/';
% filePath='/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/Donders_Data/';

%% structural analysis

for nSub=5
if numel(num2str(nSub))==1; Subject=['S0' num2str(nSub)]; else; Subject=['S' num2str(nSub)]; end %Create subject name-Tjerk's
% if numel(num2str(nSub))==1; Subject=['subject0' num2str(nSub)]; else; Subject=['subject' num2str(nSub)]; end %Create subject name-Cecilia's
cmd=['/usr/local/fsl/bin/fsl_anat -i ' filePath 'Raw_Data/' Subject '_mri.nii -o ' filePath 'Processed_Data/' Subject]; %Tjerk's
% cmd=['/usr/local/fsl/bin/fsl_anat -i ' filePath 'Raw_Data/' Subject '/defaced_mri.nii.gz -o ' filePath 'Processed_Data/' Subject]; %Cecilia's
disp(['segmenting ' Subject])
system(cmd)
fprintf('%s done\n',Subject)
end

% cmdCheck=['/usr/local/fsl/bin/first_roi_slicesdir ' filePath '/' Subject
% '_mri.nii ' Subject '-L_Amyg_first.nii.gz']; I dunno how it works yet.
