%checks the orientation of the MRI image

clear;clc

% set FSL environment -- check without the last 4 lines
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE','NIFTI_GZ'); % this to tell what the output type would be
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

addpath '/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/';
filePath='/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/Processed_Data';
cd(filePath)

% addpath '/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed Data';
% filePath='/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed Data';

% filePath='/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data';
% addpath '/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data';

allSubs=35;
%% Calculating volumes
for nSub=1:allSubs
   
    if numel(num2str(nSub))==1; Subject=['S0' num2str(nSub)]; else; Subject=['S' num2str(nSub)]; end %Create subject name
    if isfolder([Subject '.anat']) %check if subject's MRI has been segmented
        
                
            cmd=['/usr/local/fsl/bin/fslorient ' Subject '.anat/T1.nii.gz'];
            [nothing,VoxVol]=system(cmd);
            
            disp(['saving subject' Subject 'structure:' structures(lables==low+.5)])            
            save([filePath filesep 'Processed_Data' filesep Subject '.SubVol/volume'  num2str(low+.5)],'VoxVol')
            disp(['structure saved:' structures(lables==low+.5)])
            
        fprintf('%s done\n',Subject)
    else
        warning(['no' Subject '_all_fast_firstseg found, continuing to the next sub'])
    end
end
