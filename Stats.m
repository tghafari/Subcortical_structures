clear;clc

% set FSL environment -- check without the last 4 lines
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE','NIFTI_GZ'); % this to tell what the output type would be
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

addpath '/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/';
filePath='/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/';


% addpath '/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed Data';
% filePath='/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed Data';

% filePath='/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data';
% addpath '/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data';

lables=[10,11,12,13,16,17,18,26,49,50,51,52,53,54,58];
structures={'L-Thal','L-Caud','L-Puta','L-Pall','BrStem /4th Ventricle',...
    'L-Hipp','L-Amyg','L-Accu','R-Thal','R-Caud','R-Puta',...
    'R-Pall','R-Hipp','R-Amyg','R-Accu'};

%% Calculating volumes
for nSub=23
   
    if numel(num2str(nSub))==1; Subject=['S0' num2str(nSub)]; else; Subject=['S' num2str(nSub)]; end %Create subject name
    cd([filePath 'Processed_Data/' Subject '.anat/first_results'])
    if isfile('T1_first_all_fast_firstseg.nii.gz') %check if subject's MRI has been segmented
        
    mkdir([filePath filesep 'Processed_Data' filesep],[Subject '.SubVol'])
%   cd([filePath Subject '.SubVol'])
    
        for low=lables-.5
            
            cmd=['/usr/local/fsl/bin/fslstats T1_first_all_fast_firstseg -l ' num2str(low) ' -u ' num2str(low+1) ' -V'];
            disp(['Volumetring:' Subject 'structure:' structures(lables==low+.5)])
            [nothing,VoxVol]=system(cmd);
            
            disp(['saving subject' Subject 'structure:' structures(lables==low+.5)])            
            save([filePath filesep 'Processed_Data' filesep Subject '.SubVol/volume'  num2str(low+.5)],'VoxVol')
            disp(['structure saved:' structures(lables==low+.5)])
            
        end
        fprintf('%s done\n',Subject)
    else
        warning(['no' Subject '_all_fast_firstseg found, continuing to the next sub'])
    end
end
