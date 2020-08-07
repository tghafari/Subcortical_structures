%Put all the volumes in one matrix
clc;clear

addpath '/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/';
filePath='/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/';

% addpath '/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed Data/';
% filePath='/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed Data/';

% filePath='/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data/';
% addpath '/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data/';

lables=[10,11,12,13,16,17,18,26,49,50,51,52,53,54,58];
structures={'L-Thal','L-Caud','L-Puta','L-Pall','BrStem /4th Ventricle',...
    'L-Hipp','L-Amyg','L-Accu','R-Thal','R-Caud','R-Puta',...
    'R-Pall','R-Hipp','R-Amyg','R-Accu'};
allSubs=35;
rows=cell(1,allSubs);
volumes=zeros(allSubs,numel(lables));

for nSub=1:allSubs
    if numel(num2str(nSub))==1; Subject=['S0' num2str(nSub)]; else; Subject=['S' num2str(nSub)]; end %Create subject name
    rows{nSub}=Subject;
    if isfolder([filePath 'Processed_Data/' Subject '.SubVol']) %check if subject's MRI has been volumed
        cd([filePath 'Processed_Data/' Subject '.SubVol'])
        for low=lables-.5
            if isfile(['volume' num2str(low+.5) '.mat']) %check if subject's structures have been volumed
                load(['volume' num2str(low+.5) '.mat'])
                vol=str2num(VoxVol); %#ok<ST2NM>
            else
                vol=[NaN, NaN];
                warning(['no ' Subject ' structure ' num2str(low+.5) ' found, continuing to the next structure'])
            end
            volumes(nSub,lables==low+.5)=vol(2);
        end
    else
        warning(['no' Subject ' volume found, continuing to the next sub'])
    end
end

volumeTable=array2table(volumes,'RowNames',rows,'VariableNames',structures);
save([filePath 'Processed_Data/AllVolumes'],'volumes','volumeTable')