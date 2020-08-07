%analyse the volumes
clear;clc

addpath '/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/';
filePath='/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/Processed_Data/';

% addpath '/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed_Data/';
% filePath='/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed_Data/';

% filePath='/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data/';
% addpath '/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data/';

lables=[10,11,12,13,16,17,18,26,49,50,51,52,53,54,58];
structures={'L-Thal','L-Caud','L-Puta','L-Pall','BrStem /4th Ventricle',...
    'L-Hipp','L-Amyg','L-Accu','R-Thal','R-Caud','R-Puta',...
    'R-Pall','R-Hipp','R-Amyg','R-Accu'};
HemisComp={'Thal','Caud','Puta','Pall','Hipp','Amyg','Accu'};

load([filePath 'AllVolumes'])
hemisVol=volumes(:,setxor(1:length(lables),5));
%% Calculate LV index
LV=zeros(35,7);
for sub=1:35
    for i=1:7
        % [h(i),p(i)]=ranksum(hemisVol(:,i),hemisVol(:,i+7)); %statistical test
        
        LV(sub,i)=(hemisVol(sub,i+7)-hemisVol(sub,i))/(hemisVol(sub,i)+hemisVol(sub,i+7));
    end
end

%% Throw out the outliers
h=zeros(1,7);p=zeros(1,7);
format short
for ii=1:7
    LV(abs(LV(:,ii))>=(nanmean(LV(:,ii))+(2.5*(nanstd(LV(:,ii))))),ii)=NaN;
    LV(abs(LV(:,ii))<=(nanmean(LV(:,ii))-(2.5*(nanstd(LV(:,ii))))),ii)=NaN;
    [p(ii),h(ii)]=signrank(LV(:,ii));
end

%% histogram
% LV(isnan(LV))=[]; %omit the nans--disorganizes the rows
colormap={'r','g','b','c','m','y','k'};
figure;
for his=1:7
    subplot(2,4,his)
    hold on
    histfit(LV(:,his),5,'normal')
    histogram(LV(:,his),5,'FaceColor',colormap{his})
    txt=sprintf('p= %.2f',p(his));
    ylabel('Nr. of Subjects'); xlabel(['LV-' (HemisComp{his})]);
    ylim([1 25]);
    text(-.015,18,txt)
end

