%analyse the volumes
clear;clc

% addpath '/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/';
% filePath='/Users/Tara/Documents/Projects/CHBHVisit/PerceptualLoad/Data/MRI_Data/Processed_Data/';

% filePath='/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data/'; %Portal
% addpath '/Volumes/rdsprojects/j/jenseno-avtemporal-attention/Load/MRI_data/Processed Data/'; %Portal

% addpath '/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed_Data/';
% filePath='/Volumes/jenseno-avtemporal-attention/Load/MRI_data/Processed_Data/';

addpath 'Z:\Projects\Load\MRI_data\Processed_Data'
filePath   = 'Z:\Projects\Load\MRI_data\Processed_Data';
% saveFolder = 'Z:\Projects\Load\Results\FieldTrip Plots\matFiles';


lables=[10,11,12,13,16,17,18,26,49,50,51,52,53,54,58];
structures={'L-Thal','L-Caud','L-Puta','L-Pall','BrStem /4th Ventricle',...
    'L-Hipp','L-Amyg','L-Accu','R-Thal','R-Caud','R-Puta',...
    'R-Pall','R-Hipp','R-Amyg','R-Accu'};
HemisComp={'_{Thalamus}','_{CN}','_{Putamen}','_{GP}','_{Hippocampus}','_{Amygdala}','_{n. Accumbens}'};
structures={'Thalamus', 'Caudate n.', 'Putamen', 'Globus Pallidus', 'Hippocampus', 'Amygdala', 'n. Accumbens'};

load([filePath filesep 'AllVolumes'])
hemisVol=volumes(:,setxor(1:length(lables),5));
%% Calculate LV index
LV=zeros(35,7);
for sub=1:35
    for i=1:7
        % [h(i),p(i)]=ranksum(hemisVol(:,i),hemisVol(:,i+7)); %statistical test
        LV(sub,i)=(hemisVol(sub,i+7)-hemisVol(sub,i))/(hemisVol(sub,i)+hemisVol(sub,i+7)); %R-L/R+L as in Cecilia's paper
    end
end

% save(['Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices' filesep 'LV_all'],'LV')
%% Throw out the outliers + calculate p-value
h=zeros(1,7);p=zeros(1,7);
format short
for ii=1:7
%     LV(abs(LV(:,ii))<=(nanmean(LV(:,ii))-(2.5*(nanstd(LV(:,ii))))),ii)=NaN;
    [p(ii),h(ii)]=signrank(LV(:,ii));
end

% calculate mean and std for all structures
mean_LV = nanmean(LV);
std_LV = nanstd(LV);
%% histogram

% colormap={'gold', 'blue violet', 'midnight blue', 'dark red', 'olive drab', 'indigo', 'light blue'}; %names according to https://www.rapidtables.com/web/color/RGB_Color.html
% colormap=[255,215,0;138,43,226;25,25,112;139,0,0;107,142,35;75,0,130;173,216,230];  % rgd codes
colormap={'#FFD700', '#8A2BE2', '#191970', '#8B0000', '#6B8E23', '#4B0082', '#ADD8E6'};  % Hexcode
figure;
for his=1:7
    ax = subplot(2,5,his);
    hold on
    h = histfit(LV(:,his),5,'normal');
    h(1).FaceColor = colormap{his};
    h(2).Color = 'black';
    h(2).LineWidth = 0.25;
    title(ax, structures{his},'FontSize', 10)
    ylabel('Nr. of Subjects'); xlabel(['LV' (HemisComp{his})]);
    ylim([1 25]);
    line([0, 0], ylim, 'LineWidth', .25, 'Color', 'k', 'LineStyle', ':')
    txt = ['\itp\rm\bf', sprintf(' = %.3f', p(his))];
    t = text(min(xlim)+.2*max(xlim), max(ylim)-3, txt);
    t.FontSize = 8;
    t.FontWeight = 'bold';
end

