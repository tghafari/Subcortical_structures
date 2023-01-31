function [] = C4_select_indiv_RFT_ROI()

% Topoplot ROIs for individual subs and select the best roi

roi_topoplot_path = 'Z:\Load\Results\FieldTripPlots\ROI\';
save_path = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\RFT\ROI_RFT\';

chos_roi_L = nan(35,1);
chos_roi_R = nan(35,1);

%plot ranks of colors
clmp = colormap('jet');
figure;for roi=1:5; text(.2*roi,.1*roi,num2str(roi),'FontSize',20,'Color',clmp(roi*50,:)); end
xlim([0 1.5]); ylim([0 0.75]);

badSubs = 28;
for subj = setxor(1:35,badSubs)
    figure;imshow([roi_topoplot_path num2str(subj) '_left_RFT_ROI.jpg'])
    chos_roi_L(subj,:) = input('Enter which roi?');
    figure;imshow([roi_topoplot_path num2str(subj) '_right_RFT_ROI.jpg'])
    chos_roi_R(subj,:) = input('Enter which roi?');
    close; close
end
save([save_path 'selected_indiv_ROIs'],'chos_roi_L','chos_roi_R');

end