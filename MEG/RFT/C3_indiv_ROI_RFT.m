function [] = C3_indiv_ROI_RFT(ROI_num)
% find the best ROI for each subject

saveFolder = 'Z:\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\';

load([saveFolder 'RFT' filesep 'TFR_all_subs_nl_pow_RFT.mat'])
load([saveFolder 'RFT' filesep 'normalized_cue_dist_diff_RFT.mat'])
load([saveFolder 'RFT' filesep 'comb_planar_labels.mat'])

RFT_right_cue_right_dist_right_cntrst = cue_right_f1f2_right_nl_power ...
                                      - dist_right_f1f2_right_nl_power;
RFT_left_cue_left_dist_left_cntrst = cue_left_f1f2_left_nl_power ...
                                   - dist_left_f1f2_left_nl_power;
                               
[~,sorted_RFT_right_idx] = sort(RFT_right_cue_right_dist_right_cntrst,'descend');
[~,sorted_RFT_left_idx] = sort(RFT_left_cue_left_dist_left_cntrst,'descend');

ROI_all_subs_left_sens_lbl = labels(sorted_RFT_right_idx(1:ROI_num,:));
ROI_all_subs_left_sens_idx = sorted_RFT_right_idx(1:ROI_num,:);
ROI_all_subs_right_sens_lbl = labels(sorted_RFT_left_idx(1:ROI_num,:));
ROI_all_subs_right_sens_idx = sorted_RFT_left_idx(1:ROI_num,:);

save([saveFolder 'RFT' filesep 'ROI_RFT/ROI_L_sens_indiv'],'ROI_all_subs_left_sens_lbl','ROI_all_subs_left_sens_idx')
save([saveFolder 'RFT' filesep 'ROI_RFT/ROI_R_sens_indiv'],'ROI_all_subs_right_sens_lbl','ROI_all_subs_right_sens_idx')
save([saveFolder 'RFT' filesep 'cue_dist_diff_RFT_subs_individually'],'RFT_right_cue_right_dist_right_cntrst','RFT_left_cue_left_dist_left_cntrst')

end