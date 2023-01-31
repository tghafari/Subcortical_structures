
for ii = 1:102
    alpha_labels{ii,2} = sum(strncmp(alpha_labels{ii,1},right_sensors_planar,7))+1; %left = 1, righ = 2;
    if alpha_labels{ii,2}==1; alpha_labels{ii,3} = 'left';
    elseif alpha_labels{ii,2}==2; alpha_labels{ii,3} = 'right';
    end
end

% %this is just a reference for ordering sensors according to the right and
% %lef corresponding locations:
% 
% load([saveFolder filesep 'comb_planar_labels_alpha']) %loads labels and a matrix containing indices of sensors in the order of accordance in left and right
% % with normalization
% [~,id_ord]=sort(ord_sens); %get the index of each row (use as reference to sort all channels)
% sensors_in_ord = alpha_labels{:,1}(id_ord,:); %sort powspctrms in a way that first left is in accordance with first right and so on
