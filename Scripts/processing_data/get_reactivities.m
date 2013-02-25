function [ normalized_reactivity,area_peak_corrected,attenuation_corrected,reactionProb] = get_reactivities( saturated_array,diluted_array,sd_cutoff, bkg_col, ref_peak)
% Fully automated reactivity referencing workflow, starting from
% measurements of area_peak.  See "unsaturate" for more details on
% saturated_array, unsaturated_array, and sd_cutoff.
% 
% bkg_col contains references to the nomod lanes.  Input reactionProbs of lanes
% corresponding to a single nomod.  Nomod is the first lane of the reactionProb.
% e.g. for nomods in lanes 1, 5, and 10, bkg_col = [1:4;5:9;10:12]  Enter
% the lane numbers as they appear in saturated array; even if it is lane 6
% in your quick_look output, enter it as 1 if it is first in
% saturated_array.

% ref_peak is the nucleotide to which you want to normalize reactivites.
% keep in mind that area_peak is oriented 3' to 5'.
% 

% NORMALIZED REACTIVITY VALUES ARE RETURNED 5' to 3'.

% Thomas Mann, November 2012.

area_peak_corrected = [];
attenuation_corrected = [];
reactionProb = [];
normalized_reactivity = [];
reactionProb = {};

area_peak_corrected = unsaturate(saturated_array,diluted_array,sd_cutoff);

[~, num_cols] = size(area_peak_corrected);

num_nomod = length(bkg_col); %number of different background conditions to subtract.

%contains all the steps to process data after peak alignment, assignment,
%and dilution scaling.

for i = 1:num_cols;
    attenuation_corrected(:,i) = attenuation_corrector_v2(area_peak_corrected(:,i))
end;

[~,num_probs] = size(bkg_col);

for i = 1:num_nomod;
    for j = 1:num_probs;
    reactionProb{i} = attenuation_corrected(:,bkg_col(i,j)) - attenuation_corrected(:,bkg_col(i,1));
    end;
end;

for i = 1:length(reactionProb);
    reactivity(:,i) = reactionProb{i}(:) / reactionProb{i}(ref_peak);
end;

normalized_reactivity = []
[~,react_cols] = size(reactivity)
for i = 1:react_cols;
    normalized_reactivity(:,i) = transpose(sequence_reversed(reactivity(:,i)));
end;

end

