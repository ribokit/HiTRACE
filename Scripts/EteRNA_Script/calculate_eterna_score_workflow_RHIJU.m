clear;

load Workspace/R51_workspace.mat  

which_sets = 1:length(sequence);

%% generating area_pred
for j = which_sets  
  seqpos = length(sequence{j})-20 - [1:(length(sequence{j})-20)] + 1;
  [ marks{j}, all_area_pred{j}, mutpos{j} ] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence{j}, 0 , seqpos, data_types );
end

% area_pred for switched structure
if(~isempty(alt_structure))
    alt_lane = [2 4];
    for j = which_sets  
      seqpos = length(sequence{j})-20 - [1:(length(sequence{j})-20)] + 1;
      [ alt_marks{j}, alt_area_pred{j}, alt_mutpos{j} ] = get_predicted_marks_SHAPE_DMS_CMCT( alt_structure, sequence{j}, 0 , seqpos, data_types );
      all_area_pred{j}(:,alt_lane) = alt_area_pred{j}(:,alt_lane);
    end
end

% peak fitting
%overmodlength( sequence{j} )
for j = which_sets
    [area_peak{j}, prof_fit{j}] = do_the_fit_fast( d_bsub{j}, xsel{j}', 0.0, 0);
end

backgd_sub_col = find(strcmp(data_types, 'nomod'));

% This is new -- fixing the modification rate so that we don't add noise. 
% Note that this parameter was poorly constrained anyway.
fixed_overmod_correct = 1.0 * (length( sequence{j} ) - 20)/80;  % this is approximate -- we shoot for single hit over 100 residues.
for j = which_sets
  [ area_bsub{j}, darea_bsub{j}] = overmod_and_background_correct_logL( area_peak{j}, backgd_sub_col, [4:size(area_peak{j},1)-4], all_area_pred{j}, [], fixed_overmod_correct);
end

% These parameters were the 'END' and 'START' of previous scripts. How far to 'cut in' from the end in scoring.
inset_from_5prime = 4;
inset_from_3prime = 4;
ignore_points = [ 10:15  28:32]; % These are the nucleotides that bind FMN directly -- they were constant sequences in the puzzles

% This script is in charge of (1) normalizing the data, (2) making calls on which nucleotide 'switch' appropriately, (3) giving back 
%  graphical display of this scoring, (4) compiling the overall 'switch score'. 
[ switch_score, data_to_output, data_to_output_err ] = calc_switch_score_RHIJU( inset_from_5prime, inset_from_3prime, ignore_points, sequence, seqpos, area_bsub, darea_bsub, all_area_pred, design_names );
fprintf( 'Hit return to continue...\n');
pause;
print('R51_switch_score.png', '-dpng', '-r300');

% This script figures out the EteRNA score, but uses fixed thresholds. So mix/max/threshold are returned as 0, 1, and 0.5.
[ETERNA_score, min_SHAPE, max_SHAPE, threshold_SHAPE] = calc_eterna_score_RHIJU( inset_from_5prime, inset_from_3prime, data_types, data_to_output, sequence, seqpos, area_bsub, all_area_pred, design_names );
print('R51_eterna_score.png', '-dpng', '-r300');

%[name path] = uiputfile('Output.rdat', 'Save to RDAT file');
%outfile = strcat(path,name);

nres = length( data_to_output{1} );
goodbins = [ (nres-inset_from_5prime) : -1 : inset_from_3prime];
outfile = 'test.rdat'

% This script has been modified so that it doesn't output min/max/threshold if  min is zero -- simpler output of RDAT file.
% eterna_create_rdat_files_GUI( outfile, target_names{1}, structure, sequence, ...
% 			      data_to_output, ids, target_names, subrounds, sequence, ...
% 			      design_names, seqpos, goodbins, data_types, [], [], ...
% 			      min_SHAPE, max_SHAPE, threshold_SHAPE, ETERNA_score, ...
% 			      switch_score,[], data_to_output_err,trace_data );