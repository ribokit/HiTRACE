%% STEP 1
% make sure you are at the right path, i.e. the following folders are available
% use command pwd to check where you are
filenames = { ...
    'pfl_1D_QC&CM_2015-05-02_1802', ...
    'pfl_1D_QC&CM_2015-05-02_1803', ...
    'pfl_1D_QC&CM_2015-05-02_1804', ...
    'pfl_1D_QC&CM_2015-05-02_1805', ...
};
ylimit = [2501, 6000];
% reorder is for reorganizing data lanes into MATLAB
reorder = [1:8, 17:64, 9, 13, 10, 14, 11, 15, 12, 16];
[d_align, d_ref_align, ylimit, labels] = quick_look(filenames, ylimit, reorder);

%% STEP 2
% the numbers in align_blocks correspond to columns in d_align
% it has nothing to do with reorder
% only align groups of similar nature, e.g. SHAPE
% ddNTPs should never be included here
align_blocks = {1:8, 9:24, 25:40, 41:56, 57:58, 59:60, 61:62, 63:64};
d_align_before_more_alignment = d_align;
d_align_dp_fine = align_by_DP_fine(d_align_before_more_alignment, align_blocks);
% check for quality before accept
d_align = d_align_dp_fine;


%% STEP 3
sequence = 'GGAGACCTCGAGTAGAGGTCAAAAGGGTCGTGACTGGCGAACAGGTGGGAAACCACCGGGGAGCGACCCCGGCATCGATAGCCGCCCGCCTGGGCAAACAACTCGAGTAGAGTTGACAACAAAGAAACAACAACAACAAC';
sequence = strrep(sequence, 'T', 'U');
% a crystallographic structure is recommended
% untested predictions should not be used; use empty string instead
structure = '...((((((.....))))))....(((((((....[[[[....(((((....))))).....)))))))...........(((..]]]]...)))...((((((.....)))))).........................';
% data_types should reflect your experiment setup
data_types = [ ...
    repmat({'nomod'}, 1, 8), repmat({'DMS'}, 1, 16), repmat({'CMCT'}, 1, 16), repmat({'SHAPE'}, 1, 16), ...
    repmat({'ddATP'}, 1, 2), repmat({'ddTTP'}, 1, 2), repmat({'ddCTP'}, 1, 2), repmat({'ddGTP'}, 1, 2), ...
];
% offset = - length_of_5_flanking_region + N_first_nt_ROI_should_be - 1
offset = -24;
% this formula for first_RT_nucleotide is typically not changed
first_RT_nucleotide = length(sequence) - 20 + offset;


%% STEP 4
% only run this line once, it clears xsel
xsel = []; clf;
% this line is safe to run multiple times, for updating or checkign your assignment
[xsel, seqpos, area_pred] = annotate_sequence(d_align, xsel, sequence, offset, data_types, first_RT_nucleotide, structure, [], 2);

%% STEP 5
[area_peak, darea_peak] = fit_to_gaussians(d_align, xsel);


%% STEP 6
ref_segment = 'GAGUA';
% be careful if your ROI has a GAGUA somewhere
ref_peak = get_ref_peak(sequence, ref_segment, offset);
sd_cutoff = 1.5;
% pull out undiluted lanes from area_peak
saturated_idx = [9, 10, 13, 14, 17 ,18 ,21, 22, 25, 26, 29, 30, 33, 34, 37, 38, 41, 42, 45, 46, 49, 50, 53, 54];

% the mean() is because of the example has 2 nomod replicates for each condition (-/+ ligand)
% simple cases should use saturated_array = area_peak(:, saturated_idx);
saturated_array = [ ...
    mean(area_peak(:, 1:2), 2), ...
    mean(area_peak(:, 5:6), 2),...
    area_peak(:, saturated_idx), ...
];
% in the example, indices for diluted lanes are derived from undiluted ones
diluted_array = [ ...
    mean(area_peak(:, 3:4), 2), ...
    mean(area_peak(:, 7:8), 2),...
    area_peak(:, saturated_idx + 2), ...
];

saturated_error = [ ...
    mean(darea_peak(:, 1:2), 2), ...
    mean(darea_peak(:, 5:6), 2),...
    darea_peak(:, saturated_idx), ...
];
diluted_error = [ ...
    mean(darea_peak(:, 3:4), 2), ...
    mean(darea_peak(:, 7:8), 2),...
    darea_peak(:, saturated_idx + 2), ...
];

% the background lane (nomod) for each lane in saturated_array
% in this example, 1 and 2 correspond to (-/+ ligand) nomods
% the numbering is the same as saturated_array, not area_peak
bkg_col  = [1, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2];
% make sure you have the data_types correctly
% in this example, [1, 2, saturated_idx] prepends to 'nomod' strings to
% reflect the (-/+ ligand) backgrounds averaged separately
% simple cases should use data_types(saturated_idx) instead
[normalized_reactivity, normalized_error, seqpos_out] = get_reactivities( ...
    saturated_array, diluted_array, saturated_error, diluted_error, ...
    bkg_col, ref_peak, seqpos, [], data_types([1, 2, saturated_idx]), sequence, offset, sd_cutoff);

%% STEP 7
% indices of normalized_reactivity should be the same as saturated_array
[d_DMS_minus, da_DMS_minus, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 3:4), normalized_error(:, 3:4), [], seqpos_out, sequence, offset); 
[d_DMS_plus, da_DMS_plus, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 7:8), normalized_error(:, 7:8), [], seqpos_out, sequence, offset); 

[d_CMCT_minus, da_CMCT_minus, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 11:12), normalized_error(:, 11:12), [], seqpos_out, sequence, offset); 
[d_CMCT_plus, da_CMCT_plus, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 15:16), normalized_error(:, 15:16), [], seqpos_out, sequence, offset); 

[d_SHAPE_minus, da_SHAPE_minus, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 19:20), normalized_error(:, 19:20), [], seqpos_out, sequence, offset); 
[d_SHAPE_plus, da_SHAPE_plus, flags] = average_data_filter_outliers( ...
    normalized_reactivity(:, 23:24), normalized_error(:, 23:24), [], seqpos_out, sequence, offset); 


%% STEP 8
labels = data_types;
for i = 1:length(labels) - 8;
    if i <= 8;
        if mod(i - 1, 8) < 4;
            labels{i} = [labels{i}, ' (-)'];
        else
            labels{i} = [labels{i}, ' (+)'];
        end;
    else
        if mod(i - 9, 16) < 8;
            labels{i} = [labels{i}, ' (-)'];
        else
            labels{i} = [labels{i}, ' (+)'];
        end;
    end;
end;
% you can leave out the 2nd and 3rd line optionals
% append a 'X' to sequence allow printing the final additional band location
print_xsel_split(d_align, xsel, seqpos_out, ['X', sequence], offset, area_pred, labels, ...
    [3 0 0 50 75 100 4 8 52 0 8 56], [0 1 0 1 1], ...
    {'RNA Puzzle #13: ZMP Aptamer', '', '', 'T47', 'May 2015'}, [25 15 9 8 11 1 1 2]);


% this is a complicated case for plotting
% simple cases only need commands of plot and errorbar
% commands of make_lines, legend, title, set are decorative
figure();
set_print_page(gcf, 0);
% do not directly copy these commands for your own dataset
subplot(3, 1, 1);
plot(seqpos_out, d_DMS_minus, 'b', 'linewidth', 2); hold on;
plot(seqpos_out, d_DMS_plus, 'r', 'linewidth', 2); hold on;
errorbar(seqpos_out, d_DMS_minus, da_DMS_minus, 'b'); hold on;
errorbar(seqpos_out, d_DMS_plus, da_DMS_plus, 'r'); hold on;
make_lines_horizontal(-0.5, 'k');
make_lines([min(seqpos_out):0, 72:(max(seqpos_out) + 20)] - 0.5, [0.4 0.4 0.4], 0.5, 1, 0);
make_lines(ref_peak - 0.5, 'y', 1.5, 1, 0);
axis([-25 100 -0.5 3.5]);
legend('(-) ZMP', '(+) 10 uM ZMP');
title('1D DMS', 'fontweight', 'bold', 'fontsize', 15);
set(gca, 'xgrid', 'off', 'ygrid', 'on');
set(gca, 'xtick', [seqpos_out, first_RT_nucleotide + 1:20], 'xticklabel', sequence', 'fontsize', 8);
make_lines([-20:20:100] - 0.5, 'k', 0.5, 0, 1);

subplot(3, 1, 2);
plot(seqpos_out, d_CMCT_minus, 'b', 'linewidth', 2); hold on;
plot(seqpos_out, d_CMCT_plus, 'r', 'linewidth', 2); hold on;
errorbar(seqpos_out, d_CMCT_minus, da_CMCT_minus, 'b'); hold on;
errorbar(seqpos_out, d_CMCT_plus, da_CMCT_plus, 'r'); hold on;
make_lines_horizontal(-0.5, 'k');
make_lines([min(seqpos_out):0, 72:(max(seqpos_out) + 20)] - 0.5, [0.4 0.4 0.4], 0.5, 1, 0);
make_lines(ref_peak - 0.5, 'y', 1.5, 1, 0);
axis([-25 100 -0.5 5]);
legend('(-) ZMP', '(+) 10 uM ZMP');
title('1D CMCT', 'fontweight', 'bold', 'fontsize', 15);
set(gca, 'xgrid', 'off', 'ygrid', 'on');
set(gca, 'xtick', [seqpos_out, first_RT_nucleotide + 1:20], 'xticklabel', sequence', 'fontsize', 8);
make_lines([-20:20:100] - 0.5, 'k', 0.5, 0, 1);

subplot(3, 1, 3);
plot(seqpos_out, d_SHAPE_minus, 'b', 'linewidth', 2); hold on;
plot(seqpos_out, d_SHAPE_plus, 'r', 'linewidth', 2); hold on;
errorbar(seqpos_out, d_SHAPE_minus, da_SHAPE_minus, 'b'); hold on;
errorbar(seqpos_out, d_SHAPE_plus, da_SHAPE_plus, 'r'); hold on;
make_lines_horizontal(-0.5, 'k');
make_lines([min(seqpos_out):0, 72:(max(seqpos_out) + 20)] - 0.5, [0.4 0.4 0.4], 0.5, 1, 0);
make_lines(ref_peak - 0.5, 'y', 1.5, 1, 0);
axis([-25 100 -0.5 3]);
legend('(-) ZMP', '(+) 10 uM ZMP');
title('1D SHAPE', 'fontweight', 'bold', 'fontsize', 15);
set(gca, 'xgrid', 'off', 'ygrid', 'on');
set(gca, 'xtick', [seqpos_out, first_RT_nucleotide + 1:20], 'xticklabel', sequence', 'fontsize', 8);
make_lines([-20:20:100] - 0.5, 'k', 0.5, 0, 1);

print_save_figure(gcf, '1D_plot_compare', './');


%% STEP 9
filename = 'pfl_1D.rdat';
name = 'RNA Puzzle 13';
comments = { ...
    'Preliminary data, averaged over 2 replicates', ...
    'RNA was heat to 90 C for 2 min and cool on ice for 2 min, then incubated at 37 C for 20 min at pH 8.0 with 10 mM MgCl2 to aid folding.', ...
    'Additional sequences at 5'' and 3'' end not shown; two different flanking sequences were used and averaged.', ...
    'Data are normalized based on GAGUA pentaloops in flanking sequences (not shown).', ...
    'DMS, CMCT, and SHAPE normalized so that As, Us, or all 5 loop residues give mean reactivity of 1.0.', ...
    'Note: ', ...
    'Taken in early May 2015 for thirteenth RNA Puzzle community-wide modeling trial.', ...
};

annotations = { ...
    'experimentType:StandardState', ...
    'chemical:Na-HEPES:50mM(pH8.0)', 'chemical:MgCl2:10mM', ...
    'temperature:24C', ...
    'processing:backgroundSubtraction', 'processing:overmodificationCorrection', ...
};
data_annotations{1} = {'modifier:1M7'};
data_annotations{2} = {'modifier:1M7', 'chemical:ZMP:10uM'};
data_annotations{3} = {'modifier:DMS'};
data_annotations{4} = {'modifier:DMS', 'chemical:ZMP:10uM'};
data_annotations{5} = {'modifier:CMCT'};
data_annotations{6} = {'modifier:CMCT', 'chemical:ZMP:10uM'};

% only write ROI to RDAT, excluding flanking sequences
rdat_idx = 25:95;
rdat_seqpos = seqpos_out(rdat_idx);  % 1:71
rdat_sequence = sequence(rdat_idx);
rdat_structure = structure(rdat_idx);
rdat_offset = 0;

% put final profiles together and trim
rdat_d_final = [d_SHAPE_minus, d_SHAPE_plus, d_DMS_minus, d_DMS_plus, d_CMCT_minus, d_CMCT_plus];
rdat_da_final = [da_SHAPE_minus, da_SHAPE_plus, da_DMS_minus, da_DMS_plus, da_CMCT_minus, da_CMCT_plus];
rdat_d_final = rdat_d_final(rdat_idx, :);
rdat_da_final = rdat_da_final(rdat_idx, :);

output_workspace_to_rdat_file(filename, name, rdat_sequence, rdat_offset, rdat_seqpos, ...
    rdat_d_final, rdat_structure, annotations, data_annotations, rdat_da_final, ...
    [], [], [], comments);
% a preview figure of RDAT file is saved
d_rdat = show_rdat(filename);


%% STEP 11
% prediction with no data, only thermodynamic parameters
% no need to bootstrap since there is no data to resample
[structure_NA, bpp_NA] = rna_structure(sequence, [], offset, seqpos_out, [], 0);

% run Fold with 100x bootstrap on 1D SHAPE data
% 100x bootstrap means 100x independent resampling from the original data
% the original (complte) data is used once, before bootstrapping
% so total of 101 Fold calls are made
[structure_1D_Fold_SHAPE_minus, bpp_1D_Fold_SHAPE_minus] = rna_structure(sequence, d_SHAPE_minus, offset, seqpos_out, [], 100, 0);
[structure_1D_Fold_SHAPE_plus, bpp_1D_Fold_SHAPE_plus] = rna_structure(sequence, d_SHAPE_plus, offset, seqpos_out, [], 100, 0);

% run ShapeKnot with 100x bootstrap on 1D SHAPE data
% time a single run before kicking off 100x bootstrap, it may be very long
% use tic and toc commands, with bootstrap = 0 (see above why)
[structure_1D_Spkt_SHAPE_minus, bpp_1D_Spkt_SHAPE_minus] = rna_structure(sequence, d_SHAPE_minus, offset, seqpos_out, [], 100, 1);
[structure_1D_Spkt_SHAPE_plus, bpp_1D_Spkt_SHAPE_plus] = rna_structure(sequence, d_SHAPE_plus, offset, seqpos_out, [], 100, 1);

%% STEP 12
output_varna_html('pfl_NA.html', sequence, structure_NA, structure, structure_NA, offset, [], [], [], bpp_NA);

% the 3 structure-related arguments (#3-5) are
% display_structure: RNA secstr for display (drawing backbone in)
% native_structure: reference/native secstr used for secstr comparison.
% modeled_structure: predicted/new secstr used for secstr comparison.
output_varna_html('pfl_1D_Fold_SHAPE_minus.html', sequence, structure_1D_Fold_SHAPE_minus, structure, structure_1D_Fold_SHAPE_minus, offset, [], [], [d_SHAPE_minus; nan(20, 1)], bpp_1D_Fold_SHAPE_minus);
output_varna_html('pfl_1D_Spkt_SHAPE_minus.html', sequence, structure_1D_Spkt_SHAPE_minus, structure, structure_1D_Spkt_SHAPE_minus, offset, [], [], [d_SHAPE_minus; nan(20, 1)], bpp_1D_Spkt_SHAPE_minus);

% the reactivity vector should be the same length with sequence and structure
% thus 20 of NaN are appended and it will be displayed as gray
output_varna_html('pfl_1D_Fold_SHAPE_plus.html', sequence, structure_1D_Fold_SHAPE_plus, structure, structure_1D_Fold_SHAPE_plus, offset, [], [], [d_SHAPE_plus; nan(20, 1)], bpp_1D_Fold_SHAPE_plus);
output_varna_html('pfl_1D_Spkt_SHAPE_plus.html', sequence, structure_1D_Spkt_SHAPE_plus, structure, structure_1D_Spkt_SHAPE_plus, offset, [], [], [d_SHAPE_plus; nan(20, 1)], bpp_1D_Spkt_SHAPE_plus);

