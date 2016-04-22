% prepare reactivity data
d = load('pfl_1D.mat');

d_SHAPE_plus = d.d_SHAPE_plus(25:95);
d_SHAPE_minus = d.d_SHAPE_minus(25:95);

for i = [1:length(d.sequence)] - d.offset;
    if d.sequence(i + d.offset) == 'A' || d.sequence(i + d.offset) == 'C';
        d_DMS_CMCT_plus(i + d.offset) = d.d_DMS_plus(i);
        d_DMS_CMCT_minus(i + d.offset) = d.d_DMS_minus(i);
    else
        d_DMS_CMCT_plus(i + d.offset) = d.d_CMCT_plus(i);
        d_DMS_CMCT_minus(i + d.offset) = d.d_CMCT_minus(i);
    end;
end;


% MATLAB setup
mat_name = 'pfl_1D_color.mat';
imagex = imread('pfl_secstr.tif');

residue_locations = [];
base_locations = [];
square_width = 28;

sequence = 'GGGUCGUGACUGGCGAACAGGUGGGAAACCACCGGGGAGCGACCCCGGCAUCGAUAGCCGCCCGCCUGGGC';
seqpos = 1:71;
offset = 0;


% mark box coordinates
residue_locations = pick_points(imagex, offset, residue_locations, square_width);

% mark stub coordinates
direction_cell = { ...
    [25, 29:33, 46:56, 61:62], ...          % top
    [18:24, 28, 60], ...                    % bottom
    [8:17, 26:27, 39:45, 57:59, 67:68], ... % left
    [1:7, 34:38, 63:66, 69:71], ...         % right
};
base_locations = move_base_locations(residue_locations, direction_cell, offset, square_width / 2, 5);

base_locations = pick_bases(imagex, offset, residue_locations, base_locations, square_width);
save(mat_name);


% get color scheme
whichres = seqpos - offset;
color_scheme = 17;

whattoplot = d_SHAPE_plus;
color_profile = color_palette(whattoplot, 1.5, 0, color_scheme, seqpos, sequence);

% color boxes
imagex_1 = color_residues(imagex, residue_locations, whichres, whattoplot, color_profile, square_width);


% color stubs
whattoplot = d_DMS_CMCT_plus;
color_profile = color_palette(whattoplot, 1.5, 0, colorscheme, seqpos, sequence);
imagex_2 = color_bases(imagex_1, base_locations, residue_locations, whichres, whattoplot, color_profile, square_width / 2);

% add legend
orient_legend = [1 0];
pos_legend = [1 3 5 4];
labels = {'1.5', '0', 'SHAPE  '};
font_size = 0.8;

imagex_3 = color_legend(imagex_2, color_profile, square_width, orient_legend, pos_legend, labels, font_size);


% output diagram
image_output(imagex_3, 'pfl_secstr_color_plus', 300);
save(mat_name);


% more fun
whattoplot = d_SHAPE_minus;
color_profile = color_palette(whattoplot, 1.5, 0, color_scheme, seqpos, sequence);
imagex_1 = color_residues(imagex, residue_locations, whichres, whattoplot, color_profile, square_width);

whattoplot = d_DMS_CMCT_minus;
color_profile = color_palette(whattoplot, 1.5, 0, color_scheme, seqpos, sequence);
imagex_2 = color_bases(imagex_1, base_locations, residue_locations, whichres, whattoplot, color_profile, square_width/2);

imagex_3 = color_legend(imagex_2, color_profile, square_width, orient_legend, pos_legend, labels, font_size);
image_output(imagex_3, 'pfl_secstr_color_minus', 300);


color_scheme = 1;

whattoplot = d_SHAPE_plus - d_SHAPE_minus;
color_profile = color_palette(whattoplot, 1, -1, color_scheme, seqpos, sequence);
imagex_1 = color_residues(imagex, residue_locations, whichres, whattoplot, color_profile, square_width);

whattoplot = d_DMS_CMCT_plus - d_DMS_CMCT_minus;
color_profile = color_palette(whattoplot, 1, -1, color_scheme, seqpos, sequence);
imagex_2 = color_bases(imagex_1, base_locations, residue_locations, whichres, whattoplot, color_profile, square_width/2);

labels = {'1', '-1', 'SHAPE  '};
imagex_3 = color_legend(imagex_2, color_profile, square_width, orient_legend, pos_legend, labels, font_size);
image_output(imagex_3, 'pfl_secstr_color_diff', 300);
save(mat_name);


% vectorized figures
color_scheme = 12   % publication-use red-yellow-white
whattoplot = d_SHAPE_plus;
color_profile = color_palette(whattoplot, 1.5, 0, color_scheme, seqpos, sequence);

color_circles(imagex, residue_locations, whichres, whattoplot, color_profile, square_width, 'circ_SHAPE_plus');


color_scheme = 13   % publication-use red-white-blue
whattoplot = d_SHAPE_plus - d_SHAPE_minus;
color_profile = color_palette(whattoplot, 1, -1, color_scheme, seqpos, sequence);

color_circles(imagex, residue_locations, whichres, whattoplot, color_profile, square_width, 'circ_SHAPE_diff');

