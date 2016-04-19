%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Medloop -- DMS
%
% Data set is described in 
% Kladwang, W., Cordero P., and Das, R. (2011) "A mutate-and-map strategy
% accurately infers the base pairs of a 35-nucleotide model RNA".
% RNA 17:522-534.
%
% Updated in Feb. 2013, R. Das to match HiTRACE v2.0
%
% A more detailed run-through of the demo will be made available soon in:
%
% https://docs.google.com/a/stanford.edu/Doc?docid=0AfZX7iE4SSSOZGM4MmpmbXZfMTI0cjJwMjJrZmg&hl=en
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundaries of data.

%Directories to readin.
filenames =  { ...
    '063010_Ann_Medloop_L1_DMS1_Run_3100_2010-06-30_291',...
    '063010_Ann_Medloop_L1_DMS1_Run_3100_2010-06-30_292',...
    '063010_Ann_Medloop_L1_DMS1_Run_3100_2010-06-30_293',...
    '062310_Ann_Medloop_redoQC_Run_3100_2010-06-23_268',...
	     };

d_align = quick_look( filenames );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL: align by piecewise-linear mapping (dynamic programming)
align_blocks = [1:40]; % these are 'mutate/map columns'
d_align_before_more_alignment = d_align;
d_align  = align_by_DP_fine( d_align_before_more_alignment, align_blocks );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definitions of what's in this data set.
sequence = 'GGAACGACGAACCGAAAACCGAAGAAAUGAAGAAAGGUUUUCGGUACCGACCUGAAAACCAAAGAAACAACAACAACAAC';
% the value you add to the sequence position to get the 'conventional numbering' (here the numbering used in the Medloop paper)
offset = -10;  
first_RT_nucleotide = length( sequence ) - 20 + offset; % primer binds to last 20 nucleotides
structure= '..........((((((((((...............))))))))))...................................';

% definition of what to expect in each trace.
for i = 1:40; data_types{i} = 'DMS'; end; % DMS mutate-and-map data.
data_types(41:46) = { 'nomod','nomod','ddGTP','ddATP','ddCTP','ddTTP'}; % controls & sequence-ladders
data_types(47:52) = { 'nomod','nomod','ddGTP','ddATP','ddCTP','ddTTP'}; % controls & sequence-ladders

% option not taken here -- mark mutation positions instead.
%mutpos = [ NaN 1:35, NaN 1 2 3];
%for i = 1:length( mutpos );   data_types{i} = mutpos(i); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interactive sequence annotation. Can actually run auto-annotate --
% select top (unextended) and bottom (extended) bands and hit 'x'.
xsel = []; clf;
[xsel,seqpos] = annotate_sequence( d_align, xsel, sequence, offset, data_types, first_RT_nucleotide, structure );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit to Gaussians
area_peak = fit_to_gaussians( d_align, xsel );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save as an RDAT file.
filename = 'ExampleAnalysisMedloop_NEW.rdat';
name = 'Medloop';
annotations = {'experimentType:MutateMap','chemical:Na-HEPES:50mM(pH8.0)','temperature:24C','modifier:DMS','processing:overmodificationCorrection','processing:backgroundSubtraction'};
for i = 1:length(data_types); data_annotations{i} = { ['modifier:',data_types{i}] }; end;
darea_peak = [];
xsel_refine = [];
mutpos = [];
comments = {'Created with the commands in ExampleAnalysisMedloop_script.m'};

output_workspace_to_rdat_file( filename, name, sequence, offset, seqpos, area_peak, ...
			       mutpos, structure, ...
			       annotations, data_annotations, ...
			       darea_peak, ...
			       d_align, xsel, xsel_refine, comments );

rdat = show_rdat( filename );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a VARNA figure
DMS_data = boxplot_normalize( area_peak(:,1) );
easy_varna_fig( 'test_varna.html', sequence, structure, seqpos, offset, DMS_data)

% now open test_varna.html in firefox or other Java-enabled browser.

