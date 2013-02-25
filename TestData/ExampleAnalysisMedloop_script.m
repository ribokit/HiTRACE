%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Medloop -- DMS
%
% Data set is described in 
% Kladwang, W., Cordero P., and Das, R. (2011) "A mutate-and-map strategy
% accurately infers the base pairs of a 35-nucleotide model RNA".
% RNA 17:522-534.
%
% Updated in Jan. 2012, R. Das.
% Updated again in Feb. 2013, R. Das to match HiTRACE v2.0
%
% Check out more detailed run-through of the demo in:
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

% advanced: you can specify the top and bottom of the profiles by ymin and ymax,
%  as well as specify a subset of profiles (perhaps reordered). For example:
%
%ymin = 2200; 
%ymax = 4200;  
%reorder = [ 1:52 ]; 
%[d_align, d_ref_align] = quick_look( filenames, ymin, ymax, reorder );
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL: baseline subtract removes any smooth offsets in signals.
d_align_nobsub = d_align;
d_align = baseline_subtract_v2( d_align_nobsub, 1, size( d_align_nobsub,1) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL: align by piecewise-linear mapping (dynamic programming)
align_blocks = [1:40];
d_align_before_more_alignment = d_align;

d_align1  = align_by_DP_fine( d_align_before_more_alignment, align_blocks );
d_align = d_align1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL manual alignment -- kept here for completeness.
%d_old = d_bsub_cut; anchorlines = [];
%[d_manual_align, anchorlines] = align_userinput_saveanchorlines(100*d_old,1,1,anchorlines);
%d_bsub_cut = d_manual_align;
%d_align  = align_by_DP( d_bsub_cut, align_blocks );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definitions of what's in this data set.
sequence = 'GGAACGACGAACCGAAAACCGAAGAAAUGAAGAAAGGUUUUCGGUACCGACCUGAAAACCAAAGAAACAACAACAACAAC';
% the value you add to the sequence position to get the 'conventional numbering' (here the numbering used in the Medloop paper)
structure= '..........((((((((((...............))))))))))...................................';
offset = -10;  
primer_binding_site = length( sequence ) - 20 + 1 + offset;

% What is in this data set?
for i = 1:40; data_types{i} = 'DMS'; end;
data_types(41:46) = { 'nomod','nomod','ddGTP','ddATP','ddCTP','ddTTP'};
data_types(47:52) = { 'nomod','nomod','ddGTP','ddATP','ddCTP','ddTTP'};

% option not taken here -- mark mutation positions instead.
%mutpos = [ NaN 1:35, NaN 1 2 3];
%for i = 1:length( mutpos );   data_types{i} = mutpos(i); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interactive sequence annotation
xsel = []; clf;
[xsel,seqpos] = annotate_sequence( d_align, xsel, sequence, offset, data_types, primer_binding_site, structure );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit to Gaussians
area_peak = do_the_fit_fast( d_align, xsel );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Background subtraction and correction for 'overmodification', a.k.a., attenuation of longer products
backgd_col = [41];
norm_bins = [10 : (size(area_peak,1)-10) ];
[area_bsub, darea_bsub]  = overmod_and_background_correct_logL( area_peak, backgd_col, norm_bins );

image( area_bsub * 40)
errorbar( area_bsub(:,1), darea_bsub(:,1) ); ylim([-0.5 2] )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save as an RDAT file.
filename = 'ExampleAnalysisMedloop_NEW.rdat';
name = 'Medloop';
annotations = {'experimentType:MutateMap','chemical:Na-HEPES:50mM(pH8.0)','temperature:24C','modifier:DMS','processing:overmodificationCorrection','processing:backgroundSubtraction'};
for i = 1:length(data_types); data_annotations{i} = { ['modifier:',data_types{i}] }; end;
xsel_refine = [];
comments = {'Created with the commands in ExampleAnalysisMedloop_script.m'};
	    
output_workspace_to_rdat_file( filename, name, sequence, offset, seqpos, area_bsub, ...
			       [], structure, ...
			       annotations, data_annotations, ...
			       darea_bsub, ...
			       d_align, xsel, xsel_refine, comments );

rdat = show_rdat( filename );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Averaging across data sets.

rdat  = show_rdat( 'ExampleAnalysisMedloop.rdat' );
wt_mean = mean(rdat.reactivity(:,[1 37])');
wt_err  = std(rdat.reactivity(:,[1 37])');

rdat0 = show_rdat( 'ExampleAnalysisMedloop_previousDMSreplicate.rdat' );
wt_mean0 = mean(rdat0.reactivity(:,[1 37])');
wt_err0  = std(rdat0.reactivity(:,[1 37])');

clf;
errorbar( rdat.seqpos, wt_mean, wt_err ); hold on
errorbar( rdat0.seqpos, wt_mean0, wt_err0, 'r'); 
d_mean = {wt_mean, wt_mean0 };
d_err = {wt_err, wt_err0 };
[d_mean_final, d_err_final] = get_average_standard_state( d_mean, d_err );

errorbar( rdat0.seqpos, d_mean_final, d_err_final, 'k'); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a VARNA figure
easy_varna_fig( 'test_varna.html', sequence, structure, seqpos, offset, d_mean_final )