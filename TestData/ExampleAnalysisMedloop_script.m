%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Medloop -- DMS
%
% Data set is described in 
% Kladwang, W., Cordero P., and Das, R. (2011) "A mutate-and-map strategy
% accurately infers the base pairs of a 35-nucleotide model RNA".
% RNA 17:522-534.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundaries of data.
ymin = 2200; 
ymax = 4200;

%Directories to readin.
filenames =  { ...
    '063010_Ann_Medloop_L1_DMS1_Run_3100_2010-06-30_291',...
    '063010_Ann_Medloop_L1_DMS1_Run_3100_2010-06-30_292',...
    '063010_Ann_Medloop_L1_DMS1_Run_3100_2010-06-30_293',...
    '../../Data/062310_Ann_Medloop_redoQC_Run_3100_2010-06-23_268',...
	     };

reorder = [ 1:52 ]; % or could just not specify -- by default all profiles get read in.
[d,da] = quick_look( filenames, ymin, ymax, reorder );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% baseline subtract removes any smooth offsets in signals.
d_bsub = baseline_subtract_v2( d, ymin, ymax, 2e6, 2e4 );
d_bsub_cut = d_bsub( ymin:ymax, : );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align by warping (dynamic programming)
align_blocks = { [41 42 1], [1:40] };
d_align  = align_by_DP( d_bsub_cut, align_blocks );

% optional manual alignment -- kept here for completeness.
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

% What is in this data set?

% option not taken here -- mark mutation positions
%data_type_mutpos = [ NaN 1:35, NaN 1 2 3];
%for i = 1:length( data_type_mutpos);   data_type{i} = num2str( data_type_mutpos( i ) ); end;

for i = 1:40; data_type{i} = 'DMS'; end;
data_type(41:46) = { 'nomod','nomod','ddGTP','ddATP','ddCTP','ddTTP'};
data_type(47:52) = { 'nomod','nomod','ddGTP','ddATP','ddCTP','ddTTP'};

seqpos = [ (length(sequence) - 20): -1 : 1] + offset;
[ marks, area_pred, mutpos] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence, offset , seqpos, data_type );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interactive sequence annotation
xsel = []; period = 1;
sequence_actual = sequence( 1: end-20 );
xsel = mark_sequence( d_align, xsel, sequence_actual, 0, offset, period, marks,mutpos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit to Gaussians
area_peak = do_the_fit_FAST( d_align, xsel );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Background subtraction and correction for 'overmodification', a.k.a., attenuation of longer products
backgd_col = [41];
norm_bins = [10 : (size(area_peak,1)-10) ];
[area_bsub, darea_bsub]  = overmod_and_background_correct_by_LP( area_peak, backgd_col, norm_bins, area_pred );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save as an RDAT file.
filename = 'ExampleAnalysisMedloop.rdat';
name = 'Medloop';
annotations = {'experimentType:MutateMap','chemical:Na-HEPES:50mM(pH8.0)','temperature:24C','modifier:DMS','processing:overmodificationCorrection','processing:backgroundSubtraction'};
for i = 1:length(data_type); data_annotations{i} = { ['modifier:',data_type{i}] }; end;
xsel_refine = [];
comments = {'Created with the commands in ExampleAnalysisMedloop_script.m'};
	    
output_workspace_to_rdat_file( filename, name, sequence, offset, seqpos, area_bsub, ...
			       mutpos, structure, ...
			       annotations, data_annotations, ...
			       darea_bsub, ...
			       d_align, xsel, xsel_refine, comments );

rdat = show_rdat( filename );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Averaging across data sets.

rdat  = show_rdat( 'ExampleAnalysisMedloop.rdat' );
wt_mean = mean(rdat.area_peak(:,[1 37])');
wt_err  = std(rdat.area_peak(:,[1 37])');

rdat0 = show_rdat( 'ExampleAnalysisMedloop_previousDMSreplicate.rdat' );
wt_mean0 = mean(rdat0.area_peak(:,[1 37])');
wt_err0  = std(rdat0.area_peak(:,[1 37])');

errorbar( rdat.seqpos, wt_mean, wt_err ); hold on
errorbar( rdat0.seqpos, wt_mean0, wt_err0, 'r'); hold off
d_mean = {wt_mean, wt_mean0 };
d_err = {wt_err, wt_err0 };
[d_mean_final, d_err_final] = get_average_standard_state( d_mean, d_err );

 
errorbar( rdat0.seqpos, d_mean_final, d_err_final, 'k'); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a VARNA figure
easy_varna_fig( 'test_varna.html', sequence, structure, seqpos, offset, d_mean_final )