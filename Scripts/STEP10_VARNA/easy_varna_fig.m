function easy_varna_fig( filename, sequence, structure, seqpos, offset, DATA );

DATA_TO_PLOT = NaN * ones(1, length(sequence ) );
DATA_TO_PLOT = DATA( seqpos - offset );
varna_fig( filename, sequence, structure, DATA_TO_PLOT, 2 );
