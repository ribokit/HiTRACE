function easy_varna_fig( filename, sequence, structure, seqpos, offset, DATA );
% VARNA_FIG: Create html file with secondary structure -- double click to get VARNA visualizer in a web browser 
%
%  varna_fig(filename,sequence,structure,seqpos, offset, DATA)
%
% filename  = output filename [e.g., 'my_rna_in_varna.html']
% sequence  = RNA sequence
% structure = structure in dot/bracket notation. length should match sequence.
% 
% optional input parameters:
% seqpos            = which residue numbers do the DATA correspond to?
% offset            = value to  add to sequence index to get your favorite sequence numbering
% DATA              = data for coloring residues. will only use values between 0 and 1, 
%                        so normalize ahead of time. Give [] if no data.
%                        If you supply two columns, they will be assumed to be DMS and CMCT,
%                        and the plot will use DMS for A/C, and CMCT for G/U.
%
% (C) R. Das 2011

if length( seqpos ) ~= size( DATA, 1 ); fprintf( 'WARNING! length(seqpos) is not the same size as DATA\n' ); return;end;

DATA_TO_PLOT = NaN * ones( length(sequence ), 1 );

if size( DATA, 2) == 2
  fprintf( 'I am assuming you gave DMS & CMCT!\n' );
  for i = 1:length( seqpos )
    switch sequence( seqpos(i)-offset )
     case {'A','C'}
      DATA_TO_PLOT( seqpos(i)-offset ) = DATA( i, 1 );
     case {'G','U','T'}
      DATA_TO_PLOT( seqpos(i)-offset ) = DATA( i, 2 );     
    end  
  end
else
  DATA_TO_PLOT( seqpos - offset ) = DATA;
end

varna_fig( filename, sequence, structure, DATA_TO_PLOT, 2, offset );
