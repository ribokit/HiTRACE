function []=quick_varna(filename,sequence,structure,DATA,offset,seqpos,colorscheme);
% VARNA_FIG: Create html file with secondary structure -- double click to get VARNA visualizer in a web browser 
%
%  varna_fig(filename,sequence,structure,DATA,seqpos);
%
% filename  = output filename [e.g., 'my_rna_in_varna.html']
% sequence  = RNA sequence
% structure = structure in dot/bracket notation. length should match sequence.
% 
% optional input parameters:
% DATA               = data for coloring residues. will only use values between 0 and 1, 
%                        so normalize ahead of time. Give [] if no data.
% offset             = integer offset to get 'conventional numbering' -- first position of sequence should be this number+1.
% seqpos             = which nucleotide numbers each datum in DATA correspond to.
%
% (C) R. Das 2012

if ~exist( 'colorscheme' ); colorscheme = 1; end;
if ~exist( 'offset' ); offset = 1; end;
if ~exist( 'seqpos' ); seqpos = [1:length(data)]+offset;end

data_new_numbering = zeros(1, length( sequence ) );
for i = 1:length( DATA );
  seqnum = seqpos(i)-offset;
  if ( seqnum < 1 | seqnum > length(sequence) ); fprintf( 'Bad seqpos %d, corresponds to %d-th nucleotide in sequence?\n', seqpos(i), seqnum )  ; end;
  data_new_numbering( seqnum ) = DATA(i);
end
varna_fig( filename, sequence, structure, data_new_numbering, colorscheme);
