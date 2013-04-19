function [d_correct,modification_fractions] = correct_for_attenuation( d_in, seqpos )
% [d_correct,modification_fractions] = correct_for_attenuation( d_in )
%
% d_in = array of band intensities, ordered from 5' to 3'. 
%        First position should correspond to fully extended cDNA.
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'seqpos', 'var' ) seqpos = [1 : size( d_in, 1 )]; end;

for k = 1:size(d_in,2);
  
  d = d_in(:,k);

  % On calculating modification fractions, 
  % skip first value, as it might be unextended primer
  % also ignore data right at end [leave out last 3 nts]
  modification_fractions(k) = sum(d(5 : end-1)) / sum( d(1:end-1) );
  fprintf( '%2d modification fraction: %8.4f\n', k, modification_fractions(k));
  
  d_ratios(:,k) = 0 * d;
  for j = 1:length( d );
    d_ratios(j,k) =  d(j)/sum( d(1:end) );
  end;

  d_correct(:,k) = 0 * d;
  d_correct(end,k) = 1;
  for j = 1:length(d_correct) - 1;
    d_correct(j,k) = d(j) / sum(d(1:j));
  end;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some visual feedback
clf
subplot(1,2,1);
image( [1:size(d_ratios,2)], seqpos, d_ratios*2000 );
title( 'Before  correct for attenuation' );
make_lines;
subplot(1,2,2);

image( [1:size(d_ratios,2)], seqpos, d_correct*2000 );
title( 'After  correct for attenuation' );
colormap( 1 - gray(100) );
make_lines;

return;
