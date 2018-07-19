function [d_correct,d_correct_err,modification_fractions] = correct_for_attenuation( d_in, d_in_err, seqpos )
% CORRECT_FOR_ATTENUATION: takes band intensities from reverse transcription and gives back
%                           fraction of each position that was modified. 
%
% [d_correct, d_correct_err, modification_fractions] = correct_for_attenuation( d_in, d_in_err )
%
% d_in = array of band intensities, ordered from 5' to 3'. 
%        First position should correspond to fully extended cDNA.
% d_in_err = [optional] array of errors corresponding to d_in
% seqpos   = [optional] sequence positions (used in plotting).
%
% Correction formula is exact:
%
%  R(i) = F(i) / [ F(0) + F(1) + ... F(i) ]
%
% i.e. the fraction of reverse transcription that stops at residue i, compared
% to the total number of cDNAs that have been reverse transcribed up to 
% residue i or beyond
%
% (C) T. Mann, S. Tian, R. Das, 2012-2013.

if nargin < 1;  help( mfilename ); return; end;

if ~exist( 'd_in_err', 'var' ) d_in_err = 0 * d_in; end;
if ~exist( 'seqpos', 'var' ) seqpos = [1 : size( d_in, 1 )]; end;

for k = 1:size(d_in,2);
  
  d     = d_in(:,k);
  d_err = d_in_err(:,k);

  % On calculating modification fractions, 
  % skip last few values, as they might be unextended primer
  % also ignore data right at beginning [large peak from fully extended primer]
  minpos = min( 5, length(d) );
  maxpos = max( length(d)-5, 1 );
  sum_mod = sum( d(minpos : maxpos));
  sum_all = sum( d(1:maxpos) );
  modification_fractions(k) = sum_mod / sum_all;
  fprintf( '%2d modification fraction: %8.4f\n', k, modification_fractions(k));
  d_ratios(:,k) =  d/sum_all;

  % note that modification fractions aren't actually used to calculated
  % reactivities.
  d_correct(:,k) = 0 * d;
  d_correct_err(:,k) = 0 * d;
  for j = 1:length(d);
    d_correct(j,k)     = d(j)     / sum(d(1:j));
    d_correct_err(j,k) = d_err(j) / sum(d(1:j));
  end;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some visual feedback
clf;
set_print_page(gcf, 0, [0 0 500 800], 'Attenuation Correction');

subplot(1,2,1);
image( [1:size(d_ratios,2)], seqpos, d_ratios*2000 );
title( 'BEFORE', 'FontSize', 11, 'FontWeight', 'Bold' );
if ( size( d_correct, 2 ) < 100 ) make_lines; end;
subplot(1,2,2);

image( [1:size(d_ratios,2)], seqpos, d_correct*2000 );
title( 'Attenuation Corrected', 'FontSize', 11, 'FontWeight', 'Bold' );
colormap( 1 - gray(100) );
if ( size( d_correct, 2 ) < 100 ) make_lines; end;

return;
