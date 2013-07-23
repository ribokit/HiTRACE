function  [subtracted, subtracted_err ] = subtract_array( signal, background, signal_err, background_err, seqpos, image_scalefactor )
% SUBTRACT_ARRAY: Subtracts two arrays, and also progates errors. Useful for background subtaction.
%
%  [subtracted, subtracted_err ] = subtract_array( signal, background, signal_err, background_err, seqpos, image_scalefactor )
%
% Inputs:
%  signal     = array with signal areas/reactivities.
%  background = array with background areas/reactivities.
%
% Optional inputs:
%  signal_err     = array with errors on signal
%  background_err = array with errors on background
%  seqpos         = residue numberings that go with positions [for plotting; default is 0, 1, 2,...].
%  image_scalefactor = factor by which to multiply array by plotting [default is calculated based on mean intensity]
%
% (C) R. Das, 2013.
%

if ( nargin < 2 ); help( mfilename); return; end;

if ~exist( 'signal_err', 'var' ); signal_err = 0 * signal; end;
if ~exist( 'background_err', 'var' ); background_err = 0 * background; end;
if ~exist( 'seqpos', 'var' ) seqpos = [0 : size( signal, 1 ) - 1]; end;

subtracted     = signal - background;
subtracted_err = sqrt( signal_err.^2 + background_err.^2 );
  
ntrace = size( signal, 2 );
if ~exist( 'image_scalefactor', 'var' ); 
  image_scalefactor = 40/mean(mean( subtracted ) );
end;

image( 1:ntrace, seqpos, subtracted * image_scalefactor );
title( 'Background Corrected', 'FontSize', 11, 'FontWeight', 'Bold');
if ( ntrace < 100 ) make_lines; end

fprintf('Press any key to continue ...\n'); 
pause;

