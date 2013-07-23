function [xsel, best_score, msg] = auto_assign_sequence( image_x, sequence, area_pred, ideal_spacing, input_bounds, data_types )
% AUTO_ASSIGN_SEQUENCE: (still experimental) automatic assignment of bands, given expected locations of marks
%
% [xsel, best_score, msg] = auto_assign_sequence( image_x, sequence, offset, area_pred, ideal_spacing, input_bounds, data_types )
%
%
%  image_x        = matrix of aligned electrophoretic traces.
%  sequence       = sequence corresponding to each band [from top to bottom! may be reverse of actual sequence]
%  area_pred      = strength of bands expected at each position [in a range from 0 to 1 for no band to visible band]
%  ideal_spacing  = [optional] peak-to-peak spacing. [default 0 will trigger estimation of ideal_spacing from auto correlation function]
%  input_bounds   = [optional] start band position. If given, the last band position. If number of position corresponds to sequence size, 
%                     the scripts will not reassign sequence -- will just compute score. 
%  data_types     = [optional] cell of data type strings [e.g., {'DMS','ddATP', etc.} -- was used to upweight no mod lanes, and
%                     to downweight reference ladders, but this is not the case anymore].
%
% (C) R. Das, Hanjoo Kim, Sungroh Yoon, 2010-2011
% (C) R. Das, 2013.

if nargin == 0;  help( mfilename ); return; end;

[num_pixels, num_lanes] = size( image_x );
if ~exist('ideal_spacing','var') | isempty(ideal_spacing)  | ideal_spacing == 0

  max_ideal_spacing = floor(num_pixels / (size( area_pred, 1 )-1) );  
  ideal_spacing = guess_ideal_spacing( image_x, max_ideal_spacing );
  if ( ideal_spacing == 0 );  
    fprintf( 'Guess of ideal_spacing failed\n' );
    ideal_spacing = max_ideal_spacing;
    fprintf( 'Using ideal_spacing: %8.3f\n',  ideal_spacing );
  end;
end;
if ~exist('input_bounds','var'), input_bounds = []; end
if ~exist( 'data_types', 'var' ), data_types = []; end;


exist_talepeak = 1;
%FALSEPEAK_CUT = 8;
%[exist_talepeak, image_x, msg] = remove_tale_and_false_peaks( image_x, area_pred, FALSEPEAK_CUT );

% note that this was from a time when we used to list sequence 
% positions backward. Probably should fix this...
s = area_pred;
nres = length( sequence );
s(1,:)    = 1;
%s(nres,:) = 1;
sequence_at_bands = sequence( end:-1:1 );
if exist_talepeak
  s(nres+1,:) = 2;
  sequence_at_bands = [sequence_at_bands, 'X'];
end
%s( s == 0.0 ) = 0.01;
s( s == 0.0 ) = 0.2;

tic
[xsel, best_score] = solve_xsel_by_DP( image_x, s, sequence_at_bands, ideal_spacing, input_bounds, data_types );
toc % read out time required for fit.

if exist_talepeak
  fprintf( 'Assigned fully extended band... one before nominal beginning of sequence.\n' );
  %    xsel(end) = []; % keep this band!
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ideal_spacing = guess_ideal_spacing( image_x, max_spacing )
% look for closely spaced peaks using autocorrelation
max_spacing_for_corr = floor( 1.5*max_spacing );
possible_spacings = [];
for i = 1:size( image_x, 2 )
  autocorr_profile(:,i) = autocorr( image_x(:,i), max_spacing_for_corr );
  spacings = localMaximum( autocorr_profile(:,i) );

  % the first peak is always at 1 -- ignore that one and go to the next tightest spacing
  if length( spacings ) > 1 & spacings(2) < max_spacing;
    possible_spacings = [ possible_spacings, spacings(2) ];
  end
end

if length( possible_spacings ) == 0
  ideal_spacing = 0; % no clue.
else
  %ideal_spacing = median( possible_spacings ); 
  ideal_spacing = median( possible_spacings ); 
end

fprintf( 'Trying to guess ideal peak-to-peak spacing: %8.1f\n', ideal_spacing );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [exist_talepeak, image_x, msg] = remove_tale_and_false_peaks( image_x, area_pred, FALSEPEAK_CUT );

exist_talepeak = true;
msg = [];


% getting rid of data at and beyond 'tale peaks' [fully extended band] which
% can dominate fit.
[num_pixels, num_lanes] = size( image_x );
gee = zeros(1,num_lanes);
spread = floor(num_pixels / (size( area_pred, 1 )-1) )/2;
for i = 1:num_lanes
  % scale image
  scalefactor = mean( image_x(:,i) );
  image_x(:,i) = image_x(:,i)/scalefactor/2;
    
  % find all peaks.
  peaks = localMaximum( image_x(:,i), spread );
  peaks = peaks( image_x(peaks,i) > median(image_x(peaks,i))/2 );

  % get maximum peak of last five peaks. originally went to peak right before this one [ - 1].
  tmp = find( image_x(peaks,i) == max( image_x(peaks(end-5:end),i)) ) - 1;
  %tmp = find( image_x(peaks,i) == max( image_x(peaks(end-5:end),i)) );

  if (image_x(peaks(tmp),i) < FALSEPEAK_CUT/12 && num_pixels > 10000)
    gee(i) = num_pixels;
    msg = 'a tale peak is not found';
    exist_talepeak = false;
  else
    gee(i) = peaks(tmp);
  end
end

if (max(gee) < num_pixels)
    talepeakrange = (1:num_pixels) > (max(gee)+1);
    image_x(talepeakrange,:) = 0.0;
end

% Remove 'false peaks'
n = num_lanes - 1;
falsepeaks = find(image_x(peaks(1:tmp),n) > FALSEPEAK_CUT);
falsepeaks((peaks(falsepeaks) < 0.25*num_pixels) | (peaks(falsepeaks) > 0.8*num_pixels)) = [];
if ~isempty(falsepeaks)
    falsepeakindex = peaks(falsepeaks(end));
    falsepeakrange = intersect(find(image_x(:,n) > FALSEPEAK_CUT/6), falsepeakindex-30:falsepeakindex+30);
    if ~isempty(msg)
        msg = [msg ' & '];
    end
    msg = [msg 'detected false peak at (' num2str(min(falsepeakrange)) '~' num2str(max(falsepeakrange)) ') in baseline signal'];
    image_x(falsepeakrange,:) = 0;
end
