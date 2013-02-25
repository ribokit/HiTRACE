function [xsel,DP_SCORE, choice] = auto_assign_sequence_OLD( image_x, xsel, ...
				      nres, offset, ...
				      marks, mutpos );
% AUTO_ASSIGN_SEQUENCE: (still experimental) automatic assignment of bands, given expected locations of marks
%
% [xsel,DP_SCORE, choice] = auto_assign_sequence( image_x, xsel,nres, offset, marks, mutpos );
%
% (C) R. Das, 2010-2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: this assumes that first residues in sequence show up
% *later* in the profiles. Important when making use of "marks" as
% fiducial markers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numpeaks = nres;
numpixels = size( image_x, 1 );

% Scaling of intensity (peak-finding) contribution to score.
alpha = 0.001; 

% Scaling of landmark finding contribution to score.
beta = 0.01;

% to avoid sensitivity to input normalization:
image_x = image_x / mean( mean( abs( image_x ) ) );
%image_x = max( min( image_x, 2.0) , 0.0 );
%image_x = peakify( 100*image_x );
%image_x = image_x / mean( mean( abs( image_x ) ) );


STEP_SIZE = 1;

% Ideal peak spacing.
% Later need to make this computed on the fly.
IDEAL_DELTAX = 24.0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize.
%intensity_score = sum( image_x, 2 );
intensity_score = image_x(:,2);

BIG_NEGATIVE_NUMBER = -9999999999999;

DP_SCORE = BIG_NEGATIVE_NUMBER * ones( numpixels, numpeaks );

DP_SCORE( :, 1 ) = alpha * intensity_score + ...
    beta * landmark_score( 1, nres, offset, image_x, marks, mutpos );

% Disallow points that are too close to end -- this should decrease the
% time of computation by a lot.
min_pos = 1;
max_pos = numpixels;

min_pos1 = min_pos;
max_pos1 = max_pos - (IDEAL_DELTAX*0.9) * nres;

if length( xsel ) >= 2
  min_pos = round( xsel(1) );
  max_pos = round( xsel(end) );
  IDEAL_DELTAX = (xsel(end) - xsel(1) )/ nres;
  min_pos1 = xsel(1);
  max_pos1 = xsel(1);
end


min_pos1 = max( min( round(min_pos1), max_pos-1 ), min_pos );
max_pos1 = max( min( round(max_pos1), max_pos ),   min_pos+1 );

DP_SCORE( (max_pos1+1):end, 1 ) = BIG_NEGATIVE_NUMBER;
DP_SCORE( 1:(min_pos1-1) ,  1 ) = BIG_NEGATIVE_NUMBER;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill rest of dynamic programming matrix

for j = 2:numpeaks
  
  fprintf( 1, 'Doing peak: %d\n', j );

  prev_points_good = find( DP_SCORE(:,j-1) > BIG_NEGATIVE_NUMBER );
  
  minpixel_prev = min(prev_points_good);
  minpixel = round( minpixel_prev + IDEAL_DELTAX * 0.75 );
  minpixel = max( minpixel, min_pos );

  maxpixel_prev = max(prev_points_good);
  maxpixel = round( maxpixel_prev + IDEAL_DELTAX * 1.5 ); 		 
  maxpixel = min( maxpixel, max_pos );
    
  good_range = [ minpixel:STEP_SIZE:maxpixel ]';
  %good_range = [2:numpixels]';
  
  %[ j minpixel maxpixel ]
  
  % Following could probably be sped up using meshgrid.
  for i = good_range'

    min_possible_previous_pos = max( minpixel_prev+1, ...
    				     round( i - IDEAL_DELTAX * 1.5)   );
    %min_possible_previous_pos = minpixel_prev + 1;
    
    good_range_prev = [ min_possible_previous_pos:(i-1) ]';
    %good_range_prev = [ 1:(i-1) ]';

    score_from_previous_peaks = DP_SCORE( good_range_prev, j-1 );
    score_separation = -1.0 * (( i - good_range_prev - IDEAL_DELTAX ) / IDEAL_DELTAX).^2;
    score_from_previous_peaks = score_from_previous_peaks + score_separation;
  
    % also subtact intensity score at a location *between* this peak and previous one
    midpoints = 0.5 * ( i + good_range_prev );

    %score_from_previous_peaks = score_from_previous_peaks - ...
    %   alpha * interp1( [1:numpixels]', intensity_score, midpoints );

    score_from_previous_peaks = score_from_previous_peaks - alpha * intensity_score( round(midpoints) );
  
    [max_score_from_previous_peaks, max_index ] = max( score_from_previous_peaks );
      
    DP_SCORE( i, j ) = max_score_from_previous_peaks;
    choice( i, j ) = good_range_prev( max_index );
    
  end
  
  %DP_SCORE( good_range,j ) = DP_SCORE( good_range,j ) + alpha * ...
  %    intensity_score( good_range );
  
  ls = landmark_score( j, nres, offset, image_x, marks, mutpos );
  DP_SCORE( good_range,j ) = DP_SCORE( good_range,j ) + alpha * intensity_score( good_range ) + ...
      beta * ls( good_range );


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backtrack.

if length( xsel ) >= 2
  xsel( numpeaks ) = round( xsel( end ) );
else
  [ dummy, xsel( numpeaks ) ] = max( DP_SCORE( :, numpeaks ) );
end

for j = (numpeaks-1):-1:1
  xsel( j ) = choice( xsel(j+1), j+1 );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results...
cla;
image( image_x * 20 );
colormap( 1 - gray(100 ) );

hold on
for j = 1:numpeaks
  plot( [0.5  size(image_x,2)+0.5], [xsel(j) xsel(j)],'r');
end

if(~isempty(marks))
    for m = 1:size(image_x,2);
        if( m <= length(mutpos))
          goodpos = find( marks(:,1) ==  mutpos(m)  );
          for k = goodpos'
            which_xsel = nres - marks( k, 2 ) + 1 + offset;
            if ( which_xsel >= 1 & which_xsel <= length(xsel) )
              plot( m, xsel(which_xsel),'ro' );
            end
          end
        end
    end
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score = landmark_score( peak_number, nres, offset, image_x, ...
				 marks, mutpos )

seqpos = (nres - peak_number + 1) + offset;
if(~isempty(marks))
    gp = find( seqpos == marks(:,2) );
else
    gp = [];
end
score = zeros( size( image_x, 1 ), 1 );

if length( gp ) > 0

  relevant_mutpos = marks( gp,1 );
  
  for m = 1: length( relevant_mutpos )
    relevant_lanes = find( mutpos == relevant_mutpos( m ) );
    relevant_lanes = intersect( relevant_lanes, [1:size(image_x,2)] );
    if length( relevant_lanes ) > 0 
      score = score + abs( mean( image_x( :, relevant_lanes ), 2 ) - ...
			   mean( image_x,2 ) )./std(image_x,0,2);

    end
    %score = smooth( score, 10 );
  end
  
end
