function [ xsel_all, d_align ] = guess_all_peaks_2ndTRY( d, xsel );

ymin = round(min(min( xsel))) - 50;
ymax = round(max(max( xsel))) + 50;


numlanes = size( d, 2 );

xsel_guess = xsel;

% There may be a guess for every peak...
if( size( xsel_guess,2) ~= numlanes )
  if size( xsel,1) == 1
    xsel = xsel';
  end
  xsel_guess = repmat( xsel, numlanes );
end
num_xsel = size( xsel, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop.

gap_probability = 0 * xsel;

[xsel_all, insertions_all, gaps_all, gap_index_all ] = align_all_peaks( ...
    d, xsel, gap_probability );

[xsel, gap_probability] = refine_after_first_round( xsel, xsel_all, gap_index_all );

% One more iteration, since we now have a sense for where gaps
% might be.
[xsel_all, insertions_all, gaps_all, gap_index_all ] = align_all_peaks( ...
    d, xsel, gap_probability );


%%%%%%%%%%%%%%
% Plot results.
%clf
subplot(1,2,1);
image( d );
colormap( 1 - gray(100))
hold on
plot( [xsel_all]','r' );
for n = 1:numlanes

  if length( gaps_all{n} > 0 )
    plot( n, gaps_all{n},'mx');
  end
  if length( insertions_all{n} > 0 )
    plot( n, insertions_all{n},'bo' );
  end
end
hold off
axis( [ 0.5 numlanes+0.5 ymin ymax ] );
title( 'Peak-assignments to data' );

subplot(1,2,2);
fprintf(1,'Generating realigned image... \n');
d_align = realign_data( d, xsel_all );
title( 'Realigned data');

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xsel_new, insertions, gaps, gap_index] = ...
    align_peaks_dynamic_programming( ...
	xsel, peaks, gap_probability, ...
	d_ref, d_ali );

xsel_new = 0 * xsel;

num1 = length( xsel );
num2 = length( peaks );

DP_scores = zeros( num1+1, num2+1 );
last_align1 = zeros( num1+1, num2+1 );
last_align2 = zeros( num1+1, num2+1 );

width = mean( abs(xsel(1:(num1-1)) - xsel(2:num1)) );
GAP_PENALTY       = 2.0;
INSERTION_PENALTY = 2.0;
INSERTION_PENALTY_EDGES = 0.2;

% Fill in dynamic programming score matrix.
% Note the offset of +1 since the initial row and column is all zeros.
for i = 1 : num1
  for j = 1 : num2

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aligned?
    score_aligned = DP_scores( i, j );
    score_aligned = score_aligned + ...
	(abs(( xsel(i) - peaks(j) )/width))^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try to figure out a penalty based on desired peak-to-peak separation
    % Following seems to slow down the code a lot, but maybe worth it?
    pos1 = last_align1( i, j); 
    pos2 = last_align2( i, j);
    if ( pos1 > 0 & pos2 > 0 )
      score_aligned = score_aligned + ...
	  2 * abs( ( xsel(i) - xsel(pos1)  )  - ...
		     ( peaks(j) - peaks(pos2) ) )  / ...
	  ( width * ( i - pos1 ) )^2;
      
      % And also try a score based on matching peak intensity.
      score_aligned = score_aligned + ...
	  0.2 * abs( log(  d_ref( round(xsel(i)) ) / d_ref( round(xsel( pos1)) ) ) - ...
		     log(  d_ali( round(peaks(j)) ) / ...
			   d_ali( round(peaks( pos2 ) ) ) ) )^2;
    
    end

    %score_aligned = score_aligned + ...
%	min( 0.5*abs( d_ref( round(xsel(i)) ) -...
%		       d_align( round( peaks(j)) ) )/d_ref(round(xsel(i))), 0.5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gap?
    score_gap = DP_scores( i, j + 1) + ...
	GAP_PENALTY * ( 1 - gap_probability(i) ) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insertion?
    if (  ( peaks(j) > xsel(1) )  & ( peaks(j) < xsel(end-1) ) )
      insertion_penalty_local = INSERTION_PENALTY;
    else
      insertion_penalty_local = INSERTION_PENALTY_EDGES;
    end
    score_insertion = DP_scores( i + 1, j) + insertion_penalty_local;

    %if ( j > 1 & j < num2)
    %  score_insertion = score_insertion + min( 0.2 * insertion_penalty_local * ...
    %					       ((peaks(j+1) - peaks(j-1) )/width)^2, 1.0);
    %    end
    
    [ DP_scores(i+1, j+1 ), choice ]   = min( [ score_aligned, score_gap, ...
		    score_insertion ] );
    which_choice( i+1, j+1 ) = choice;
  
    switch choice
     case 1
      last_align1( i+1, j+1 ) = i;
      last_align2( i+1, j+1 ) = j;
     case 2
      last_align1( i+1, j+1 ) = last_align1( i, j+1);
      last_align2( i+1, j+1 ) = last_align2( i, j+1);
     case 3
      last_align1( i+1, j+1 ) = last_align1( i+1, j);
      last_align2( i+1, j+1 ) = last_align2( i+1, j);
    end
    
  end
end

for i = 1:num1
  which_choice( i, 1 ) = 2; % gap
end
for j = 1:num2
  which_choice( 1, j ) = 3; % insertion
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Back-track.
pos1 = num1+1;
pos2 = num2+1;

gap_index = [];
insertions = [];

%image( which_choice*10 );
%image( DP_scores );
%%hold on
while( pos1 > 1  & pos2 > 1 )

  choice = which_choice( pos1, pos2 );
  %  [pos1 pos2 choice]

  switch choice
   case 1
      xsel_new( pos1 - 1 ) = peaks( pos2 - 1 ); % aligned!
      pos1 = pos1 - 1;
      pos2 = pos2 - 1;
   case 2
      % its a gap -- need to fill this in later.
      gap_index = [ gap_index (pos1-1) ];
      pos1 = pos1 - 1;
   case 3
      % its an insertion
      insertions = [ insertions peaks( pos2 - 1 ) ];
      pos2 = pos2 - 1;
  end
  %plot( pos1, pos2,'x' );
end
%hold off
%pause;

%image( DP_scores )
%pause;

% Interpolate linearly to estimate positions of desired bands
% that were not found by automated peak fitting
assigned_index = setdiff( 1:num1, gap_index );

if (length( assigned_index ) > 4 )
  gaps = interp1( assigned_index, xsel_new( assigned_index ), gap_index, ...
		  'linear','extrap' );
  xsel_new( gap_index ) = gaps;
else
  fprintf( 1, ['Problem assigning peaks -- reverting to user-defined' ...
	       ' xsel!\n']);
  
  xsel_new = xsel;
  gaps = xsel;
  gap_index = 1: length( xsel );
  insertions = [];
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xsel_all, insertions_all, gaps_all, gap_index_all ] = ...
    align_all_peaks( d, xsel, gap_probability );

ymin = round(min(min( xsel))) - 100;
ymax = round(max(max( xsel))) + 100;

numlanes = size( d, 2 );

for n = 1:numlanes

  buff = d(:,n);

  % Find peaks in profile.
  peaks = localMaximum(buff(ymin:ymax),10, false);
  peaks = peaks + ymin - 1;

  PLOT_STUFF = 0;
  if PLOT_STUFF
    plot( buff );
    hold on    
    for k = 1:length( peaks )
      plot( [peaks(k) peaks(k)],[0 buff(round(peaks(k)))],'k-' );
    end
    for k = 1:length( xsel )
      plot( [xsel(k) xsel(k)],[0 buff(round(xsel(k)))],'r-' );
    end
    xlim( [ymin ymax]);
    hold off;
    pause;
  end
  
  fprintf( 1, 'Aligning peaks for lane: %d\n', n );
  [ xsel_fit, insertions, gaps, gap_index ] = align_peaks_dynamic_programming(...
      xsel, peaks, gap_probability, d(:,1), d(:,n) );  
  
  xsel_all(:,n) = xsel_fit;
  insertions_all{ n } = insertions;
  gaps_all{ n } = gaps;
  gap_index_all{ n } = gap_index;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xsel, gap_probability] = refine_after_first_round( xsel, xsel_all, gap_index_all );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xsel = mean( xsel_all, 2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_xsel = size( xsel_all, 1 );
gap_probability_sum = ones(1, num_xsel );

num_lanes = length( gap_index_all );
for n = 1: num_lanes
  gap_probability = zeros(1, num_xsel );
  gap_probability( gap_index_all{n} ) = 1;
  gap_probability_sum = gap_probability_sum + gap_probability;
end

gap_probability = gap_probability_sum/num_lanes;
