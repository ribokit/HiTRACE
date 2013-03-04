function [data_align, x_realign] = align_to_first_NEW( data, PLOT_STUFF, refcol ); 

PLOT_STUFF = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PEAK_SEP = 20; %minimum peak separation, in pixels.
MIN_CUTOFF = 0.5 * std( data(:,refcol) ); % minimum intensity for a peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_pixels = size( data, 1);
x = 1:num_pixels;

% No need to look at data if lower than a certain cutoff.
[pixels_ref, d_ref] = get_profile( data, refcol, MIN_CUTOFF ) ;
% Find peaks.
peak_ref = get_peaks( d_ref, PEAK_SEP, MIN_CUTOFF );
if PLOT_STUFF  & 0
  subplot(2,1,1);  make_peak_plot( d_ref, peak_ref );
end

% Main loop.
numlanes = size( data, 2);
for n = 1: numlanes

  fprintf(1,'Peak-based alignment for ... %d\n',n);
  % Generate similar peak list for data to be aligned.
  [pixels_ali, d_ali] = get_profile( data, n, MIN_CUTOFF ) ;

  pixels_ali = pixels_ali+30;
  d_ali = data(pixels_ali,n);
  
  peak_ali = get_peaks( d_ali, PEAK_SEP, MIN_CUTOFF );

  if PLOT_STUFF & 0
    subplot(2,1,2);  make_peak_plot( d_ali, peak_ali );
    pause;
  end

  % Align by dynamic programming!
  %peak_ref = peak_ref(1:2);
  %peak_ali = peak_ali(1:4);
  %peak_ref(2) - peak_ref(1)
  %peak_ali(3:4) - peak_ali(2)
    [ corresponding_peaks, ...
    unassigned_ref, ...
    unassigned_ali ] = ...
      align_peaks( d_ref, d_ali, peak_ref, peak_ali );

  if PLOT_STUFF 
    plot_alignment( d_ali, d_ref, corresponding_peaks, unassigned_ref, unassigned_ali,...
		    peak_ref, peak_ali);
  end
  
  peak_pos1 =  corresponding_peaks(1,:) + min( pixels_ref) - 1;
  peak_pos2 =  corresponding_peaks(2,:) + min( pixels_ali) - 1;
  
  p = polyfit( peak_pos1, peak_pos2, 1 );
  new_x = p(1) * x + p(2);
  %new_x = interp1( peak_pos1, peak_pos2, x, 'linear','extrap');
  
  if PLOT_STUFF 
    plot( peak_pos1, peak_pos2, 'x');
    hold on
    plot( pixels_ref, new_x(pixels_ref),'r-' );
    hold off
    pause;
  end
  
  da = interp1( x, data(:,n), new_x, 'linear',0);

  data_align(:,n) = da;
  x_realign( :, n ) = new_x;
  
end
subplot( 1,2,1)
image( data(pixels_ref,:) )
subplot(1,2,2);
image( data_align(pixels_ref,:) )

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pixels_ref, d_ref ] = get_profile( data, refcol, MIN_CUTOFF )

d_ref = data(:,refcol );
goodpoints = find( d_ref > MIN_CUTOFF );

pixels_ref = min( goodpoints ): max(goodpoints );
d_ref = max( d_ref( pixels_ref ), MIN_CUTOFF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function peak_ref = get_peaks(  d_ref, PEAK_SEP, MIN_CUTOFF );
peak_ref = localMaximum( d_ref, PEAK_SEP, true );
peak_ref = peak_ref( find( d_ref(peak_ref) > 1.2 * MIN_CUTOFF ) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_peak_plot( d_ref, peak_ref );
plot( d_ref );
hold on
for k = 1: length( peak_ref )
  plot( peak_ref(k)*[1 1], [0 d_ref( round(peak_ref(k)))],'r-' );
end
hold off  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_alignment( d_ali, d_ref, corresponding_peaks, ...
			 unassigned_ref, unassigned_ali , ...
			 peak_ref, peak_ali);

clf
imagex = zeros( max( length(d_ali), length(d_ref)), 2);
imagex(1:length(d_ref),1) = d_ref;
imagex(1:length(d_ali),2) = d_ali;
image( imagex );
colormap( 1 - gray(500) );
hold on
for k = 1:size( corresponding_peaks,2 )
  plot( [1 2], [corresponding_peaks(1,k) corresponding_peaks(2,k)],'r' )
end
for k = 1:length( unassigned_ref )
    plot( 1, peak_ref(unassigned_ref(k)),'mx' )
end
for k = 1:length( unassigned_ali )
  plot( 2, peak_ali(unassigned_ali(k)),'mx' )
end
hold off
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [ corresponding_peaks, unassigned1, unassigned2 ] = align_peaks( d_ref, d_ali, peak_ref, peak_ali );

%for k = 1:min(length( peak_ref ),length( peak_ali ))
%  corresponding_peaks(:,k) = [ peak_ref(k) peak_ali(k) ];
%end

num1 = length( peak_ref );
num2 = length( peak_ali );

DP_scores = zeros( num1+1, num2+1);

prev_align1 = zeros( num1+1, num2+1 );
prev_align2 = zeros( num1+1, num2+1 );

prev_prev_align1 = zeros( num1+1, num2+1 );
prev_prev_align2 = zeros( num1+1, num2+1 );

which_choice = zeros( num1+1, num2+1 );

GAP_PENALTY = 1.0;
INSERTION_PENALTY = GAP_PENALTY;

% initial values.
for i = 1:num1
  DP_scores( i+1, 1 ) = DP_scores( i, 1 ) + GAP_PENALTY;
  which_choice( i+1, 1 ) = 2; % gap
end
for j = 1:num2
  DP_scores( 1, j+1 ) = DP_scores( 1, j ) + INSERTION_PENALTY;
  which_choice( 1, j+1 ) = 3; % insertion
end

% Fill in rest of dynamic programming matrix.
for i = 1: num1
  for j = 1: num2
    
    %%%%%%%%%%%%%%%%
    % Align? The most important metric is distance to previous peak.
    score_aligned = DP_scores( i, j ) ;

    %For debugging
    %score_aligned = score_aligned +  ...
    %	  0.0002* abs( peak_ref(i) -   peak_ali(j) );

    prev_pos1 = prev_align1(i,j);
    prev_pos2 = prev_align2(i,j);

    if (prev_pos1 > 0 & prev_pos2 > 0 )	

      peak_sep1 = peak_ref(i) - peak_ref(prev_pos1);
      peak_sep2 = peak_ali(j) - peak_ali(prev_pos2);

      score_aligned = score_aligned +  ...
      	  0.1 * abs( peak_sep1 - peak_sep2 );
      %score_aligned = score_aligned +  ...
      %	  2.0*abs( log( peak_sep1/peak_sep2 ) );
      
      prev_prev_pos1 = prev_prev_align1(i,j);
      prev_prev_pos2 = prev_prev_align2(i,j);
      
      if (prev_prev_pos1 > 0 & prev_prev_pos2 > 0 & 0 )	
	prev_peak_sep1 = peak_ref(prev_pos1) - peak_ref(prev_prev_pos1);
	prev_peak_sep2 = peak_ali(prev_pos2) - peak_ali(prev_prev_pos2);
	
      	score_aligned = score_aligned +  ...
	    0.2 *abs( log( peak_sep1/prev_peak_sep1 ) - ...
		      log( peak_sep2/prev_peak_sep2 ));
      end

    end
    

    %%%%%%%%%%%%%%%%%%%%
    % Gap
    score_gap = DP_scores( i, j + 1) + GAP_PENALTY;    

    %%%%%%%%%%%%%%%%%%%%
    % Insertion
    score_insertion = DP_scores( i + 1, j) + INSERTION_PENALTY;    

    
    what_choices = [ score_aligned, score_gap,  score_insertion ];
    
    [ best_val, choice ]   = min( what_choices );
    DP_scores(i+1, j+1 ) = best_val;
    which_choice( i+1, j+1 ) = choice;

    if ( length(find( what_choices == best_val )) > 1)
      fprintf( 'Warning! ambiguity! %d %d\n',i,j);
    end
    
    switch choice
     case 1
      prev_align1( i+1, j+1 ) = i;
      prev_align2( i+1, j+1 ) = j;
      prev_prev_align1( i+1, j+1 ) = prev_align1( i, j);
      prev_prev_align2( i+1, j+1 ) = prev_align2( i, j);
     case 2
      prev_align1( i+1, j+1 ) = prev_align1( i, j+1);
      prev_align2( i+1, j+1 ) = prev_align2( i, j+1);
      prev_prev_align1( i+1, j+1 ) = prev_prev_align1( i, j+1);
      prev_prev_align2( i+1, j+1 ) = prev_prev_align2( i, j+1);
     case 3
      prev_align1( i+1, j+1 ) = prev_align1( i+1, j);
      prev_align2( i+1, j+1 ) = prev_align2( i+1, j);
      prev_prev_align1( i+1, j+1 ) = prev_prev_align1( i+1, j);
      prev_prev_align2( i+1, j+1 ) = prev_prev_align2( i+1, j);
    end

  
  end
end
 
%image( 100*DP_scores )
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Back-track.
pos1 = num1+1;
pos2 = num2+1;

corresponding_peaks = [];
unassigned1 = [];
unassigned2 = [];
assigned1 = [];
assigned2 = [];
while( pos1 > 1  & pos2 > 1 )

  choice = which_choice( pos1, pos2 );
  %fprintf( 1, '%d %d %d \n', pos1, pos2, choice);

  switch choice
   case 1
    pos1 = pos1 - 1;
    pos2 = pos2 - 1;

    corresponding_peaks = [ corresponding_peaks; ...
		    peak_ref( pos1) peak_ali( pos2 ) ];
    assigned1 = [assigned1 pos1];
    assigned2 = [assigned2 pos2];
   
   case 2
    % its a gap 
    pos1 = pos1 - 1;
    unassigned1 = [ unassigned1 pos1 ];
   case 3
    % its an insertion
    pos2 = pos2 - 1;
    unassigned2 = [ unassigned2 pos2 ];
  end

end


corresponding_peaks = corresponding_peaks';

return;

