function [ switch_score_combined, data_to_output, data_to_output_err ] = calc_switch_score_GUI( inset_from_5prime, inset_from_3prime, ignore_points, sequence, seqpos, area_bsub, darea_bsub, all_area_pred, design_names );
% This script is in charge of (1) normalizing the data, (2) making calls on which nucleotide 'switch' appropriately, (3) giving back 
%  graphical display of this scoring, (4) compiling the overall 'switch score'.
%
% The data are normalized so that their mean is 0.5. (For DMS, the mean at A & C is 0.5.)
%

%threshold_is_a_change = 0.2;
threshold_is_a_change = 0.1; % this parameter controls the cutoff for what is a switch
threshold_not_a_change = 0.5; % this is not actually in use. 

START = inset_from_3prime; % where to start data for eterna input
END   = inset_from_5prime;   % where to end data for eterna input.

% clf;
which_sets = 1:length( area_bsub );
num_sets = length( which_sets );

for j = which_sets

  nres = size( area_bsub{j}, 1 );
  goodbins = START:(nres-END);
  for m = 1:length( ignore_points)
      goodbins = setdiff( goodbins, find(seqpos == ignore_points(m)) );
  end

  fprintf( 'Sequence %d\n',j);

  data_switch = [];  
  
  s = zeros( nres, 2 );
  for a = 1:2
    str_on = all_area_pred{j}(:, (a-1) * 2 + 2);
    str_off = all_area_pred{j}(:, (a-1) * 2 + 1);
    % Max/Min to avoid 'noisy' points. -- Rhiju
    % Also note that we use area_peak instead of area_bsub -- the backsub adds noise to the difference comparison. -- Rhiju

    %d_on_norm  = min( max( SHAPE_normalize(area_peak{j}(:, (a-1) * 2 + 2)), 0), 2.0);
    %d_off_norm = min( max( SHAPE_normalize(area_peak{j}(:, (a-1) * 2 + 1)), 0), 2.0 );

    %d_on_norm  = min( max( SHAPE_normalize(area_bsub{j}(:, (a-1) * 2 + 2)), 0), 2.0);
    %d_off_norm = min( max( SHAPE_normalize(area_bsub{j}(:, (a-1) * 2 + 1)), 0), 2.0 );

    [d_on_norm,  d_on_norm_err  ] = get_std_norm( area_bsub{j}(:, (a-1) * 2 + 2), darea_bsub{j}(:, (a-1) * 2 + 2), goodbins,all_area_pred{j}(:, (a-1) * 2 + 2));
    [d_off_norm, d_off_norm_err ] = get_data_norm( area_bsub{j}(:, (a-1) * 2 + 1), darea_bsub{j}(:, (a-1) * 2 + 1), goodbins );

    % normalize filter for removing outliers
    
    ordered_on = sort(d_on_norm);
    ordered_off = sort(d_off_norm);
  
    x25_on = ordered_on(round(0.25 * length(ordered_on)));
    x75_on = ordered_on(round(0.75 * length(ordered_on)));
    
    x25_off = ordered_off(round(0.25 * length(ordered_off)));
    x75_off = ordered_off(round(0.75 * length(ordered_off)));
      
    target_idx_on = find(d_on_norm <= 1.5 * (x75_on -x25_on) + x75_on);
    target_idx_off = find(d_off_norm <= 1.5 * (x75_off -x25_off) + x75_off);
    
    target_idx_on = intersect(target_idx_on, goodbins);
    target_idx_off = intersect(target_idx_off, goodbins);
    
    % One more normalization to try to bring the data close to each other.
    rescale_on_off = mean( d_off_norm( target_idx_off ) ) / mean( d_on_norm( target_idx_on ) );
    d_on_norm     = d_on_norm * rescale_on_off;
    d_on_norm_err = d_on_norm_err * rescale_on_off;

    data_switch(:,  (a-1)*2+2  ) = d_on_norm;
    data_switch(:,  (a-1)*2+1  ) = d_off_norm;

    data_switch_err(:,  (a-1)*2+2  ) = d_on_norm_err;
    data_switch_err(:,  (a-1)*2+1  ) = d_off_norm_err;
    
    min_th = 0;

    s_tot = 0; n_tot = 0;	

    % Keep track of which data points really switch
    % In the end, I just focused on the places that *should* switch. -- Rhiju
    switch_bin = 0;    
    for k = goodbins
      switch_bin(k) = 1;     

      if(str_on(k) == str_off(k)) % case 1,4 -- not actually used. See below.
	switch_bin(k) = 0;
	if(abs(d_on_norm(k)-d_off_norm(k)) < threshold_not_a_change )
	  s(k,a) = 1;
	else
	  s(k,a) = 0;
	end
      
      elseif(str_on(k) > str_off(k)) % case 2
	dif = d_on_norm(k) - d_off_norm(k);
	if(dif > threshold_is_a_change)
	  s(k,a) = 1;
	elseif ( dif < min_th)
	  s(k,a) = 0;
	else
	  s(k,a) = dif/threshold_is_a_change;
	end
      
      else % case 3
	dif = d_on_norm(k) - d_off_norm(k);
	if(dif < -threshold_is_a_change)
	  s(k,a) = 1;
	elseif ( dif > min_th)
	  s(k,a) = 0;
	else
	  s(k,a) = -dif/threshold_is_a_change;
	end
      end
    end
    
    % only save switch calls at switch positions!
    not_switch_points = find( ~switch_bin );
    s(not_switch_points,a) = 0;

    if (a==1)
      switch_bin_SHAPE = switch_bin; % placed where we are evaluating whether switch occurred.
    end

    switch_score(j,a) = 100 * sum( s( find(switch_bin),a ) )/ sum( switch_bin ) ;
    fprintf( 1, 'Switch score %d: %8.1f\n ', a, switch_score(j,a) );
  end

  s_combine = max(s,[],2); 
  switch_score_combined(j) = 100 * sum( s_combine( find(switch_bin_SHAPE) ) )/sum( switch_bin_SHAPE );
  fprintf( 1, 'Switch score SHAPE/DMS: %8.1f\n ', switch_score_combined(j) );
  %pause;

%   subplot(num_sets,1,j);
  num_lanes = size( data_switch, 2);
  nres = size( data_switch, 1 );
%   image( seqpos, [1:num_lanes], data_switch'*80 );
% 
  data_to_output{j}     = data_switch;
  data_to_output_err{j} = data_switch_err;
% 
%   % draw some grid lines.
%   hold on
%   for m = 1:(num_lanes+1); plot( [0.5 nres+0.5 ], [m-0.5 m-0.5], 'color', [0.1 0.1 0.1], 'linew', 0.5 ); end;
%   for m = 1:(nres+1); plot( [m-0.5 m-0.5], [0.5 num_lanes+0.5 ], 'color', [0.1 0.1 0.1], 'linew', 0.5 ); end;
%   % show sequence.
%   for m = 1:nres; text( m, 0.5, sequence{j}(m),'verticalalign','bottom','horizontalalign','center' ); end;
% 
%   % show annotations of where switch should occur
%   for n = 1:nres
%     for m = 1:num_lanes
%       if ( all_area_pred{j}(n,m) > 0)
% 	plot( seqpos(n), m, 'ro'); 
%       end;
%     end
%   end
% 
%   % show annotations of whether switch occurred.
%   % gray square -- not evaluated (excluded position)
%   str_on  = all_area_pred{j}(:, 2);
%   str_off = all_area_pred{j}(:, 1);
%   for n = 1:nres
%     if  isempty( find(n==goodbins) )
%       plot( seqpos(n),  num_lanes+1, 's','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'clipping','off' );
%     else
%       % black square -- not a switch position
%       if ( str_on(n) == str_off(n) )
% 	plot( seqpos(n),  num_lanes+1, 's','color','k','markerfacecolor','k','clipping','off' );
%       else
% 	if (s_combine(n) == 0 )
% 	  % red x -- no switch where there should be one
% 	  plot( seqpos(n),  num_lanes+1, 'x','color','r','linew',2,'clipping','off' );	  
% 	else
% 	  % green circle -- OK. Strength of color indicates how strong the switch is.
% 	  colorcode = 1 - ( 1 - [0, 0.5, 0] ) * s_combine(n);
% 	  plot( seqpos(n),  num_lanes+1, 'o','color','k','markerfacecolor',colorcode,'clipping','off' );	  
% 	end
%       end
%     end      
%   end
%   ylim([-0.5 num_lanes+0.5])
%   xlim([-0.5 nres+0.5])
%   text( nres+1, 0.5, ['Sequence ',num2str(j),'\newline', put_newlines_in_name(design_names{j}) ],'fontsize',8,'fontwe','bold','verticalalign','top' );
%   text( -0.5, 0.5, ['Switch\newlinescore: \newline',...
% 		  num2str(sum( s_combine( find(switch_bin_SHAPE) ) ),'%3.1f'),'/',...
% 		  num2str(sum(switch_bin_SHAPE),'%d'), '\newline',...
% 		  num2str(switch_score_combined(j),'%8.1f')] ,'fontsize',8,'fontwe','bold','horizontalalign','right','verticalalign','top' );
%   hold off
%   axis off
end
% 
% colormap( 1 - gray(100) );
% set(gcf, 'PaperPositionMode','auto','color','white');
% 


function [d_norm, d_norm_err] = get_std_norm(d, d_err, goodbins, area_pred);
d_norm =d;
d_norm_err = d_err;

data = d(goodbins)';
pred = area_pred(goodbins);

n = length( goodbins );

% penalties for:
% scale-factor  baseline    dev>0       dev<0      a little unpenalized slop    
f  = [   0,           0,    ones(1,n),  ones(1,n), zeros( 1, n )           ];

% Force all coefficients to be positive.
% Force coefficient of data to be at least 0.5...
LB(1) = 0.60 / mean( max(data,0.0) );
%LB(1) = 0.70 / mean( max(data,0.0) );

data_sort=sort(data);
data_range = abs( data_sort( floor(n*0.1)+1) - data_sort( floor(n*0.9)+1 ) );
%LB(1) = 0.3 / data_range
%UB(1) = 1.3 / data_range;
%LB(1) = 0.0;
UB(1) = inf;

% baseline
LB(2) = -0.1;
UB(2) =  0.1;

% 'slack' variable  -- positive deviations from data.
LB( 2 + [1:n] ) = zeros(1,n);
UB( 2 + [1:n] ) = inf * ones(1,n);

% 'slack' variable  -- negative deviations from data.
LB( 2 + n + [1:n] ) = zeros(1,n);
UB( 2 + n + [1:n] ) = inf * ones(1,n);

% 'slop' variables -- some range is allowed.
% for unpaired region, allow deviations up to 2-fold above 'max', and down to halfway point.
pred_high = find(  pred );

%LB( 2 + 2*n + pred_high) = -2.0;
%UB( 2 + 2*n + pred_high) =  0.5;

LB( 2 + 2*n + pred_high) =  0.0;
UB( 2 + 2*n + pred_high) =  0.0;

% for protected region, ask that we stay close to zero.
pred_low  = find( ~pred);
LB( 2 + 2*n + pred_low) = -0.0;
UB( 2 + 2*n + pred_low) = +0.0;

% how to transform from variables to prediction
Aeq = [ data', ones(n,1), eye(n,n), -eye(n,n), eye(n,n) ]; 

% what we're trying to match -- note that we're looking for equality
beq = pred;

options = optimset('Display','off');
% linear program asking for equality. (The [], [] would define greater-than, less-than constraints).
params = linprog( f, [], [], Aeq, beq, LB, UB,[],options);

scale_factor =  params(1);
baseline     =  params(2);

d_norm = d_norm * scale_factor + baseline;
d_norm_err = d_norm_err * scale_factor + baseline;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [d_norm, d_norm_err] = get_data_norm( d, d_err, goodbins );

[dummy, scalefactor] = quick_norm(d, goodbins);
d_norm     = 0.5 * d * scalefactor;
d_norm_err = 0.5 * d_err * scalefactor; 
d_norm  = min( max( d_norm, 0), 2.0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name_out =  put_newlines_in_name( name )

buffsize = 15;
name_out = '';
for i = 1:length( name )
  name_out = [name_out, name(i) ];
  if ( mod(i,buffsize) == 0 )
    name_out = [name_out, '\newline'];
  end
end
