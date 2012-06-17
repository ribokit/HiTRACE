function [ switch_score_combined, data_to_output, data_to_output_err ] = calc_switch_score_RHIJU( inset_from_5prime, inset_from_3prime, ignore_points, sequence, seqpos, area_bsub, darea_bsub, all_area_pred, design_names );
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

clf;
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

    [d_on_norm,  d_on_norm_err  ] = get_data_norm( area_bsub{j}(:, (a-1) * 2 + 2), darea_bsub{j}(:, (a-1) * 2 + 2), goodbins );
    [d_off_norm, d_off_norm_err ] = get_data_norm( area_bsub{j}(:, (a-1) * 2 + 1), darea_bsub{j}(:, (a-1) * 2 + 1), goodbins );

    % One more normalization to try to bring the data close to each other.
    rescale_on_off = mean( d_off_norm( goodbins ) ) / mean( d_on_norm( goodbins ) );
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

  if ( j <= length(which_sets) / 2);
    plot_id = 2 * j - 1;
  else
    plot_id = j - length(which_sets) / 2;
    plot_id = 2 * plot_id;
  end
  subplot(length(which_sets) / 2,2, plot_id);
  
  if j == 1
      title('warning:badQuality');
  end
  num_lanes = size( data_switch, 2);
  nres = size( data_switch, 1 );
  image( seqpos, [1:num_lanes], data_switch'*80 );

  data_to_output{j}     = data_switch;
  data_to_output_err{j} = data_switch_err;

  % draw some grid lines.
  hold on
  for m = 1:(num_lanes+1); plot( [0.5 nres+0.5 ], [m-0.5 m-0.5], 'color', [0.1 0.1 0.1], 'linew', 0.5 ); end;
  for m = 1:(nres+1); plot( [m-0.5 m-0.5], [0.5 num_lanes+0.5 ], 'color', [0.1 0.1 0.1], 'linew', 0.5 ); end;
  % show sequence.
  for m = 1:nres; text( m, 0.5, sequence{j}(m),'verticalalign','bottom','horizontalalign','center' ); end;

  % show annotations of where switch should occur
  for n = 1:nres
    for m = 1:num_lanes
      if ( all_area_pred{j}(n,m) > 0)
	plot( seqpos(n), m, 'ro'); 
      end;
    end
  end

  % show annotations of whether switch occurred.
  % gray square -- not evaluated (excluded position)
  str_on  = all_area_pred{j}(:, 2);
  str_off = all_area_pred{j}(:, 1);
  for n = 1:nres
    if  isempty( find(n==goodbins) )
      plot( seqpos(n),  num_lanes+1, 's','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'clipping','off' );
    else
      % black square -- not a switch position
      if ( str_on(n) == str_off(n) )
	plot( seqpos(n),  num_lanes+1, 's','color','k','markerfacecolor','k','clipping','off' );
      else
	if (s_combine(n) == 0 )
	  % red x -- no switch where there should be one
	  plot( seqpos(n),  num_lanes+1, 'x','color','r','linew',2,'clipping','off' );	  
	else
	  % green circle -- OK. Strength of color indicates how strong the switch is.
	  colorcode = 1 - ( 1 - [0, 0.5, 0] ) * s_combine(n);
	  plot( seqpos(n),  num_lanes+1, 'o','color','k','markerfacecolor',colorcode,'clipping','off' );	  
	end
      end
    end      
  end
  ylim([-0.5 num_lanes+0.5])
  xlim([-0.5 nres+0.5])
  if j == 1
      text( nres+1, 0.5, ['Sequence ',num2str(j),'\newline', design_names{j}, num2str(j),'\newline', 'warning:badQuality'],'fontsize',8,'fontwe','bold','verticalalign','top' );
  else
      text( nres+1, 0.5, ['Sequence ',num2str(j),'\newline', design_names{j}],'fontsize',8,'fontwe','bold','verticalalign','top' );
  end
  text( -0.5, 0.5, ['Switch\newlinescore: \newline',...
		  num2str(sum( s_combine( find(switch_bin_SHAPE) ) ),'%3.1f'),'/',...
		  num2str(sum(switch_bin_SHAPE),'%d'), '\newline',...
		  num2str(switch_score_combined(j),'%8.1f')] ,'fontsize',8,'fontwe','bold','horizontalalign','right','verticalalign','top' );
  hold off
  axis off
end

colormap( 1 - gray(100) );
set(gcf, 'PaperPositionMode','auto','color','white');



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
