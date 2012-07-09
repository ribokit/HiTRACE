function [ETERNA_score, min_SHAPE, max_SHAPE, threshold_SHAPE] = calc_eterna_score_GUI( inset_from_5prime, inset_from_3prime, data_types, data_to_output, sequence, seqpos, area_bsub, all_area_pred, design_names );

min_SHAPE = {};
max_SHAPE = {};
threshold_SHAPE = {};
ETERNA_score = {};

COMPARE_TO_OLD_SCORE = 0;

START = inset_from_3prime; % where to start data for eterna input
END   = inset_from_5prime;   % where to end data for eterna input.

which_sets = 1:length( area_bsub );
num_sets = length( which_sets );
for j = which_sets


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %create image for visual feedback
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  num_lanes  = size( data_to_output{j}, 2 );
  nres       = size( data_to_output{j}, 1);

%   data_image = zeros( nres, 3*num_lanes);
%   for n = 1:num_lanes
%     data_image(:,3*(n-1)+1) = data_to_output{j}(:,n);
%     data_image(:,3*(n-1)+2) = data_to_output{j}(:,n);
%   end
%   subplot(length(which_sets),1,j);
%   image( seqpos, [1 : 3*num_lanes], 80*data_image' );
%   % draw some grid lines.
%   hold on

%   for q = 1:(num_lanes); 
%     plot( [0.5 nres+0.5 ], (3*(q-1)+0.5)*[1 1], 'color', [0.1 0.1 0.1], 'linew', 0.5 );
%     plot( [0.5 nres+0.5 ], (3*(q-1)+2.5)*[1 1], 'color', [0.1 0.1 0.1], 'linew', 0.5 );
%     for m = 1:(nres+1); plot( [m-0.5 m-0.5], 3*(q-1)+[0.5 2.5], 'color', [0.1 0.1 0.1], 'linew', 0.5 ); end;
%   end
  % show sequence.
%   for m = 1:nres; text( m, 0.5, sequence{j}(m),'verticalalign','bottom','horizontalalign','center' ); end;
%   ylim([-0.5 3*num_lanes+0.5]);
%   xlim([-0.5 nres+0.5]);
%   axis off
%   text( nres+1, 0.5, ['Sequence ',num2str(j),'\newline', design_names{j}],'fontsize',8,'fontwe','bold','verticalalign','top' );
%   hold on
%   
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculate eterna score.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n = 1:num_lanes
        
    goodbins = [(nres-END):-1:START]; % have to go backwards. Silly convention switch
				
    % In case of DMS, only 'A' and 'C' are handled
    if( strcmp(data_types{n},'DMS') )      
      seq = sequence{j}( seqpos );
      good_pos = find(seq == 'A' | seq == 'C');
      goodbins = intersect(goodbins, good_pos);
    end
    
    % average over appropriate columns...
    data = data_to_output{j}(:,n); % don't renormalize
    pred =  all_area_pred{j}(:,n); 

    % fixed thresholds!
    min_SHAPE_fixed = 0.0;
    max_SHAPE_fixed = 1.0;
    threshold_SHAPE_fixed = 0.50;

    correct_hit = zeros( 1, nres );
    for k = goodbins
      if ( pred(k) == 1 )
	if ( data(k) > (0.25*threshold_SHAPE_fixed + 0.75*min_SHAPE_fixed ) )
	  correct_hit(k) = 1;
	end
      else
	if ( data(k) < threshold_SHAPE_fixed)
	  correct_hit(k) = 1;
	end
      end
    end

    %subplot(2,1,1);
    %plot( data(goodbins), 'rx-' ); hold on
    %plot( pred(goodbins),'k','linew',2 ); hold off
    %ylim([-0.5 2]);

    eterna_fixed_score = sum( correct_hit( goodbins ) )/ length( goodbins ) * 100

    min_SHAPE{j,n} = min_SHAPE_fixed;
    max_SHAPE{j,n} = max_SHAPE_fixed;
    threshold_SHAPE{j,n} = threshold_SHAPE_fixed;
    ETERNA_score{j,n} = eterna_fixed_score;
%     text( -0.5,  3*(n-1)+1.5, num2str(eterna_fixed_score, '%8.1f' ) ,'clipping','off','fontweight','bold','fontsize',8,'horizontalalign','right' );
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % show annotations of what was scored as 'correct'
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for q = 1:nres
%       if ( pred(q) > 0)
% 	plot( seqpos(q), 3*(n-1)+1.5, 'ro'); 
%       end;
%     end
% 
%     %
%     for q = 1:nres
%       if  isempty( find(q==goodbins) )
% 	plot( seqpos(q),  3*(n-1)+3, 's','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'clipping','off' );
%       else
% 	if (correct_hit(q) == 0 )
% 	  % red x -- no switch where there should be one
% 	  plot( seqpos(q), 3*(n-1)+3, 'x','color','r','linew',2,'clipping','off' );	  
% 	else
% 	  % green circle -- OK. Strength of color indicates how strong the switch is.
% 	  colorcode =[0 0.5 0];
% 	  plot( seqpos(q), 3*(n-1)+3, 'o','color','k','markerfacecolor',colorcode,'clipping','off' );	  
% 	end
%       end
%     end      
%     
%     if COMPARE_TO_OLD_SCORE
%       % To compute EteRNA score, need predicted paired/unpaired for 'perfect' design.
%       % this was computed above to aid in sequence annotation, background subtraction, etc.
%       data_norm = data( goodbins)';
%       pred      = pred( goodbins);
%       [min_SHAPE{j,n}, max_SHAPE{j,n}, threshold_SHAPE{j,n}, ETERNA_score{j,n} ] = determine_thresholds_and_ETERNA_score( data_norm, pred );
%       pause;
%     end
%     
  end
%   hold off

end