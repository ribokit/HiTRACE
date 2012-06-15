%% description
% parameter
%     ADAPTIVE = 0: fixed threshold
%     ADAPTIVE = 1: adaptive threshold

% output
%     switch_score: paired-unpaired comparison based scoring
%     score2      : alternative scoring
%     FIXED_SCORE : for checking correlation with traditional ETERNA score
%     score_edit(i, j , k)    : score based edit distance and eterna score
%                  (i-th sequence, j-th weight composite, k - SHAPE/DMS)
%     score_hamming(i, j , k) : score based hamming distance and eterna score
%                  (i-th sequence, j-th weight composite, k - SHAPE/DMS)

%% initialization
addpath(genpath('../'));

alt_structure = '.....((((......((((....)))).....))))....................................';
structure = '..........................(((((((............)))))))....................';

load 'data.mat';
load 'area_peak.mat';

alt_lane = [2 4];

offset = 0;
dist = 20;
which_sets = 1:12;

ADAPTIVE = 0; 

%START = 6;
%END = 6;

% I included an additional point near the 3' end to get a better read on the switch score. -- Rhiju
START = 5;
END = 5;

% These nucleotides are in the FMN aptamer region -- they go from paired to 'unpaired' but
% really can get protected as they structure around the FMN molecule. We should ignore the data 
% there in computing the switch score. -- Rhiju
ignore_points = [10:15 28:32];

ETERNA_score = [];

%% prepare area_pred matrix from structure
for j = which_sets  
   seqpos = length(sequence{j})-dist - [1:(length(sequence{j})-dist)] + 1 + offset;
   [ marks{j}, all_area_pred{j}, mutpos{j} ] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence{j}, offset , seqpos, data_types );
end
    
for j = which_sets  
  seqpos = length(sequence{j})-dist - [1:(length(sequence{j})-dist)] + 1 + offset;
  [ alt_marks{j}, alt_area_pred{j}, alt_mutpos{j} ] = get_predicted_marks_SHAPE_DMS_CMCT( alt_structure, sequence{j}, offset , seqpos, data_types );
  all_area_pred{j}(:,alt_lane) = alt_area_pred{j}(:,alt_lane);
end

DO_ETERNA_SCORE_CALC = 0;

NRES = size( area_bsub{1}, 1); % NRES = 52

if DO_ETERNA_SCORE_CALC
  %% EteRNA Score calculation
  for j = 1:12
    for a = find(strcmp(data_types, 'SHAPE') | strcmp(data_types, 'DMS'));
      
      nres = size( area_bsub{j}, 1);
      
      goodbins = [(nres-END):-1:START+1]; % have to go backwards. Silly convention switch.
      data_cols = [a];
      
      % average over appropriate columns...
      data = mean(area_bsub{j}( goodbins, data_cols ), 2)';
      
      % normalize.
      [data_norm, scalefactor] = SHAPE_normalize( data );
      area_bsub_norm{j}(:,data_cols)  = area_bsub{j}(:,data_cols)/scalefactor;
      darea_bsub_norm{j}(:,data_cols) = darea_bsub{j}(:,data_cols)/scalefactor;
      
      % To compute EteRNA score, need predicted paired/unpaired for 'perfect' design.
      % this was computed above to aid in sequence annotation, background subtraction, etc.
      pred = all_area_pred{j}( goodbins, data_cols ); 
      
      %  I optimized this scaling to match the previous EteRNA score with adaptive thresholds -- Rhiju
      scalefactor = 2.0;

      others = 0;
      for k = 1:40
	if(pred(k) == 1)
          if( data_norm(k) > 0.125/scalefactor)
	    others = others + 1;
          end
	else
          if( data_norm(k) < 0.5/scalefactor)
	    others = others + 1;
          end
	end
      end
      
      subplot(2,1,1);
      plot( data_norm*scalefactor, 'rx-' ); hold on
      plot( pred,'k','linew',2 ); hold off
      ylim([-0.5 2]);
      %pause;
      
      % ETERNA score from fixed threshold.
      FIXED_SCORE(j,a) = others / (NRES - START - END) * 100; 
      
      % traditional ETERNA score.
      [min_SHAPE{j,a}, max_SHAPE{j,a}, threshold_SHAPE{j,a}, ETERNA_score(j,a), d_bin{j}(:,a)] = determine_thresholds_binarization_and_ETERNA_score( data_norm, pred );
      subplot(2,1,1); ylim([-0.5 2]);
      
      %pause;
    end
  end
end

%% pre-calculation for switch score from base-pair comparison -- option 1 (adaptive/fixed threshold)
% case 1: paired (off) - paired (on)
% case 2: paired (off) - unpaired (on)
% case 3: unpaired (off) - paired (on)
% case 4: unpaired (off) - unpaired (on)

% remove nucleotide from FMN-binding region for consideration in switch score.
goodbins = (START+1):(NRES-END);
for m = 1:length( ignore_points)
  goodbins = setdiff( goodbins, find(seqpos == ignore_points(m)) );
end

for i = 1:12
  fprintf( 'Sequence %d\n',i);
  clf;

  % this keeps track of whether a switch occurred at the right place or not.
  % according to SHAPE (a=1) or DMS (a=2).
  s = zeros( NRES, 2 );
  for a = 1:2

    subplot(2,1,a);
    str_on  = all_area_pred{i}(:, (a-1) * 2 + 2);
    str_off = all_area_pred{i}(:, (a-1) * 2 + 1);
    
    % Max/Min to avoid 'noisy' points. -- Rhiju
    % Also note that we use area_peak instead of area_bsub -- the backsub adds noise to the difference comparison. -- Rhiju
    d_on_norm  = min( max( SHAPE_normalize(area_peak{i}(:, (a-1) * 2 + 2)), 0), 2.0);
    d_off_norm = min( max( SHAPE_normalize(area_peak{i}(:, (a-1) * 2 + 1)), 0), 2.0 );
    
    % One more normalization to try to bring the data close to each other.
    d_on_norm = d_on_norm * mean( d_off_norm( goodbins ) ) / mean( d_on_norm( goodbins ) );
        
    if(ADAPTIVE)
      % Adaptive threshold calculation
      threshold = mean([threshold_SHAPE{i,(a-1) * 2 + 1} threshold_SHAPE{i,(a-1) * 2 + 2}]);
      min_th = mean([min_SHAPE{i,(a-1) * 2 + 1} min_SHAPE{i,(a-1) * 2 + 2}]);
      
      threshold_is_a_change = threshold;
      threshold_not_a_change = threshol;
    else
      % fixed threshold -- like above, give the design the 'benefit of the doubt' by using different thresholds when looking for a change  -- Rhiju
      %  compared to when looking for fixed.
      threshold_is_a_change = 0.2;
      threshold_not_a_change = 0.5;
      min_th = 0;
    end
    
    s_tot = 0; n_tot = 0;	

    % Keep track of which data points really switch
    % In the end, I just focused on the places that *should* switch. -- Rhiju
    switch_bin = 0;    
    for j = goodbins
      switch_bin(j) = 1;     
      if(str_on(j) == str_off(j)) % case 1,4
	switch_bin(j) = 0;
	if(abs(d_on_norm(j)-d_off_norm(j)) < threshold_not_a_change )
	  s(j,a) = 1;
	else
	  s(j,a) = 0;
	end
      elseif(str_on(j) > str_off(j)) % case 2
	dif = d_on_norm(j) - d_off_norm(j);
	if(dif > threshold_is_a_change)
	  s(j,a) = 1;
	elseif ( dif < min_th)
	  s(j,a) = 0;
	else
	  s(j,a) = dif/threshold_is_a_change;
	end
      else % case 3
	dif = d_on_norm(j) - d_off_norm(j);
	if(dif < -threshold_is_a_change)
	  s(j,a) = 1;
	elseif ( dif > min_th)
	  s(j,a) = 0;
	else
	  s(j,a) = -dif/threshold_is_a_change;
	end
      end
    end

    % only save switch calls at switch positions!
    not_switch_points = find( ~switch_bin );
    s(not_switch_points,a) = 0;
    
    if (a==1)
      switch_bin_SHAPE = switch_bin; % placed where we are evaluating whether switch occurred.
    end
    
    switch_score(i,a) = 100 * sum( s( find(switch_bin),a ) )/ sum( switch_bin ) ;
    fprintf( 1, 'Switch score %d: %8.1f\n ', a, switch_score(i,a) );
    
    plot( seqpos, str_on,  'k', 'linew', 2 ); hold on;
    plot( seqpos, str_off, 'b', 'linew', 2 ); 
    plot( seqpos, d_on_norm, 'kx-' );
    plot( seqpos, d_off_norm, 'bx-' );
    
    plot_points = find( switch_bin );
    plot( seqpos(plot_points), s( plot_points,a ), 'go','markerfacecolor','g' );
    hold off
    ylim([-0.5 3]);
    xlim([0 50])
    set(gca,'xtick',[0:5:50],'xgrid','on');
    
  end

  %plot_points = find( switch_bin_SHAPE );
  %plot( plot_points, s( plot_points,1 ), 'ko','markerfacec','k' );hold on
  %plot( plot_points, s( plot_points,2 ), 'bx','markerfacec','b' );hold off

  % Let either SHAPE or DMS give evidence of switch.  
  s_combine = max(s,[],2); 
  switch_score_combined(i) = 100 * sum( s_combine( find(switch_bin_SHAPE) ) )/sum( switch_bin_SHAPE );
  fprintf( 1, 'Switch score SHAPE/DMS: %8.1f\n ', switch_score_combined(i) );

  fprintf( '\n');
  pause
end

switch_score_combined


% I didn't optimize the following -- Rhiju

%% calculation switch score from EteRNA score and exp with hamming distance -- option 2
f_shape = ETERNA_score(:,1) / 100;
g_shape = ETERNA_score(:,2) / 100;
f_dms = ETERNA_score(:,3) / 100;
g_dms = ETERNA_score(:,4) / 100;
for j = 1:12
  for a = 1:2
        ideal_on_str = [];
        ideal_off_str = [];

        actual_on_str = [];
        actual_off_str = [];

        for k = 1:(length(all_area_pred{j})-START-END)
            ideal_on_str = [ideal_on_str num2str(all_area_pred{j}(START + k,(a -1) * 2 + 2))];
            ideal_off_str = [ideal_off_str num2str(all_area_pred{j}(START + k,(a -1) * 2 + 1))];
            
            actual_on_str = [actual_on_str num2str(d_bin{j}(k,(a -1) * 2 + 2))];
            actual_off_str = [actual_off_str num2str(d_bin{j}(k,(a -1) * 2 + 1))];
        end
        f = ETERNA_score(j, (a -1) * 2 + 1) / 100;
        g = ETERNA_score(j, (a -1) * 2 + 2) / 100;
        
        d_f = 1 - f;
        d_g = 1 - g;
        
        ideal_on_str = double(ideal_on_str) - '0';
        ideal_off_str = double(ideal_off_str) - '0';
        actual_on_str = double(actual_on_str) - '0';
        actual_off_str = double(actual_off_str) - '0';
        
        d = pdist([ideal_on_str; ideal_off_str],'hamming');
        e = pdist([actual_on_str; actual_off_str],'hamming');
        x = pdist([ideal_off_str; actual_on_str],'hamming');
        y = pdist([ideal_on_str; actual_off_str],'hamming');

        w_mat = [];
        
        
        % method 2
        for k = 0:63       
            weight = dec2bin(k,6);
            
            for l = 1:6
                w(l) = str2double(weight(l));
            end
            score_hamming(j, k+1, a) = geomean([f g]) * exp(-( w(1) * ((d-e) / length(ideal_on_str))^2 + w(2) * ((x - y) / length(ideal_on_str))^2 + ...
                w(3) * ((d - y) / length(ideal_on_str))^2 + w(4) * ((e - y) / length(ideal_on_str))^2 + ...
                w(5) * ((d - x) / length(ideal_on_str))^2 + w(6) * ((e - x) / length(ideal_on_str))^2 )) * 100;
            
            w_mat = [w_mat w'];
        end
        
        % method 1
        score2(j,a) = geomean([f g switch_score(j,a)/100]) * 100;
    end
end

score2
%score_hamming

%% calculation switch score from EteRNA score and exp with edit distance --option 3
for j = 1:12
    for a = 1:2
        ideal_on_str = [];
        ideal_off_str = [];

        actual_on_str = [];
        actual_off_str = [];

        for k = 1:(length(all_area_pred{j})-START-END)
            ideal_on_str = [ideal_on_str num2str(all_area_pred{j}(START + k,(a -1) * 2 + 2))];
            ideal_off_str = [ideal_off_str num2str(all_area_pred{j}(START + k,(a -1) * 2 + 1))];
            
            actual_on_str = [actual_on_str num2str(d_bin{j}(k,(a -1) * 2 + 2))];
            actual_off_str = [actual_off_str num2str(d_bin{j}(k,(a -1) * 2 + 1))];
        end
        f = ETERNA_score(j, (a -1) * 2 + 1) / 100;
        g = ETERNA_score(j, (a -1) * 2 + 2) / 100;
        
        d_f = 1 - f;
        d_g = 1 - g;
        
        d = EditDist(ideal_on_str, ideal_off_str);
        e = EditDist(actual_on_str, actual_off_str);
        
        x = EditDist(ideal_off_str, actual_on_str);
        y = EditDist(ideal_on_str, actual_off_str);

        w_mat = [];
        
        for k = 0:63       
            weight = dec2bin(k,6); % perhaps, only 0, 32, 37, 38, 41, 42, 44 and 48 are meaningful
            
            for l = 1:6
                w(l) = str2double(weight(l));
            end
            score_edit(j, k+1, a) = geomean([f g]) * exp(-( w(1) * ((d-e) / length(ideal_on_str))^2 + w(2) * ((x - y) / length(ideal_on_str))^2 + ...
                w(3) * ((d - y) / length(ideal_on_str))^2 + w(4) * ((e - y) / length(ideal_on_str))^2 + ...
                w(5) * ((d - x) / length(ideal_on_str))^2 + w(6) * ((e - x) / length(ideal_on_str))^2 )) * 100;
            
            w_mat = [w_mat w'];
        end
    end
end

%score_edit
