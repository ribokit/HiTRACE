clear all;

%% script path (HiTRACE_HOME/Scripts)
addpath(genpath(strcat('../')));

%% pick your own workspace filename;
uiload;

which_sets = 1:length(sequence);

%% generating area_pred
for j = which_sets  
  seqpos = length(sequence{j})-20 - [1:(length(sequence{j})-20)] + 1;
  [ marks{j}, all_area_pred{j}, mutpos{j} ] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence{j}, 0 , seqpos, data_types );
end

% area_pred for switched structure
if(~isempty(alt_structure))
    alt_lane = [2 4];
    for j = which_sets  
      seqpos = length(sequence{j})-20 - [1:(length(sequence{j})-20)] + 1;
      [ alt_marks{j}, alt_area_pred{j}, alt_mutpos{j} ] = get_predicted_marks_SHAPE_DMS_CMCT( alt_structure, sequence{j}, 0 , seqpos, data_types );
      all_area_pred{j}(:,alt_lane) = alt_area_pred{j}(:,alt_lane);
    end
end

% peak fitting
for j = which_sets
    [area_peak{j}, prof_fit{j}] = do_the_fit_fast( d_bsub{j}, xsel{j}', 0.0, 0);
end

backgd_sub_col = find(strcmp(data_types, 'nomod'));
penalize_negative_weight = 2.0;

for j = which_sets
  [ area_bsub{j}, darea_bsub{j}] = overmod_and_background_correct_logL( area_peak{j}, backgd_sub_col, [4:size(area_peak{j},1)-4], all_area_pred{j});
end


% Average data over Ann/DC replicates, and output to screen in nice format that can be copy/pasted into EteRNA.
% Also, automated EteRNA scoring.
START = 6; % where to start data for eterna input
END   = 6;   % where to end data for eterna input.


min_SHAPE = {}; max_SHAPE={};threshold_SHAPE={};ETERNA_score={};
area_bsub_norm  = area_bsub;

for j = which_sets
  for a = find(strcmp(data_types, 'SHAPE') | strcmp(data_types, 'DMS'));

  nres = size( area_bsub{j}, 1);

  goodbins = [(nres-END):-1:START+1]; % have to go backwards. Silly convention switch.
  data_cols = [a];

  % In case of DMS, only 'A' and 'C' are handled
%               if(strcmp(data_types{a},'DMS'))
%                   seq = sequence{i}(end-dist:-1:1);
%                   good_pos = find(seq == 'A' | seq == 'C');
%                   goodbins = intersect(goodbins, good_pos);
%               end

  % average over appropriate columns...
  data = mean(area_bsub{j}( goodbins, data_cols ), 2)';

  % normalize.
  [data_norm, scalefactor] = SHAPE_normalize( data );
  area_bsub_norm{j}(:,data_cols)  = area_bsub{j}(:,data_cols)/scalefactor;
  darea_bsub_norm{j}(:,data_cols) = darea_bsub{j}(:,data_cols)/scalefactor;

  % To compute EteRNA score, need predicted paired/unpaired for 'perfect' design.
  % this was computed above to aid in sequence annotation, background subtraction, etc.
  pred = all_area_pred{j}( goodbins, data_cols ); 
  [min_SHAPE{j,a}, max_SHAPE{j,a}, threshold_SHAPE{j,a}, ETERNA_score{j,a} ] = determine_thresholds_and_ETERNA_score( data_norm, pred );
  end

  ignore_points = [10:15 28:32];

  goodbins = (START+1):(nres-END);
  for m = 1:length( ignore_points)
      goodbins = setdiff( goodbins, find(seqpos == ignore_points(m)) );
  end

  fprintf( 'Sequence %d\n',j);

  s = zeros( nres, 2 );
  for a = 1:2
    str_on = all_area_pred{j}(:, (a-1) * 2 + 2);
    str_off = all_area_pred{j}(:, (a-1) * 2 + 1);
    % Max/Min to avoid 'noisy' points. -- Rhiju
    % Also note that we use area_peak instead of area_bsub -- the backsub adds noise to the difference comparison. -- Rhiju
    d_on_norm  = min( max( SHAPE_normalize(area_peak{j}(:, (a-1) * 2 + 2)), 0), 2.0);
    d_off_norm = min( max( SHAPE_normalize(area_peak{j}(:, (a-1) * 2 + 1)), 0), 2.0 );

    % One more normalization to try to bring the data close to each other.
    d_on_norm = d_on_norm * mean( d_off_norm( goodbins ) ) / mean( d_on_norm( goodbins ) );

    threshold_is_a_change = 0.2;
    threshold_not_a_change = 0.5;
    min_th = 0;

    s_tot = 0; n_tot = 0;	

    % Keep track of which data points really switch
    % In the end, I just focused on the places that *should* switch. -- Rhiju
    switch_bin = 0;    
    for k = goodbins
      switch_bin(k) = 1;     
      if(str_on(k) == str_off(k)) % case 1,4
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
end

[name path] = uiputfile('Output.rdat', 'Save to RDAT file');

eterna_create_rdat_files_GUI( strcat(path,name), target_names{1}, structure, sequence, ...
            area_bsub_norm, ids, target_names, subrounds, sequence, ...
            design_names , seqpos, goodbins, data_types, [], [], ...
            min_SHAPE, max_SHAPE, threshold_SHAPE, ETERNA_score, ...
            switch_score,[], darea_bsub_norm );