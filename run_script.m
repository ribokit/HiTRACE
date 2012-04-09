%% output
% FIXED_ETERNA - for checking correlation with traditional ETERNA score
% score2 - score based option 1
% score_edit(i, j , k) - score based edit distance and eterna score
%                  (i-th sequence, j-th weight composite, k - SHAPE/DMS)
% score_hamming(i, j , k) - score based hamming distance and eterna score
%                  (i-th sequence, j-th weight composite, k - SHAPE/DMS)

%% initialization
addpath(genpath('../'));

alt_structure = '.....((((......((((....)))).....))))....................................';
structure = '..........................(((((((............)))))))....................';

load 'data.mat';

alt_lane = [2 4];

offset = 0;
dist = 20;
which_sets = 1:12;

START = 6;
END = 6;
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
  
  others = 0;
  for k = 1:40
      if(pred(k) == 1)
          if( data_norm(k) > 0.125)
              others = others + 1;
          end
      else
          if( data_norm(k) < 0.5)
              others = others + 1;
          end
      end
  end
  
  % ETERNA score from fixed threshold.
  FIXED_SCORE(j,a) = others / (52 - START - END) * 100; 
  
  % traditional ETERNA score.
  [min_SHAPE{j,a}, max_SHAPE{j,a}, threshold_SHAPE{j,a}, ETERNA_score(j,a), d_bin{j}(:,a)] = determine_thresholds_binarization_and_ETERNA_score( data_norm, pred );
  
  
  %pause;
  end
end

%% pre-calculation for switch score from base-pair comparison -- option 1 (adaptive/fixed threshold)
for i = 1:12
    for a = 1:2
        str_on = all_area_pred{i}(:, (a-1) * 2 + 2);
        str_off = all_area_pred{i}(:, (a-1) * 2 + 1);
        
        d_on_norm = area_bsub_norm{i}(:, (a-1) * 2 + 2);
        d_off_norm = area_bsub_norm{i}(:, (a-1) * 2 + 1);
      
        % Adaptive threshold calculation
        threshold = mean([threshold_SHAPE{i,(a-1) * 2 + 1} threshold_SHAPE{i,(a-1) * 2 + 2}]);
        min_th = mean([min_SHAPE{i,(a-1) * 2 + 1} min_SHAPE{i,(a-1) * 2 + 2}]);

        % fixed threshold
%         threshold = 0.5;
%         min_th = 0;


        s = [];

        for j = 1:52
            if(str_on(j) == str_off(j)) % case 1,4
                if(abs(d_on_norm(j)-d_off_norm(j)) < threshold)
                    s(j) = 1;
                else
                    s(j) = 0;
                end
            elseif(str_on(j) > str_off(j)) % case 2
                dif = d_on_norm(j) - d_off_norm(j);
                if(dif > threshold)
                    s(j) = 1;
                elseif ( dif < min_th)
                    s(j) = 0;
                else
                    s(j) = 2* dif;
                end
            else % case 3
                dif = d_on_norm(j) - d_off_norm(j);
                if(dif < -threshold)
                    s(j) = 1;
                elseif ( dif > min_th)
                    s(j) = 0;
                else
                    s(j) = -2* dif;
                end
            end
        end

        switch_score(i,a) = sum(s((START+1):(52-END)) / length(s((START+1):(52-END)))) * 100;
    end
end

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
        
        
        % option 2
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
        
        % option 1
        score2(j,a) = geomean([f g switch_score(j,a)/100]) * 100;
    end
end

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