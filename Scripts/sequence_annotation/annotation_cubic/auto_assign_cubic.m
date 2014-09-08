function [ xsel, escore, peak, score, bonus, peak_bonus_idx, gamma_] = auto_assign_cubic( data, prediction, sequence, param )
% AUTO_ASSIGN_CUBIC - Tool for rapid manual assignment of bands in electropherograms through cubic dynamic programming.
%
% [ xsel, peak, score, bonus, peak_bonus_idx, gamma_] = auto_assign_cubic(data, prediction, sequence, param );
%
%
% Input:
%  data           = matrix of aligned electrophoretic traces.
%  prediction     = prediction matrix.
%  sequence       = subsequence of RNA to be used.
%  param     WINDOW_SIZE = the size of window that shifts with sequence along the time series;
%            window_jump = 
%            begin = the starting position of window. 'auto' can be chosen.
%            end = the ending position of window. 'auto' can be chosen.
%            peak_bonus_weight
%            peak_range = maximum range in which a peak is allowed to be matched;
%            min_gamma
%            gamma_add
%            primary_weight
%
% Output:
% xsel      = positions of bands across all lanes.
% escore    = E-score
% peak      = information on peak matching
% etc.
%
% (C) Seungmyung Lee, 2014
%

%% automated primary_weight setting
if (param.primary_weight == 'auto')
    num_band_location = sum(max(data')');
    num_primary_band = sum(data(:,end));
    
    primary_weight = (num_primary_band / num_band_location)^0.5;
else
    primary_weight = param.primary_weight;
end

%% automated begin and end setting (only last lane)
    q23_max = 0;
    q1_max = 0;
    q4_max = 0;
    startline = 1;
    endline = round(size(data,1) * 0.9);
    
    for i = 1:round(size(data,1)/4)
        if (q1_max < data(i,end))
            q1_max = data(i,end);
            startline = i;
        end
    end
    
    for i = round(3*size(data,1)/4):size(data,1)
        if (q4_max < data(i,end))
            q4_max = data(i,end);
            endline = i;
        end
    end
    data(endline+1:end,end) = data(endline,end);
    
        
    for i = round(3*size(data,1)/8):round(3*size(data,1)/4)
        q23_max = max(q23_max, data(i,end));
    end

    if (param.begin == 'auto')
        param.begin = startline - round(param.WINDOW_SIZE/2);
    end
    param.end = endline;
    if (param.window_jump == 'auto')
        for j = 0:size(prediction,1)-1
            if (prediction(end-j,end) == 1)
                tailing_zeros = j;
                break;
            end
        end
        param.window_jump = (param.end-param.begin-param.WINDOW_SIZE)/(size(prediction,1) - j/2);
    end
        
%% modify prediction matrix
    % the first row of prediction matrix always contains only 1
    prediction (1,:) = 1;

%% build peak bonus matrix - (peak_bonus)
    % build gamma matrix
    delta1_1 = (data(3:size(data,1)-2,:) - data(1:size(data,1)-4,:))/2;
    delta1_2 = data(3:size(data,1)-2,:) - data(2:size(data,1)-3,:);
    delta2_1 = data(4:size(data,1)-1,:) - data(3:size(data,1)-2,:);
    delta2_2 = (data(5:size(data,1),:) - data(3:size(data,1)-2,:))/2;
    delta1 = max(delta1_1,delta1_2);
    delta2 = min(delta2_1,delta2_2);

    gamma = (delta1 - delta2).*(delta1_2 > 0).*(delta2_1 < 0);  % -gamma
    
    gamma = [zeros(2,size(data,2));gamma;zeros(2,size(data,2))];
        
    peak_bonus = zeros(size(gamma));
    for i = 1:size(data,2)
        optimal_peak_num=round(sum(prediction(:,i))*2);
        if (param.min_gamma == 'auto')
            sorted_gamma = sortrows(gamma, -i);
            min_gamma = sorted_gamma(optimal_peak_num, i);
        else
            min_gamma = param.min_gamma;
        end
        gamma_i = gamma(:,i);
        
        gamma_add = param.gamma_add * mean(gamma_i(gamma_i > min_gamma));
        peak_bonus(:,i) = (gamma_i + gamma_add) .* (gamma_i > min_gamma);
    end
    peak_bonus (:,end) = peak_bonus (:,end) * primary_weight * (size(peak_bonus,2)-1);
    peak_bonus (:,1:end-1) = peak_bonus (:,1:end-1) * (1-primary_weight);
    
    peak_bonus = peak_bonus * param.peak_bonus_weight;
    
    gamma_ = peak_bonus;

%% generate distance bonus - (dist_bonus)
    % based on the density function of Gaussian with mean 12 and stdv 4
    x = [1:60];
    dist_bonus = normpdf(x, 1, param.window_jump/3) / normpdf(1, 1, param.window_jump/3);
    
%% proximity component of matching peak bonus - (prox_peak)
    % based on the density function of Gaussian with mean 0 and stdv 2
    prox_peak = normpdf(1:param.peak_range+1,1.5, param.window_jump/6) / normpdf(1, 1, param.window_jump/6);
   
    
%% Preprocesesing for peak_bonus window for the last lane
    lanes = size(prediction,2);
    peak_bonus_idx = zeros(sum(peak_bonus(:,lanes)>0),2);
    idx = 1;
    for i = 1:size(peak_bonus, 1)
        if (peak_bonus(i,lanes) > 0)
            peak_bonus_idx(idx,1) = i;
            peak_bonus_idx(idx,2) = peak_bonus(i,lanes);
            idx = idx + 1;
        end
    end
    assert (idx == size(peak_bonus_idx,1) + 1, 'peak_bonus_idx dimension error');
    
    peak_bonus_window = zeros(size(data,1),2);
    for i = 1:size(peak_bonus_window,1)
        peak_bonus_window(i,1) = sum(peak_bonus_idx(:,1) < i-5) + 1;
        peak_bonus_window(i,2) = sum(peak_bonus_idx(:,1) <= i+5);
    end
    

%% DYNAMIC PROGRAMMING
    bonus = (-10000)*ones(size(prediction,1), param.WINDOW_SIZE, size(peak_bonus_idx,1)+1);
    trace1 = zeros(size(prediction,1), param.WINDOW_SIZE, size(peak_bonus_idx,1)+1);
    trace2 = zeros(size(prediction,1), param.WINDOW_SIZE, size(peak_bonus_idx,1)+1);
    
    % seq = 1
    window_begin{1} = param.begin + 1;
    for j1 = param.begin + round(param.WINDOW_SIZE/2)
        for j2 = peak_bonus_window(j1,1):peak_bonus_window(j1,2)
            if (prediction(1,lanes) == 1)
                dist = abs(peak_bonus_idx(j2,1) - j1);
                bonus(1,j1 - window_begin{1} + 1,j2) = peak_bonus_idx(j2,2) * prox_peak(dist + 1);
            end
        end
        bonus(1,j1 - window_begin{1} + 1,:) = bonus(1,j1 - window_begin{1} + 1,:) + calc_match_bonus (peak_bonus, j1, prediction(1,:), prox_peak);
    end

    % seq > 1
    for seq = 2:size(prediction,1)
        ratio = sum(prediction(seq,:))/size(prediction,2);
        if (ratio > 0)
            bonus_boost = (1/ratio)^0.5;
        else
            bonus_boost = 1;
        end
        
        ideal_dist = round(param.window_jump) + 1;%13
        dist_weight = 1;
        if (sequence(seq-1) == 'G')
            ideal_dist = round(param.window_jump*3/4) + 1;%10;
            dist_weight = 1/1;
        end
        window_begin{seq} = param.begin + round((seq-1)*param.window_jump) + 1;
        for j1 = window_begin{seq}:window_begin{seq}+param.WINDOW_SIZE-1
            if (j1 >= param.end && seq <= size(prediction,1))
                continue;
            end
            for i1 = max(j1-30,max(window_begin{seq-1},j1-param.WINDOW_SIZE)):min(window_begin{seq-1}+param.WINDOW_SIZE-1,j1-1)
                % fill in bonus(:,:,end).
                for i2 = [peak_bonus_window(i1,1):peak_bonus_window(i1,2) size(peak_bonus_idx,1)+1]
                    % 'for loop' is not necessary in fact. to be revised later...
                    discr_dist = abs(j1 - i1 - ideal_dist);
                    new_candidate = bonus(seq-1, i1 - window_begin{seq-1} + 1, i2) + dist_weight*dist_bonus(discr_dist + 1);
                    if (new_candidate > bonus(seq, j1 - window_begin{seq} + 1, end))
                        % When i2 is in the range of j1, (seq,j1,i2) must
                        % be filled instead of (seq,j1,end)
                        if (peak_bonus_window(j1,1) <= i2)
                            bonus(seq, j1 - window_begin{seq} + 1, i2) = new_candidate;
                            trace1(seq, j1 - window_begin{seq} + 1, i2) = i1;
                            trace2(seq, j1 - window_begin{seq} + 1, i2) = i2;
                        else
                            bonus(seq, j1 - window_begin{seq} + 1, end) = new_candidate;
                            trace1(seq, j1 - window_begin{seq} + 1, end) = i1;
                            trace2(seq, j1 - window_begin{seq} + 1, end) = i2;
                        end
                    end
                end
                % fill in bonus(:,:,1:end-1). no need to fill if prediction is zero
                if (prediction(seq,lanes) == 1)
                    for j2 = peak_bonus_window(j1,1):peak_bonus_window(j1,2)
                        for i2 = [peak_bonus_window(i1,1):min(peak_bonus_window(i1,2),j2-1) size(peak_bonus_idx,1)+1]
                            dist = abs(peak_bonus_idx(j2,1) - j1);
                            discr_dist = abs(j1 - i1 - ideal_dist);
                            new_candidate = bonus(seq-1, i1 - window_begin{seq-1} + 1, i2) + dist_weight*dist_bonus(discr_dist + 1) + peak_bonus_idx(j2,2) * prox_peak(dist + 1);
                            if (new_candidate > bonus(seq, j1 - window_begin{seq} + 1, j2))
                                bonus(seq, j1 - window_begin{seq} + 1, j2) = new_candidate;
                                trace1(seq, j1 - window_begin{seq} + 1, j2) = i1;
                                trace2(seq, j1 - window_begin{seq} + 1, j2) = i2;
                            end
                        end
                    end
                end
           
            end
            bonus(seq, j1 - window_begin{seq} + 1, :) = bonus(seq, j1 - window_begin{seq} + 1, :) + calc_match_bonus (peak_bonus, j1, prediction(seq,:), prox_peak);
        end
    end

    
%% BACKTRACE
    xsel = zeros(size(prediction,1), 1);
    peak = zeros(size(prediction,1), 1);
    % find max for seq = end
    score = 0;
    
    for i = 1:size(bonus,2)
        for j = 1:size(bonus, 3)
            if (bonus(end,i,j) > score)
                score = bonus(end,i,j);
                xsel(end) = i + window_begin{size(bonus,1)} - 1;
                peak(end) = j;
            end
        end
    end
    
    for seq = size(prediction,1)-1:-1:1
        xsel(seq) = trace1(seq+1, xsel(seq+1) - window_begin{seq+1} + 1, peak(seq+1));
        peak(seq) = trace2(seq+1, xsel(seq+1) - window_begin{seq+1} + 1, peak(seq+1));
    end
    

%% E-SCORE
    n1 = sum(peak == max(peak) & prediction(:,end)==1);
    n2 = sum(diff(xsel) <= param.window_jump/4); 
    n3 = sum(diff(xsel) > param.window_jump * 2);
    
    escore = 1 - max(n1/sum(prediction(:,end)), (n2 + n3)/(size(prediction,1)-1));
    
end