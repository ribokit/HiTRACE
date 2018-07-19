function [ xsel, escore, peak, score, bonus, peak_bonus_idx, gamma_, processed_data] = auto_assign_cubic( data, prediction, sequence, param, data_types )
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
%            exist_tailpeak
%  data_types     = an array of data types
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
        param.begin = startline; % - round(param.WINDOW_SIZE/2);
    end
    if (param.end == 'auto')
        param.end = endline;
        end_auto = 1;
    else
        end_auto = 0;
    end

    if (param.window_jump == 'auto')
        param.WINDOW_SIZE = round(max((param.end-param.begin)/10, 2*(param.end-param.begin)/size(prediction,1)));
    end
    param.window_jump = (param.end-param.begin-param.WINDOW_SIZE/3)/(size(prediction,1) - 1);

ideal_spacing = (param.end - param.begin) / (size(prediction,1) - 1);
    
%% delete empty lanes and choose the primary lane among the lefts
ddTTP = 0;
if exist( 'data_types', 'var' ) & length( data_types ) > 0
    for i = 1:min(size(data,2), length(data_types))
        if strcmp(data_types{i},'ddTTP')
            ddTTP = i;
            break;
        end
    end
end

if ddTTP > 0
    % if ddTTP exists, use it as primary lane
    primary_lane = ddTTP;
else
    % if not, choose the lane with greatest number of bands as primary lane
    band_num = sum(prediction);
    [ll primary_lane] = max(band_num);
    
    % primary_lane = size(data,2);
end

primary_data = data(:, primary_lane);
primary_pred = prediction(:, primary_lane);
data(:, primary_lane) = data(:, end);
prediction(:, primary_lane) = prediction(:, end);
data(:, end) = primary_data;
prediction(:, end) = primary_pred;

% remove empty lanes
peak_num = sum(prediction);
data = data(:,peak_num>0);
prediction = prediction(:,peak_num>0);

[num_pixels num_lanes] = size(data);


%% normalize and convlute data
data = window_normalize( data, ideal_spacing * 2)/2;

x = (1:num_pixels)';
xmid = round(median(x));
ideal_gaussian = get_gaussian( x, xmid, ideal_spacing/6 );

n_fft_row = 2*num_pixels - 1;
gaussian_fft =  fftn( ideal_gaussian, [n_fft_row num_lanes] );
data_fft = fftn( data, [n_fft_row num_lanes]);
f = real( ifft2( gaussian_fft .* data_fft ) );
f = f( xmid + (0:num_pixels-1), : );
data = f;
    
%% modify prediction matrix
    % the first row of prediction matrix always contains only 1
    prediction (1,:) = 1;
    if (param.exist_tailpeak == 1)
        prediction (end+1,:) = 2;
    end

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
    
    % deal with clustered peaks
    for j = 1:size(data,2)
        for i = 1:size(data,1)
            if (gamma(i,j) == 0)
                continue;
            end
            
            for k = [max(1, i-round(ideal_spacing/2)):i-1 i+1:min(size(data,1), i+round(ideal_spacing/2))]
                if (gamma(k,j) > 0 && data(k,j) < data(i,j))
                    gamma(k,j) = 0;
                end
            end
        end
    end
    
    for i = 1:size(data,2)
        optimal_peak_num=min(round(sum(prediction(:,i))*2), size(prediction,1));
        %optimal_peak_num=max(optimal_peak_num, round(sum(gamma(:,i)>0)/2));
        if (optimal_peak_num == 0)
            continue;
        end
        if (param.min_gamma == 'auto')
            sorted_gamma = sortrows(gamma, -i);
            min_gamma = sorted_gamma(optimal_peak_num, i);
        else
            min_gamma = param.min_gamma;
        end
        gamma_i = gamma(:,i);
        
        gamma_add = param.gamma_add * mean(gamma_i(gamma_i > min_gamma));
        data_i = data(:,i);
        peak_bonus(:,i) = (gamma_i + gamma_add) .* (gamma_i > min_gamma) ;%.* data_i / mean(data_i(gamma_i > min_gamma));
    end
    
    %{
    peak_bonus = zeros(size(gamma));
    peak_bonus(:,1:end-1) = data(:,1:end-1);
    for j = size(data,2):size(data,2)
        peaks = getpeaks( data(:,j), 'THRESHOLD', min([ mean(data(:,j)), median(data(:,j)), max(data(:,j))*0.05 ]) );
        for i = 1:length(peaks)
            peak_bonus(peaks(i),j) = data(peaks(i),j);
        end
    end
    %}
    
    peak_bonus (:,end) = peak_bonus (:,end) * primary_weight * (size(peak_bonus,2)-1);
    peak_bonus (:,1:end-1) = peak_bonus (:,1:end-1) * (1-primary_weight);
        
    peak_bonus = peak_bonus * param.peak_bonus_weight;
    
    gamma_ = peak_bonus;

%% generate distance bonus - (dist_bonus)
    % based on the density function of Gaussian with mean 12 and stdv 4
    x = [1:round(ideal_spacing * 3)];
    %dist_bonus = normpdf(x, 1, ideal_spacing/3) / normpdf(1, 1, ideal_spacing/3);
    dist_bonus = 1 - (2*x/ideal_spacing).^2;
    
%% proximity component of matching peak bonus - (prox_peak)
    % based on the density function of Gaussian with mean 0 and stdv 2
    if (param.peak_range == 'auto')
        param.peak_range = round(ideal_spacing/2)+1;
    end
    prox_peak = normpdf(1:ideal_spacing+1,1, ideal_spacing/6) / normpdf(1, 1, ideal_spacing/6);
   
    
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
        peak_bonus_window(i,1) = sum(peak_bonus_idx(:,1) < i-round(ideal_spacing/2)) + 1;
        peak_bonus_window(i,2) = sum(peak_bonus_idx(:,1) <= i+round(ideal_spacing/2));
    end
    

%% DYNAMIC PROGRAMMING
    bonus = (-10000)*ones(size(prediction,1), param.WINDOW_SIZE, size(peak_bonus_idx,1)+1);
    trace1 = zeros(size(prediction,1), param.WINDOW_SIZE, size(peak_bonus_idx,1)+1);
    trace2 = zeros(size(prediction,1), param.WINDOW_SIZE, size(peak_bonus_idx,1)+1);
    
    % seq = 1
    window_begin{1} = param.begin - round(param.WINDOW_SIZE/3);
    for j1 = param.begin
        % in case such j2 does not exist
        bonus(1,j1-window_begin{1} + 1, :) = 0;
        %bonus(1, 1, :) = 0;
        for j2 = peak_bonus_window(j1,1):peak_bonus_window(j1,2)
            if (prediction(1,lanes) == 1)
                dist = abs(peak_bonus_idx(j2,1) - j1);
                bonus(1,j1 - window_begin{1} + 1,j2) = peak_bonus_idx(j2,2) * prox_peak(dist + 1);
            end
        end
        bonus(1,j1 - window_begin{1} + 1,:) = bonus(1,j1 - window_begin{1} + 1,:) + calc_match_bonus (peak_bonus, j1, prediction(1,:), prox_peak, param.peak_range);
    end

    % seq > 1
    for seq = 2:size(prediction,1)
        ratio = sum(prediction(seq,:))/size(prediction,2);
        if (ratio > 0)
            bonus_boost = (1/ratio)^0.5;
        else
            bonus_boost = 1;
        end
        
        ideal_dist = round(ideal_spacing*10/9) + 0;%13
        dist_weight = 1;
        if (sequence(seq-1) == 'G')
            ideal_dist = round(ideal_spacing*3/5) + 0;%10;
            dist_weight = 1/1;
        end
        window_begin{seq} = round(window_begin{1} + (seq-1)*param.window_jump) + 1;
        for j1 = max(1,window_begin{seq}) : min( window_begin{seq}+param.WINDOW_SIZE-1,size( peak_bonus_window,1))
            if (j1 >= param.end && seq <= size(prediction,1))
                continue;
            end
            for i1 = max(1, max(j1-round(param.window_jump*2),max(window_begin{seq-1},j1-param.WINDOW_SIZE))):min(window_begin{seq-1}+param.WINDOW_SIZE-1,j1-1)
                % fill in bonus(:,:,end).
                peak_bonus_window(i1,1);
                peak_bonus_window(i1,2);
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
            bonus(seq, j1 - window_begin{seq} + 1, :) = bonus(seq, j1 - window_begin{seq} + 1, :) + calc_match_bonus (peak_bonus, j1, prediction(seq,:), prox_peak, param.peak_range);
        end
    end

    
%% BACKTRACE
    xsel = zeros(size(prediction,1), 1);
    peak = zeros(size(prediction,1), 1);
    % find max for seq = end
    score = 0;
    
    if (end_auto == 1)
        for i = 1:size(bonus,2)
            for j = 1:size(bonus, 3)
                if (bonus(end,i,j) > score)
                    score = bonus(end,i,j);
                    xsel(end) = i + window_begin{size(bonus,1)} - 1;
                    peak(end) = j;
                end
            end
        end
    else
        xsel(end) = param.end-1;
    end
    
    for j = 1:size(bonus,3)
        if (bonus(end, param.end - window_begin{size(bonus,1)}, j) > score)
            score = bonus(end, param.end - window_begin{size(bonus,1)}, j);
            peak(end) = j;
        end
    end
    
    for seq = size(prediction,1)-1:-1:1
        xsel(seq) = trace1(seq+1, xsel(seq+1) - window_begin{seq+1} + 1, peak(seq+1));
        peak(seq) = trace2(seq+1, xsel(seq+1) - window_begin{seq+1} + 1, peak(seq+1));
    end
    

%% E-SCORE
    n1 = sum(peak == max(peak) & prediction(:,end)==1);
    n2 = sum(diff(xsel) <= ideal_spacing/4); 
    n3 = sum(diff(xsel) > ideal_spacing * 2);
    
    escore = 1 - max(n1/sum(prediction(:,end)), (n2 + n3)/(size(prediction,1)-1));
    
    xsel = xsel';
    
    processed_data.data = data;
    processed_data.pred = prediction;
end