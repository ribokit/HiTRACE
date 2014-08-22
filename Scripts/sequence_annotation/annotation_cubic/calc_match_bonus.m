function [ bonus ] = calc_match_bonus( peak_bonus, pos, prediction, prox_peak )
% calculate peak matching bonus for non-primary lanes
    bonus = 0;
    for lane = 1:size(prediction,2)-1
        if (prediction(lane) == 0)
            continue
        end
        tmp = zeros(11, 1);
        for i = -5:5
            tmp(i + 6) = peak_bonus(pos + i, lane) * prox_peak(abs(i)+1);
        end
        bonus = bonus + max(tmp);
    end
    ratio = sum(prediction)/size(prediction,2);
    if (ratio > 0)
        bonus = bonus * (1/ratio)^0.5;
    end
end