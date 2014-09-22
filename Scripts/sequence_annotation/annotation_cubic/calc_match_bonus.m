function [ bonus ] = calc_match_bonus( peak_bonus, pos, prediction, prox_peak, peak_range)
    bonus = 0;
    for lane = 1:size(prediction,2)-1
        if (prediction(lane) == 0)
            continue
        end
        tmp = zeros(peak_range*2+1, 1);
        for i = -peak_range:peak_range
            if (pos + i > 0)
                tmp(i + peak_range + 1) = peak_bonus(pos + i, lane) * prox_peak(abs(i)+1);
            end
        end
        bonus = bonus + max(tmp);
    end
    ratio = sum(prediction)/size(prediction,2);
    if (ratio > 0)
        bonus = bonus * (1/ratio)^0.5;
    end
end

