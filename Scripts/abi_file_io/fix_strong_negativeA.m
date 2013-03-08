function [range] = fix_strong_negativeA(data)

% The 350-ROX ladder from GeneScan contains a strong band that migrates at 
% a smaller molecular weight than the 35-nucleotide marker. This band is
% almost always saturating and cannot be corrected by the leakage
% correction in quick_look, leading to problems with the following smooth
% baseline subtraction. Thus, use interpolation to flatten the strong
% negative peak before the smooth baseline subtraction step in quick_look.
%
% Clarence Cheng, 2013


    [local_peaks(:,1),local_peaks(:,2)] = findpeaks(data(:,2)*-1);     %find peaks and peak timepoints of negative intensity, storing in local_peaks columns 1 and 2; use lane 2 in data because that is the HEX channel
    %lp = local_peaks;

    
    stddev = std(local_peaks(:,1));
    peakmean = mean(local_peaks(:,1));       %find mean of negative peaks
    [peakmax peakindex] = max(local_peaks(:,1));     %find value and index of maximum negative peak
    local_peaks2 = local_peaks;
    local_peaks2(peakindex,1) = 0;
    %lp2 = local_peaks;
    [peakmax2 peakindex2] = max(local_peaks2(:,1));     %find value and index of second-maximum negative peak


    if peakmax - peakmean > 2*stddev
        %bias framing of interpolation range based on distance between max and second max peak to account for double-peak from saturation
        if abs(local_peaks(peakindex2,2) - local_peaks(peakindex,2)) > 30       %if maximum and second-maximum peaks are more than 30 sec apart, frame maximum peak symmetrically
            frame_pos = [local_peaks(peakindex,2)-18; local_peaks(peakindex,2)+18];
        elseif peakindex < peakindex2     %if maximum and second-maximum peaks are less than 30 sec apart and max is to the left of second max, bias asymmetric framing to the right
            frame_pos = [local_peaks(peakindex,2)-11; local_peaks(peakindex,2)+25];
        elseif peakindex > peakindex2     %if maximum and second-maximum peaks are less than 30 sec apart and max is to the right of second max, bias asymmetric framing to the left
            frame_pos = [local_peaks(peakindex,2)-25; local_peaks(peakindex,2)+11];
        end
        xi = frame_pos(1):1:frame_pos(2);    %define range to flatten by expanding around time point of negative peak by 15 ms on either side
        range = xi;
    else
        range = [0 0];      %this will prevent fix_strong_negativeB from flattening any data after leakage correction if there are no strong negative peaks in the data
    end
