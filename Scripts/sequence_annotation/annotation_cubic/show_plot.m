function [  ] = show_plot( data, xsel, lane, prediction )
%SHOW_PLOT Summary of this function goes here
%   Detailed explanation goes here
    data_lane = data(:,lane);
    pred_lane = prediction(:,lane);
    
    plot(1:size(data,1), data_lane, 'k-');
    hold;
    
    xsel_o = xsel(pred_lane == 1);
    plot(xsel_o, data_lane(round(xsel_o)), 'ro');


end

