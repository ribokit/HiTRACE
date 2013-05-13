function [d_exp, h_length, w_length] = d_expand_divisible(d_org, H, W)


h_length = floor(size(d_org, 1)/H);
if h_length * H ~= size(d_org, 1); h_length = h_length + 1; end;
d_exp((size(d_org, 1) + 1):(h_length * H), :) = 0;

w_length = floor(size(d_org, 2)/W);
if w_length * W ~= size(d_org, 2); w_length = w_length + 1; end;
d_exp(:, (size(d_org, 2) + 1):(w_length * W)) = 0;