function [f, p1_name, p2_name ] = double_exp( times, tau1, tau2 );

p1_name = 'tau1';
p2_name = 'tau2';

f = [ ones(1,length(times)); exp( -times./tau1); exp( -times./tau2 ) ];
