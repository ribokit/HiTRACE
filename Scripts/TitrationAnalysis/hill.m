function [f, p1_name, p2_name ] = hill( conc, K1, n );

p1_name = 'K';
p2_name = 'nHill';
pred = (conc/K1).^n ./ (1 + (conc/K1).^n );
f = [1-pred; pred ];
