function [f, p1_name, p2_name ] = one_two( conc, K1, K2 );

p1_name = 'K1';
p2_name = 'K2';

f = [ ones(1,length(conc));   (conc/K1);   (conc/K2).*(conc/K1) ];
[dummy, Z]  = meshgrid( ones(1,3), sum( f ) );
f = f./Z';
