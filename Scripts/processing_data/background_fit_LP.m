function [params,sum_abs_deviation] = background_fit_LP( s,y,b, penalize_negative_weight)
%  [params,sum_abs_deviation] = background_fit_LP( s,y,b, penalize_negative_weight)
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'penalize_negative_weight' ); penalize_negative_weight = 1.0; end;

b = b/mean(b);
s = s/mean(s);

% coefficients in linear programming problem. Last 2*n variables are 'slack'.
n = length( s );
f  = [   0,      0,    ones(1,n),  ones(1,n) ];

%new -- penalize 'negative' deviations more than 'positive' ones.
f  = [   0,      0,    ones(1,n),  penalize_negative_weight * ones(1,n) ];

% for 'exposed' residues, allow there to be a bigger range of values... weaker penalties.
pred_high            = find( y );
%f( 1+1+pred_high )   = 0.5;
%f( 1+1+n+pred_high ) = 0.5;

% Force all coefficients to be positive.
LB = [   0,      0,    zeros(1,n),     zeros(1,n)];
UB = [  inf,   inf, inf*ones(1,n),  inf*ones(1,n)];

%LB(1) = 0.1;
%UB(1) = 0.3;
%LB(2) = 0.10;
%UB(2) = 1.0;

Aeq = [ b, y, eye(n,n), -eye(n,n) ]; 
beq = s;

options = optimset('Display','off');
params = linprog( f, [], [], Aeq, beq, LB, UB,[],options);


plot( s, 'k','linewidth',2 );
hold on;
%plot( params(1)*s - params(2)*b  + params( 2+ [1:n] ) - params( 2+n+ [1:n]) + params( 2+2*n+[1:n]), 'kx-' ); 
%plot( params(1)*s - params(2)*b  +  params( 2+2*n+[1:n]), 'ko-','markerfacecolor','k' ); 
plot( params(1)*b + params(2)*y, 'kx-' ); 
%plot( params(1)*s,'r');
plot( params(1)*b,'b');
hold off;


sum_abs_deviation = sum(abs(params(1)*b + params(2)*y - s));

