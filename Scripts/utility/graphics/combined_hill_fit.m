function [deviation,pred_fit] = combined_hill_fit( params, conc, data )
%
% [deviation,pred_fit] = combined_hill_fit( params, conc, data )
%

if nargin == 0;  help( mfilename ); return; end;

n_hill = abs(params( 1 ));
midpt  = 10.^params( 2 );

%numdata = ( length( params ) - 2 ) /2;
numdata = size( data, 1 );
start_val = params(2+[1:numdata] );
end_val = params(2+numdata+[1:numdata]);

pred =  1 - ( 1 ./ ( 1 + (conc./midpt).^n_hill ) );
for k = 1:numdata
  pred_fit(k,:) = pred * (end_val(k) - start_val(k) ) + start_val(k);
end

deviation = sum( sum( ( pred_fit - data ).^2 ) );

% Fudge factor to induce hill coeffiecient of 2
% Removed by CVL
%+ ...    100*( n_hill - 2).^2 ;
