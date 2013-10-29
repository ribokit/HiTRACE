function  [data_average, data_average_err ] = average_data( data, data_err );
% AVERAGE_DATA: Averages multiple profiles/reactivities, taking into account their estimated errors.
%
% [data_average, data_average_err ] = average_data( data, data_err );
%
% Inputs:
%  data     = Array with multiple profiles
%  data_err = [optional] Error on the profiles
%
% (C) R. Das, 2013.
%

if ( nargin < 1 ); help( mfilename); return; end;

if ~exist( 'data_err', 'var' ); 
  data_average = mean( data, 2 );
else
  data_average = sum( data./data_err.^2, 2) ;
  sum_weights  = sum( 1./data_err.^2, 2);

  data_average = data_average ./ sum_weights;    
  data_average_err = sqrt( 1 ./sum_weights );
end
