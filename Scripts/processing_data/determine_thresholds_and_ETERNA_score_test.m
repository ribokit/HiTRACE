function [ min_SHAPE, max_SHAPE, threshold_SHAPE, ETERNA_score ] = determine_thresholds_and_ETERNA_score_test( data, pred, NEW_SETTINGS );
%  [ min_SHAPE, max_SHAPE, threshold_SHAPE, ETERNA_score ] = determine_thresholds_and_ETERNA_score_test( data, pred, NEW_SETTINGS );

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'NEW_SETTINGS' ); NEW_SETTINGS = 0; end;

figure(3)
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define a linear program
% to 'solve' for optimal scaling
% and baselines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length( data );

% penalties for:
% scale-factor  baseline    dev>0       dev<0      a little unpenalized slop    
f  = [   0,           0,    ones(1,n),  ones(1,n), zeros( 1, n )           ];

% Force all coefficients to be positive.
% Force coefficient of data to be at least 0.5...
LB(1) = 0.1 / mean( max(data,0.0) );
%LB(1) = 0.60 / mean( max(data,0.0) );
%LB(1) = 0.70 / mean( max(data,0.0) );

data_sort=sort(data);
data_range = abs( data_sort( floor(n*0.1)+1) - data_sort( floor(n*0.9)+1 ) );
%LB(1) = 0.3 / data_range
%UB(1) = 1.3 / data_range;
%LB(1) = 0.0;
UB(1) = inf;

% baseline
if NEW_SETTINGS
  LB(2) = -0.3;
  UB(2) =  0.3;
else
  LB(2) = -0.1;
  UB(2) =  0.1;
end

% 'slack' variable  -- positive deviations from data.
LB( 2 + [1:n] ) = zeros(1,n);
UB( 2 + [1:n] ) = inf * ones(1,n);

% 'slack' variable  -- negative deviations from data.
LB( 2 + n + [1:n] ) = zeros(1,n);
UB( 2 + n + [1:n] ) = inf * ones(1,n);

% 'slop' variables -- some range is allowed.
% for unpaired region, allow deviations up to 2-fold above 'max', and down to halfway point.
pred_high = find(  pred );

%LB( 2 + 2*n + pred_high) = -2.0;
%UB( 2 + 2*n + pred_high) =  0.5;

LB( 2 + 2*n + pred_high) =  0.0;
UB( 2 + 2*n + pred_high) =  0.0;

% for protected region, ask that we stay close to zero.
pred_low  = find( ~pred);
LB( 2 + 2*n + pred_low) = -0.0;
UB( 2 + 2*n + pred_low) = +0.0;

% how to transform from variables to prediction
Aeq = [ data', ones(n,1), eye(n,n), -eye(n,n), eye(n,n) ]; 

% what we're trying to match -- note that we're looking for equality
beq = pred;

options = optimset('Display','off');
% linear program asking for equality. (The [], [] would define greater-than, less-than constraints).
params = linprog( f, [], [], Aeq, beq, LB, UB,[],options);

%plot( params(1)*s - params(2)*b  + params( 2+ [1:n] ) - params( 2+n+ [1:n]) + params( 2+2*n+[1:n]), 'kx-' ); 
%plot( params(1)*s - params(2)*b  +  params( 2+2*n+[1:n]), 'ko-','markerfacecolor','k' ); 

subplot(2,1,1);
plot( pred, 'k','linewidth',2 );
hold on;
plot( params(1)*data + params(2), 'kx-' ); 
plot( params(1)*data + params(2) + params( 2+2*n+[1:n])', 'rx-' ); 
hold off;


subplot(2,1,2);
%params(1:2)
scale_factor =  params(1);
baseline     =  params(2);
max_SHAPE  = (1.0 - baseline)/scale_factor;
min_SHAPE =  (0.0 - baseline)/scale_factor;
threshold_SHAPE = ( 0.5 - baseline)/scale_factor;

fprintf( 'Auto thresh...\n')
fprintf( 'min      : %8.1f\n',min_SHAPE );
fprintf( 'max      : %8.1f\n',max_SHAPE );
fprintf( 'threshold: %8.1f\n', threshold_SHAPE );
bar( data,'facecolor',[0.7 0.7 0.7],'edgecolor',[0 0 0],'barwidth',0.5 );
hold on
plot( [1 length(data)], [max_SHAPE max_SHAPE], 'k' );
plot( [1 length(data)], [min_SHAPE min_SHAPE], 'k' );
plot( [1 length(data)], [threshold_SHAPE threshold_SHAPE], 'k' );
plot( (pred - baseline)/scale_factor,'color',[0.5 0 0] );


params(1:2)

ETERNA_score = 0.0;
badpt = [];
for k = 1:length( data )
  if ( pred(k) ) % should be unpaired
    if ( data(k) > (0.25*threshold_SHAPE + 0.75*min_SHAPE ) );    
      ETERNA_score  = ETERNA_score + 1;
    else
      badpt = [badpt, k ];
    %  [k,pred(k)]
    end
  else
    if ( data(k) < threshold_SHAPE )
      ETERNA_score  = ETERNA_score + 1;
    else
      badpt = [badpt, k ];
    %  [k,pred(k)]
    end
  end
end
plot( badpt, (pred(badpt)-baseline)/scale_factor, 'x','color',[0. 0.5 0],'linewidth',2 );
ylim([ (-0.2-baseline)/scale_factor (1.5-baseline)/scale_factor ]);
hold off

ETERNA_score = 100 * ETERNA_score / length( data );
fprintf( 'ETERNA score: %5.1f\n', ETERNA_score );



