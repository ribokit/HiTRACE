function [dsub, alpha] = background_subtracter( s, b, norm_index, M, P );

% background column may be specific, rather than background profile
if ( size(b,1) == 1 )
  b = s(:,b);
end
   
if ~exist( 'norm_index' ) | length(norm_index)==0
  nres = size( s, 1 );
  eat_in = min( (nres - 20)/2, 10 );
  norm_index = [ eat_in: nres - eat_in];
end

if ~exist( 'M' )
  M = 100;
end

if ~exist( 'P' )
  P = 1000;
end

% loop through.
for i = 1:size( s,2)
  [dsub(:,i), alpha(i)] = ...
      background_sub_individual( s(:,i), b, norm_index, M , P );
end
fprintf( '\n' )

make_image( s, dsub );

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dsub, alpha ] = background_sub_individual( s, b, norm_index, M, P );

iterations = 20;

s_orig = s;
b_orig = b;
s = s( norm_index );
b = b( norm_index );

% One term keeps background-subtracted data uniform (its
% coefficient is 1.0 ).

% Second term is a positivity penalty -- keep scaled background below signal!
P = 1000.0;

% Third term tries to maximize background level. Units on term
%  should be signal-units^2.
M_scale = M*sum( b.^2 );

% initial background scaling
%alpha = 1.0;
alpha=0.1;

for i = 1:iterations

  g = ( ( s - alpha* b ) < 0 ); %positivity penalty
  alpha = ( sum( s.*b.*( 1 + P*g ) ) - (sum( b ) * sum( s )) + M_scale ) / ...
	  ( sum( b.*b.*( 1 + P*g ) ) - sum( b )^2 );

  alpha = abs( alpha ); % stay positive!
  
end

%bar([ s, alpha*b, s - alpha*b] );
%plot( 1:length(s), s,'b' ); hold on
%plot( 1:length(s), alpha*b,'r' ); hold off

dsub = s_orig - alpha *b_orig;
fprintf( 1, 'Background scaling: %6.2f\n',alpha );

