function [d_correct,modification_fractions] = attenuation_corrector_v2( d_in )
% [d_correct,modification_fractions] = attenuation_corrector_v2( d_in )
%

if nargin == 0;  help( mfilename ); return; end;

for k = 1:size(d_in,2);
  
  d = d_in(:,k);

  % On calculating modification fractions, 
  % skip first value, as it might be unextended primer
  % also ignore data right at end [leave out last 3 nts]
  modification_fractions(k) = sum(d(2 : end-5)) / sum( d(2:end) );
  fprintf( '%2d modification fraction: %8.4f\n', k, modification_fractions(k));
  
  d_ratios(:,k) = 0 * d;
  for j = 1:length( d );
    d_ratios(j,k) =  d(j)/sum( d(1:end) );
  end;

  d_correct(:,k) = 0 * d;
  d_correct(end,k) = 1;
  for j = 1:length(d_correct) - 1;
    d_correct(j,k) = d(j) / sum(d(j:end));
  end;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some visual feedback
clf
subplot(1,2,1);
image( d_ratios*1000 );
subplot(1,2,2);
image( d_correct*1000 );
colormap( 1 - gray(100) );

return;
