function d2 = distcorr_wrapper( r1, r2 )
%  d2 = mydist( r1, r2 )
%
% Testing an alternative distance metric for dendrograms
% of deep chemical profiles.
%

for j = 1:size( r2, 1 )
  d2(j) =  1 - distcorr( r1', r2(j,:)' );
end