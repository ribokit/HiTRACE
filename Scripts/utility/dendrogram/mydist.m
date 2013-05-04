function d2 = mydist( r1, r2 )
%  d2 = mydist( r1, r2 )
%
% Testing an alternative distance metric for dendrograms
% of deep chemical profiles.
%

r1err = 0.2 + 0.2 * r1;
r2err = 0.2 + 0.2 * r2;

for j = 1:size( r2, 1 )
  r_err = sqrt( r1err.^2 + r2err(j,:).^2);
  d2(j) =  mean( (( r1 - r2(j,:) ) ./ r_err).^2 ) / 10.0;
end