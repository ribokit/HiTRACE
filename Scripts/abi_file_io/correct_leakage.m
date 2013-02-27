function alldata_correct = correct_leakage( alldata, leakage_matrix );
%  CORRECT_LEAKAGE: applies 4x4 inversion to separate channels for ABI sequencers run on BigDye settings.
%
%    alldata_correct = correct_leakage( alldata, leakage_matrix );
%
% (C) R. Das, 2008-2011

alldata_correct = {};
if nargin == 0;  help( mfilename ); return; end;

invmatrix = inv( leakage_matrix );
for j = 1: size( alldata,2 )
  data = alldata{j};
  alldata_correct{j} = data * invmatrix;
end