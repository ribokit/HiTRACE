function overlap_Z_to_Rx (reactivity, Z, is_invert)
% OVERLAP_Z_TO_RX (reactivity , Z, [is_invert]);
%
% Overlaps Z score matrix to reactivity matrix in image.
%
% =Input=
%   reactivity      Pre-saved image file of reactivity, e.g. '1.jpg'.
%                       Saving could be done by imsave. Format in string.
%   Z               Z score matrix, e.g. flipud(Z*-2). Format in double.
%   [is_invert]     Optional flag for whether to invert the value of Z 
%                       matrix. Format in double, 0 for no inverting, 1
%                       for yes, default 0.
%
% by T47, Mar 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('is_invert','var'); is_invert=0; end;
if is_invert; Z = Z*-1; end;
    
reference = imread(reactivity);
figure, imshow(reference);
hold on;

data = Z;
h = imshow(data,[]);
hold off;

colormap jet;
alphamap = zeros(size(reference,1),size(reference,2));

for i = 0:size(data,1)-1
    for j = 0:size(data,2)-1
        if(~(data(i+1,j+1) == 0))
            alphamap(i+1,j+1) = 0.75;
        end;
    end;
end;
set(h, 'AlphaData', alphamap);