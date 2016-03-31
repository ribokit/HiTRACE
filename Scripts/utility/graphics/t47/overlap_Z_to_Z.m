function overlap_im = overlap_Z_to_Z (Z1, Z2, fuse_method, im_method)
% overlap_im = OVERLAP_Z_TO_Z (Z1, Z2, [fuse_method], [im_method]);
%
% Overlaps two Z score matrices using imfuse function. Colored green
%   and purple.
%
% =Input=
%   Z1, Z2              Two Z score matrices [double] to compare.
%   [fuse_method]       Optional argument for imfuse blending method.
%                           Format in string, default 'falsecolor'.
%   [im_method]         Optional argument for image choice. Format in
%                           double, 0 for image, 1 for imagesc. 
%                           Default 0.
%
% =Output= 
%   overlap_im          [RGB image] matrix of Z score overlay.
%
% by T47, Mar 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('fuse_method','var'); fuse_method = 'falsecolor'; end;
if ~exist('im_method','var');  im_method = 0; end;

c = imfuse(Z1, Z2, fuse_method);
if im_method;
    imagesc(c);
else
    image(c);
end;

overlap_im = c;
