function bpp_mat = bpp_pos( length, helix1, helix2 )

% Quick function to calculate matrices for plotting base pairs on bpp data
% 
% Input:
%     length = total length of RNA
%     helix1 = range of 5' sides of helices with each range on separate row
%     helix2 = range of 3' sides of helices with each range on separate row
% Output:
%     bpp_mat = matrix with 1's at helix positions and 0's elsewhere
%
% Example: To get a matrix for an RNA with a helix between 1-10 and 45-54
%          and a helix between 15-18 and 23-26 in a 100-nt RNA, use:
%          bpp_mat = bpp_pos( 100, [1 10; 15 18], [45 54; 23 26] );
%          
%
% Clarence Cheng, 3/2014

bpp_mat = zeros(length, length);

for j = 1:size(helix1,1)
    max = helix2(j,2) + helix1(j,1);
    for i = helix1(j,1):helix1(j,2)
        bpp_mat(i, max-i) = 1;
        bpp_mat(max-i, i) = 1;
    end
end