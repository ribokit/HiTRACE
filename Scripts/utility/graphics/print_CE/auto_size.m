function [auto_H, auto_W] = auto_size(user_H, seqpos, user_W, mutpos, flg)

%
% [auto_H, auto_W] = AUTO_SIZE(user_H, seqpos, user_W, mutpos, flag);
%
% Returns page dimensions based on seqpos and mutpos sizes. User specified
%  page dimensions will be taken if non-zero and auto-size enabled.
%
%
% Input
% =====
%   user_H          Required        Provides the user specified height
%                                    (vertical page numbers). Feed in 0 to
%                                    allow auto-size output.
%   seqpos          Required        Provides the vertical data size as a
%                                    determinant of height. Usually use
%                                    seqpos, other arrays OK.
%   user_W          Required        Provides the user specified width
%                                    (horizontal page numbers). Feed in 0 to
%                                    allow auto-size output.
%   mutpos          Required        Provides the horizontal data size as a
%                                    determinant of width. Usually use
%                                    mutpos or size(d_align, 2), other 
%                                    arrays OK.
%   flag            Optional        Provides whether to perform auto_size. 
%                                    Default is 1 (YES). 0 to disable 
%                                    auto_size and take user specified 
%                                    dimensions.
%
% Output
% ======
%   auto_H                          Gives the height (vertical page
%                                    numbers).
%   auto_W                          Gives the width (horizontal page
%                                    numbers).
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('flg','var') || isempty(flg); flg = 1; end;

auto_H = user_H; auto_W = user_W;

is_auto_size_H = flg || (user_H == 0);
is_auto_size_W = flg || (user_W == 0);

if is_auto_size_H;
    H_temp = length(seqpos); 
    H_step = [75 175 250 300 350 400];
    auto_H = length(find(H_temp > H_step)) + 1;
end;
if is_auto_size_W;
    W_temp = length(mutpos);
    W_step = [55 130 185 225 260 300];
    auto_W = length(find(W_temp > W_step)) + 1;
end;

