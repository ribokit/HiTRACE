function dir_cell_out = fill_direction_cell(dir_cell_in, seqpos, which2fill)

%
% dir_cell_out = FILL_DIRECTION_CELL(dir_cell_in, seqpos, which2fill);
%
% Fills the empty element of direction_cell.
%
% Input
% =====
%   dir_cell_in         Required        Provides the direction_cell with
%                                        empty element to fill.
%   seqpos              Required        Provides all seqpos positions for a
%                                        complete direction_cell.
%   which2fill          Optional        Provides which element to fill.
%                                        Default is the first empty element 
%                                        in the cell.
%
% Output
% ======
%   dir_cell_out                        Gives the filled direction_cell.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

% default which2fill: find the 1st empty cell
if ~exist('which2fill','var') || isempty(which2fill);
    which2fill = find(cellfun('isempty', dir_cell_in));
    if length(which2fill) > 1; which2fill = which2fill(1); end;
    if isempty(which2fill); return; end;
end;

% fill the empty cell
dir_cell_out = dir_cell_in;
all_pos = cell2mat(dir_cell_in);
all_pos(dir_cell_in{which2fill}) = [];
seqpos(all_pos) = [];
dir_cell_out{which2fill} = seqpos; 



