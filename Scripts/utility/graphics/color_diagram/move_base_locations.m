function base_locations = move_base_locations(residue_locations, direction_cell, offset, square_width, square_width_correct)

% base_locations = MOVE_BASE_LOCATIONS(residue_locations, direction_cell, offset, ...
%                                      square_width, square_width_correct);
%
% Creates base_locations (off-sequence box coordinates) from residue_locations 
%   (on-sequence box coordinates) by moving each position relatively to one
%   of the four directions.
%
% Input
% =====
%   residue_locations      Required         Specifies residue_locations input. Format
%                                            in 2 x N double array.
%   direction_cell         Required         Specifies moving direction for each 
%                                            position. Format in 1 x 4 cell, with double 
%                                            array at each cell. First positon is an 
%                                            array of positions moving upward, then
%                                            downward, leftward, and rightward. 
%                                            Numbering is congruent with graph, including 
%                                            offset. Default is {[], [], [], []}, no 
%                                            moving applied.
%   offset                 Required         Specifies the numbering offset.
%   square_width           Required         Specifies square_width of each box. Format 
%                                            in double.
%   square_width_correct   Optional         Specifies correction value of square_width.
%                                            Format in double, default is 0. This offers 
%                                            chance to move the boxes in a value that is 
%                                            different from square_width if neccessary. 
%                                            Final moving distance will be (square_width
%                                            + square_width_correct).        
%
% Output
% ======
%   base_locations                          Specifies base_locations output. Format in 
%                                            2 x N double array.
%   
%
% by T47, Apr 2013.
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('direction_cell','var'); 
    base_locations = residue_locations;
    return; 
end;

if ~exist('offset','var') || isempty(offset); offset = 0; end;
if ~exist('square_width_correct','var'); square_width_correct = 0; end;


sq_width = square_width + square_width_correct;

move_up = direction_cell{1} - offset;
move_down = direction_cell{2} - offset;
move_left = direction_cell{3} - offset;
move_right = direction_cell{4} - offset;

% check if match
move_all = [move_up move_down move_left move_right];
move_max = max(move_all); move_min = min(move_all);
if length(move_all) ~= (move_max - move_min + 1) || length(move_all) ~= size(residue_locations, 2);
    fprintf('WARNING: direction_cell dimension mismatch. Please double check for missing residues.\n');
end;
for i = 1:(length(move_all) - 1)
    for j = (i + 1):length(move_all)
        if move_all(i) == move_all(j);
            fprintf('WARNING: direction_cell value invalid. Please double check for repeated residues.\n');
        end;
    end;
end;

% print summary
fprintf(['Move up: ', num2str(move_up), '.\n']);
fprintf(['Move down: ', num2str(move_down), '.\n']);
fprintf(['Move left: ', num2str(move_left), '.\n']);
fprintf(['Move right: ', num2str(move_right), '.\n']);

% move locations
base_locations(1, [move_up move_down]) = residue_locations(1, [move_up move_down]);
base_locations(2, move_up) = residue_locations(2, move_up) - sq_width;
base_locations(2, move_down) = residue_locations(2, move_down) + sq_width;

base_locations(2, [move_left move_right]) = residue_locations(2, [move_left move_right]);
base_locations(1, move_left) = residue_locations(1, move_left) - sq_width;
base_locations(1, move_right) = residue_locations(1, move_right) + sq_width;
