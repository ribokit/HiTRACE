function area_pred = generate_area_pred (sequence, structure, offset, data_types, numlanes)

%
% area_pred = GENERATE_AREA_PRED(sequence, structure, offset, data_types, numlanes);
%
% Generates area_pred predicted based on given structure of sequence and data_types.
%
%
% Input
% =====
%   sequence        Required        Provides the sequence.
%   structure       Required        Provides the secondary structure in dot
%                                       format.
%   offset          Required        Provides the offset of sequence numbering.
%   data_types      Required        Provides the data_types for prediction.
%   numlanes        Required        Provides the number of lanes.
%
%
% by T47, Rhiju Das, May 2013.
%

if nargin == 0; help( mfilename ); return; end;

if isempty( data_types ) % no input.
    area_pred = zeros( length(sequence), numlanes );
elseif iscell( data_types ) % input is a matrix
    for m = length( data_types )+1 : numlanes; data_types{m} = ''; end; % pad to number of lanes.
    area_pred = get_area_pred( sequence, data_types, offset, structure );
else
    if size( data_types, 1 ) == length( sequence )
        area_pred = data_types;
    else
        fprintf( 'Problem: input area_pred/data_types does not have same size [%d] as sequence [%d]\n', size(data_types,1),  length( sequence ) );
        return;
    end
    data_types = [];
end