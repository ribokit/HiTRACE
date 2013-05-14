function annotation_subset = annotation_type_finder(annotations, type_string)

% annotation_subset = ANNOTATION_TYPE_FINDER(annotations, type_string)
%
% Find and join all annotations of a certain type to a string.
%
%
% Input
% =====
%   annotations     Required    Provides the annotations string cell. Feed
%                                   in d_rdat.annotation or 
%                                   d_rdat.data_annotations here.
%   type_string     Required    Provides the token string. No colon needed.
%
% Output
% ======
%   annotation_subset           Provides a combined string of all entries
%                                   of specified type. Double spaces will
%                                   be added between entry strings.
%
% Example
% =======
%   annotations = {'experimentType:MutateAndMap',...
%                  'chemical:Na-HEPES:50mM(pH8.0)','chemical:MgCl2:10mM',...
%                  'chemical:S-adenosylmethionine:8.8mM', ...
%                  'temperature:24C','modifier:SHAPE'};
%   annotation_type_finder(annotations, 'chemical')
% gives answer:
%   'Na-HEPES:50mM(pH8.0)  MgCl2:10mM  S-adenosylmethionine:8.8mM'
%
%
% by T47, May 2013
%

if nargin == 0; help( mfilename ); return; end;

if ~exist('type_string','var') || isempty(type_string); return; end;

annotation_subset = '';
for i = 1:length(annotations)
    tok = strfind(annotations{i}, type_string);
    if ~isempty(tok);
        annotation_subset = [annotation_subset, '  ', ...
            annotations{i}((tok + length(type_string) + 1):end)];
    end;
end;

if ~isempty(annotation_subset); 
    annotation_subset = annotation_subset(3:end);
end;