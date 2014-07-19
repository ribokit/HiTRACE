function [ final_rdat,flags ] = rdat_combine( rdat_files, outfilename)
% RDAT_COMBINE - 
% [ final_rdat,flags ] = rdat_combine( rdat_files, outfilename)
%
% Combines multiple rdat files into one, whose reactivity
%   values come from an average of all profiles, weighted by uncertainty.
% Final error values come from standard deviation across traces [divided by sqrt(N)],
%   with outliers ignored
% Flags positions or entire traces that agree poorly between replicates.  Assumes that
%  all the data is the same between input rdats other than the filename,
%  reactivity, and reactivity_error.
%
%
% Inputs:
%   rdat_files  = cell of strings containing rdat filenames to combine
%   outfilename = filename for final RDAT, formatted for RMDB
%
% Outputs:
%   final_rdat  = unified .rdat file with weighted average reactivities
%   flags       = residue positions that agree poorly with others
%
% (C) T. Mann, R. Das, Stanford University, 2013.
%

if ( nargin < 1 | length( rdat_files ) < 1 ); help( mfilename); return; end;
rdats = {};
reactivities = {}; errors = {}; final_reactivity = {}; prop_err = {};

%reads in all reactivities and errors together
N_rdat = length( rdat_files );
MIN_REL_ERROR = 0.1;
N = 0;
for i = 1:N_rdat
  if ischar( rdat_files{i} )
    rdats{i} = read_rdat_file(rdat_files{i});
    rdat_tag = basename( rdat_files{i} );
  else
    rdats{i} = rdat_files{i};
    rdat_tag = rdats{i}.name;
  end
  
  reactivity = rdats{i}.reactivity;
  for j = 1:size( reactivity,2 )
    if all( reactivity(:,j) == 0 ); continue; end;
    N = N + 1;
    reactivities{N} = reactivity(:,j);
    min_error = mean( max(reactivity(:,j),0) ) * MIN_REL_ERROR;
    errors{N} = max( rdats{i}.reactivity_error(:,j), min_error );    
    rdat_file_legends{N} = rdat_tag;
  end
  
end;

for i = 1:N
  for j = 1:size( reactivities{i} )
    if ( reactivities{i}(j) < 0 ) errors{i}(j) = max( errors{i}(j), abs( reactivities{i}(j) ) ); end;
  end
  L(i) = length( reactivities{i} );
end
minL = min( L );
for i = 1:N
  reactivities{i} = reactivities{i}(1:minL);
  errors{i} = errors{i}(1:minL);
end

reactivities = cell2mat(reactivities);
errors = cell2mat(errors); %reformatting to use std function
L = size( reactivities, 1 );
sequence_full = rdats{1}.sequence;
sequence = rdats{1}.sequence(1:minL);
offset = rdats{1}.offset;
seqpos = rdats{1}.seqpos(1:minL);

 %averages data, weighted by uncertainty
[final_reactivity, final_error, flags ] = average_data_filter_outliers( reactivities, errors, seqpos, sequence, offset, trace_legends );
subplot(2,1,1);
h=title( outfilename );set(h,'interp','none','fontweight','bold');

if ~exist( 'Figures', 'dir' ) mkdir( 'Figures/' ); end;
export_fig(  ['Figures/',basename(outfilename), '.pdf'] );
save( ['Figures/',basename(outfilename), '.fig'] );


%RDAT PREPARATION using new filename, reactivity, and reactivity_error
name = rdats{1}.name;
if isempty(rdats{1}.structure); structure = []; else structure = rdats{1}.structure( 1:length( sequence_full) ); end;
if isempty(rdats{1}.annotations); annotations = {}; else annotations = rdats{1}.annotations; end;
if isempty(rdats{1}.data_annotations); data_annotations = {}; else data_annotations = rdats{1}.data_annotations; end;
if isempty(rdats{1}.trace); trace_in = []; else trace_in = rdats{1}.trace_in; end;
if isempty(rdats{1}.xsel); xsel = []; else xsel = rdats{1}.xsel; end;
if isempty(rdats{1}.xsel_refine); xsel_refine = []; else xsel_refine = rdats{1}.xsel_refine; end;
if isempty(rdats{1}.comments); comments = []; else comments = rdats{1}.comments; end;

final_rdat = output_workspace_to_rdat_file( outfilename, name, sequence_full, offset, ...
			       seqpos, final_reactivity , ...
			       structure, ...
			       annotations, data_annotations, final_error,...
			       trace_in,xsel,xsel_refine,comments );

%%%%%%%%%%%%%%%
function b = basename( tag );

remain = tag;

while ~isempty( remain )
  [token, remain ] =strtok( remain, '/' );
end

b = token;