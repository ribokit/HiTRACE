function eterna_create_rdat_files_GUI( rdat_filename, target_name, structure, sequence, area_bsub, ...
    ids, target_names, subrounds, sequences, design_names , ...
    seqpos, goodbins, data_type, added_salt, bad_lanes, min_SHAPE, max_SHAPE, threshold_SHAPE, EteRNA_score, EteRNA_switch,comments, darea_bsub, trace_data,xsel );
% eterna_create_rdat_files( rdat_tag, target_name, structure, sequence, area_bsub, eterna_sequence_files, seqpos, data_type, added_salt, bad_lanes );

% read eterna sequence files ... this is our reference for eterna IDs, design names, etc.


if(isempty(trace_data))
    trace = [];
else
    for i = 1:length( sequence)
        trace(:,[1:6] + (i - 1) * 6) = trace_data{i};
    end    
end

if(isempty(xsel))
    xsels = [];
else
    for i = 1:length( sequence)
        xsels(:,i) = xsel{i};
    end
    xsels = xsels(goodbins,:);
end


if(size(structure,1) > 1)
    str = structure(1,:);
    for i = 2:size(structure,1)
        str = [str, ' ', structure(i,:)];
    end
    structure = str;
end

% find sequence matches
match = zeros( 1, length(sequence) );
for i = 1:length( sequence )
  found_match = 0;
  for j = 1:length( sequences )
    if strcmp( sequence{i}, sequences{j} )
      match(i) = j;
      break;
    end
  end
  if match(i) == 0
    fprintf( 'Error: could not find a match of sequence %d to any of the input eterna_sequence_files\n', i );
    return;
  end
end

% get ready for output.
area_out = [];
darea_out = [];

if ~exist( 'darea_bsub')
    fprintf( 'WARNING!WARNING! DID NOT SPECIFY DAREA_bsub!\n');
    darea_bsub = area_bsub * 0.0;
end

sequence_out = '';
for m = 1:length(sequence{1}); sequence_out = [sequence_out,'X']; end;

if(isempty(trace_data))
    annotations = {'experimentType:StandardState','chemical:Na-HEPES:50mM(pH8.0)','temperature:24C','processing:overmodificationCorrection','processing:backgroundSubtraction'};
else
    annotations = {'experimentType:StandardState','chemical:Na-HEPES:50mM(pH8.0)','temperature:24C'};
end

count = 0;
offset = 0; % for EteRNA generally, there are no weird sequence offsets.
mutpos = [];

v = zeros(length(sequence),1);

v(bad_lanes) = 1;
bad_lanes = v;

for i = 1:length( sequence)
  m = match(i);

  if ( size( area_bsub{i}, 1 )  ~= length(seqpos) )
    fprintf( 'Error: Number of points in area_bsub [%d] does not match length of seqpos[%d]\n', size(area_bsub{i},1), length(seqpos) );
    return;
  end
  
  
  
  for j = 1:length( data_type )
    if (strcmp( data_type{j}, 'SHAPE' ) || strcmp( data_type{j}, 'DMS' ))
      count = count + 1;
      area_out(:,count) = area_bsub{i}(:,j);
      darea_out(:,count) = darea_bsub{i}(:,j);

      
     
      data_annotations{count} = { ...
		    ['EteRNA:design_name:', strrep( strrep( design_names{m},' ','_' ), '%', '%%' ) ],...		    
            ['EteRNA:target:',strrep(target_names{m},' ','_') ],...
	  ['EteRNA:ID:',num2str(ids(m)) ],...
	  ['EteRNA:subround:',num2str(subrounds(m))], ...
	  ['EteRNA:score:EteRNA_score:',num2str(EteRNA_score{i,j},'%6.1f')]};
      if min_SHAPE{i,j} ~= 0
	data_annotations{count} = [ data_annotations{count}, ...		    
	  { ['EteRNA:score:min_', data_type{j},':',num2str(min_SHAPE{i,j},'%7.3f')],...
	    ['EteRNA:score:max_', data_type{j},':',num2str(max_SHAPE{i,j},'%7.3f')],...
	    ['EteRNA:score:threshold_', data_type{j},':',num2str(threshold_SHAPE{i,j},'%7.3f')] } ];
      end
      data_annotations{count} = [ data_annotations{count}, ...		          
		    {['EteRNA:score:switch_score:', num2str(EteRNA_switch(i),'%6.1f')],...
		     ['modifier:',data_type{j}]} ];
      if ( ~isempty(added_salt) && length(added_salt{m}{j} > 0 ) )
          data_annotations{count} = [ data_annotations{count}, added_salt{m}{j} ]; 
      end
      data_annotations{count} = [ data_annotations{count}, ['sequence:',sequence{i}] ];
      
      mutpos( count ) = NaN;
      
      if(bad_lanes(i))
        data_annotations{count} = [ data_annotations{count}, ['warning:badQuality'] ];
      end
    elseif(~isempty(trace_data))
      count = count + 1;
      area_out(:,count) = 0;
      darea_out(:,count) = 0;
      
      data_annotations{count} = { ...
		    ['EteRNA:design_name:', strrep( strrep( design_names{m},' ','_' ), '%', '%%' ) ],...		    
            ['EteRNA:target:',strrep(target_names{m},' ','_') ],...
	  ['EteRNA:ID:',num2str(ids(m)) ],...
	  ['EteRNA:subround:',num2str(subrounds(m))]};
	      data_annotations{count} = [ data_annotations{count}, ['modifier:',data_type{j}] ];
      if ( ~isempty(added_salt) && length(added_salt{m}{j} > 0 ) )
          data_annotations{count} = [ data_annotations{count}, added_salt{m}{j} ]; 
      end
      data_annotations{count} = [ data_annotations{count}, ['sequence:',sequence{i}] ];
      
      mutpos( count ) = NaN;
      
      if(bad_lanes(i))
        data_annotations{count} = [ data_annotations{count}, ['warning:badQuality'] ];
      end
    end
  end
end

filename = rdat_filename; %[rdat_tag,'_000',num2str(rdat_index),'.rdat'];



output_workspace_to_rdat_file(filename,target_name,...
			      sequence_out, offset, seqpos( goodbins ), area_out( goodbins, : ), ...
			      mutpos, structure, ...
			      annotations, data_annotations,...
			      darea_out( goodbins, : ),trace,[],xsels', comments );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function match = match_bad_lanes( i, j, bad_lanes );
match = 0;
if isempty( bad_lanes )
    return;
end
gp = find( bad_lanes(:,1) == i );
if ~isempty( gp ) & ( sum( bad_lanes(gp,2) == j ) > 0 )
  match = 1
end
