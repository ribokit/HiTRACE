function [ids, target_names, subrounds, sequences, design_names ] = parse_EteRNA_sequence_file( eterna_sequence_file );

ids = [];
target_names = {};
subrounds = [];
sequences = {};
design_names = {};

fid = fopen( eterna_sequence_file, 'r' );

count = 0;
while ~feof( fid )
  line = fgetl( fid );
  if length(line) > 6 

    [t,r] = strtok( line, char(9) );
    id = str2num( t );
    if length(t) > 4 & ~isempty( id )
      [t,r] = strtok( r, char(9) );
      target_name = t;

      [t,r] = strtok( r, char(9) );
      subround = str2num(t);
      
      [t,r] = strtok( r, char(9) );  
      if ( t(1) == ' ' ); t = t(2:end); end;   % this is weird, but jee started adding a space to the beginning of the sequence.
      sequence = t;
      
      design_name = r(2:end); % the first character is tab because of the way strtok works

      if target_name(1) == ' '; target_name = target_name(2:end); end
      if design_name(1) == ' '; design_name = design_name(2:end); end
      
      count = count + 1;
      ids(count) = id;
      target_names{count} = target_name;
      subrounds(count) = subround;
      sequences{count} = sequence;
      design_names{count} = design_name;
      
    end    
  end
end