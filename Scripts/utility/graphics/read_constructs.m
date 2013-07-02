function construct_names = read_constructs( filename );

if nargin == 0;  help( mfilename ); return; end;

fid = fopen( filename );
k=0;
while ~feof(fid)
  k=k+1;
  construct_name = fgetl(fid);
  construct_name = strrep( construct_name, '"', '' );
  %fprintf(1,'%s\n',construct_name );
  construct_names{k} = construct_name;
end

fclose(fid);
