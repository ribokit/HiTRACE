function construct_names = read_constructs( filename );

fid = fopen( filename );
k=0;
while ~feof(fid)
  k=k+1;
  construct_name = fgetl(fid);
  %fprintf(1,'%s\n',construct_name );
  construct_names{k} = construct_name;
end

fclose(fid);
