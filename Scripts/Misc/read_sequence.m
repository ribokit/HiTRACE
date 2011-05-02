function sequence = read_sequence( filename );

fid = fopen( filename );
sequence = fgetl(fid);
fprintf(1,'%s\n',sequence );
fclose(fid);
