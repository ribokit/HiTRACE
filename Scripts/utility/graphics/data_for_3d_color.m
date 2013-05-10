function data_for_3d = data_for_3d_color( data, seqrange )
% 
% data_for_3d_color( data, seqrange )
% 
% This script is used to generate a text file or files containing
%  reactivity data to be plotted using the color_by_data command in
%  pymol_rhiju.py.
% 
% The data are formatted in the text file so that the first column is the
%  residue position and the second column is the corresponding reactivity.
%  Each column in the data array corresponds to one text file.
% 
% When run, the script will prompt a user entry for each column in the data
%  array, allowing the user to set the name of the text file that is output
%  for each column of data. For example, if the data array contains lanes
%  for DMS and SHAPE, one may enter "DMS" for the first prompt and "SHAPE"
%  for the second prompt, and the text files "DMS.txt" and "SHAPE.txt" will
%  be collected in the folder "Data_for_3D_color".
%
% 
% Inputs:
% -------
%     data          = MatLab array containing reactivity data in columns
%     seqrange      = Range covering reactivities in the data array that
%                     correspond to the sequence in the .pdb structure
%                     (e.g. [24:96])
% 
% 
% (c) Clarence Cheng, May 2013

if ~exist( 'data','var' ); fprintf('Must enter a dataset'); return; end;
if ~exist( 'seqrange','var' ); seqrange = length(data); end;


mkdir('Data for 3d color');
cd('Data for 3d color');


[x,y] = size(data);
for i = 1:y
    data_for3d = [];
	data_for3d(:,2) = data(seqrange,i);
    data_for3d(:,1) = 1:1:length(seqrange);
	data_for3dT = transpose(data_for3d);        % transpose the matrix 
    fid = -1;
    msg = '';  
    while fid < 0 
        disp(msg);
        filename = input('Condition: ', 's');   % enter desired filename (e.g. SHAPE)
        filename = strcat(filename,'.txt');
        [fid,msg] = fopen(filename,'w' );
    end
    fid = fopen(filename,'w');                  % open a text file for reading in the reactivities
    fprintf(fid,'%1.0f  %5.4f\n',data_for3dT);
    fclose(fid);
end




