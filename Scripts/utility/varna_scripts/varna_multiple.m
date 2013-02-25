function []=varna_multiple(filename,sequence,DBN,DATA)
% varna_multiple(filename,sequence,DBN,DATA)

%% DATA Prep
Z=hSHAPE(DATA);
%Z=DATA;

Z = max( Z, 0 );
Z = min( Z, 1.0 );

Z=100*Z;

%%
%filename='VARNA/test.html';
fid=fopen(filename,'w');

vers=2;


fprintf(fid,'%s\n','<HTML><HEAD><META http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></HEAD><BODY>');
fprintf(fid,'%s\n','<APPLET code="VARNA.class" codebase="http://varna.lri.fr/bin" archive="VARNA.jar" width="900" height="1200">');
fprintf(fid,'%s\n','<PARAM name="algorithm" value="radiate" />');
fprintf(fid,'%s\n','<PARAM name="flat" value="true" />');
fprintf(fid,'%s%s%s\n','<PARAM name="sequenceDBN1"  value="',sequence,'"/>');
fprintf(fid,'%s%s%s\n','<PARAM name="structureDBN1" value="',DBN,'"/>');
fprintf(fid,'%s\n','<PARAM name="rows" value="6" />');
fprintf(fid,'%s\n','<PARAM name="columns" value="3" />');
for j=1:18;
fprintf(fid,'%s%i%s','<param name="colorMap',j,'" value="');

for i=1:length(Z)
    fprintf(fid,'% 6.3f%s',Z(i,j),',');
end

fprintf(fid,'%s\n','"/>');
end

fprintf(fid,'%s\n','<param name="colorMapStyle1" value="0:#0000FF;10:#0000FF;50:#FFFFFF;90:#FF0000;100:#FF0000" />'); 
fprintf(fid,'%s\n','</applet></BODY></HTML>');
