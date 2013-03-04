function []=varna_fig_error(filename,sequence,DBN,DATA,err,Display)
%
%  varna_fig_error(filename,sequence,DBN,DATA,err,Display)
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'err' )
 err=[];
end

if ~exist( 'Display' )
 Display=[];
end


%% DATA Prep
Z=hSHAPE(DATA);
Z=DATA;

Z = max( Z, 0 );
Z = min( Z, 1.0 );

Z=100*Z;

mark_err=find(err>=Display(1) & err<=Display(end));
res_err=err(mark_err)-Display(1);
Z(res_err)=-5;

%%

%% VARNA creation
%filename='VARNA/test.html';
fid=fopen(filename,'w');

vers=1;


fprintf(fid,'%s\n','<HTML><HEAD><META http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></HEAD><BODY>');
fprintf(fid,'%s\n','<APPLET code="VARNA.class" codebase="http://varna.lri.fr/bin" archive="VARNA.jar" width="900" height="700">');
fprintf(fid,'%s%s%s\n','<PARAM name="sequenceDBN"  value="',sequence,'"/>');
fprintf(fid,'%s%s%s\n','<PARAM name="structureDBN" value="',DBN,'"/>');
fprintf(fid,'%s\n','<PARAM name="algorithm" value="radiate" />');
fprintf(fid,'%s','<param name="colorMap" value="');

for i=1:length(Z)
    fprintf(fid,'% 6.3f%s',Z(i),',');
end

fprintf(fid,'%s\n','"/>');
fprintf(fid,'%s\n','<PARAM name="colorMapStyle" value="-5:#949494;0:#0000FF;10:#0000FF;40:#FFFFFF;60:#FFFFFF;90:#FF0000;100:#FF0000";/>'); 
fprintf(fid,'%s\n','<PARAM name="flat" value="true" />');
fprintf(fid,'%s%s%s\n','<PARAM name="title" value="',filename,'"/>');
fprintf(fid,'%s\n','</applet></BODY></HTML>');
