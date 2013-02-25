function []=makepdf(xsel,d,sequence_names,FigTITLE,offset,spacing,bufferspace,xticknames,xticking)
% Function to make a file that will be ready to be saved in a pdf format.
% Takes your aligned data and presents it in a useable format with sequence
% assignments and numbering.


if ~exist( 'FigTITLE','var' )
  FigTITLE='Aligned data';
end

if ~exist( 'spacing','var' )
    spacing=5;
end

if ~exist( 'offset','var' )
    offset=0;
end

if ~exist( 'bufferspace','var' )
    bufferspace=100;
end

if ~exist( 'xticknames','var' )
    xticknames=1:size(d,2);
end

if ~exist( 'xticking','var' )
    xticking=1:size(d,2);
end

if ~exist('JUST_PLOT_SEQUENCE')
  JUST_PLOT_SEQUENCE = 0;
end

mod_shift=spacing-mod(offset,spacing);
%Creates ticks on y-axis and x-axis
sequence_ticks=(0:spacing:length(xsel)-2*spacing);
sequence_ticks=sequence_ticks(end:-1:1);
yticking=xsel(end:-1:1);
yticking=yticking(sequence_ticks+1+mod_shift);

name_form='%03.0f';
prec=ceil(log10(length(sequence_names)));

name_form(3)=num2str(prec);


for i=1:length(sequence_ticks)
    seqnum=sequence_ticks(i)+offset+mod_shift;
    sequence_names_ticks(i,:)=char([ sequence_names(sequence_ticks(i)+1+mod_shift) num2str(seqnum,name_form)]);
end

%Replot realigned image
ymin = round(min(min( xsel))) - 150;
ymax = round(max(max( xsel))) + bufferspace;

%d=5*d_align;
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/3 -100 scrsz(3)/2 scrsz(4)])
image( .5*d);
sides=max(xticking)*.01;
axis([0.5-sides max(xticking)+0.5+sides ymin ymax]);
set(gca,'xtick',xticking,'xticklabel',xticknames);
set(gca,'ytick',yticking,'yticklabel',sequence_names_ticks,'fontsize',8);
xticklabel_rotate;
title( FigTITLE,'fontsize',10,'fontweight','bold' );

colormap(  1 -gray(100) );