function [profiles_align, profiles_combine,anchorlines] = align_userinput_saveanchorlines(profiles,numfinebins,t1num,anchorlines,pointertoaxes);
% [profiles_align, profiles_combine]  = align_userinput(profiles,numfinebins,t1num);
%
% R. Das and S. Pearlman, update to make it reversible, June 2, 2004.
%

if ~exist( 'pointertoaxes' ) pointertoaxes = gca; end;

% Set the amount of change for the contrast
ContrastScale = sqrt(1.5);
global maxprof;
global grayscaleon;
global renderSqrt;
renderSqrt = 0;
grayscaleon = 1;

axes(pointertoaxes);

%The "T1 lane" is a reference lane which is kept fixed through the whole
%alignment.
if (nargin<3) t1num=1; end;
t1lane = numfinebins*(t1num - 1)+1;

[Profile_Size,numprofiles_fine] = size(profiles);
currentaxis = [1 numprofiles_fine 1 Profile_Size]; 
axis(currentaxis);zoomedin = 0;
stopalign = 0;
plottitle = [...
        'left button: start anchorline       right button (apple-click): zoom in/out \newline',...
        '1,2,C: adjust colorscale        Q: do alignment            Z:finished'];

numanchorlines = size(anchorlines,1);
numlanes = size(profiles,2);
if (numanchorlines>0)
    profiles_align = calculate_alignedprofiles(profiles,anchorlines,t1lane);
else
    profiles_align = profiles;
end

%plot the image!
hold off;
if(renderSqrt == 1)
    image(sqrt(abs(profiles_align)));
else
    image(abs(profiles_align));
end

 maxprof = squeeze(max(max(profiles_align)))/160;
 grayscaleon = 1;
 setcolormap(grayscaleon,maxprof);
set(gca,'xtick',(numfinebins+1)/2:numfinebins:numlanes,'xticklabel',1:numlanes/numfinebins,'xminortick','on');

%If there were predefined anchorlines, draw them too -- as horizontal
%lines.
for k=1: numanchorlines
    hold on ; h_old(k) = plot([1 numprofiles_fine], [anchorlines(k,t1lane) anchorlines(k,t1lane)],'r'); hold off
end

while (stopalign == 0)
    pickedpeak = 0;
    stopinput = 0; 
    count = 0;
    while (stopinput == 0)
        title(plottitle);
        [xselpick,yselpick,button] = ginput(1);
        switch(button)
            case 1
                [xsel,ysel,pickedpeak,grayscaleon,maxprof] = quickpolyxline(profiles,0,yselpick,0,grayscaleon,maxprof,ContrastScale);
                if (pickedpeak>0) 
                    count = count + 1;
                    alignmentpoints(count,:) = interp1(xsel,ysel,1:numprofiles_fine);   hold on                                     
                    laneline(count)    = plot(1:numprofiles_fine,alignmentpoints(count,:),'r');    
                    firstsymbol(count) = plot(1,alignmentpoints(count,1),'ro'); 
                    lastsymbol(count)  = plot(numprofiles_fine,alignmentpoints(count,numprofiles_fine),'ro'); hold off
                end    
                %
                % Hey Sam! Don't let the user  cross streams!
                %  when interp is called in calculate_alignedprofiles, make
                %  sure its not trying to stretch something that's zero.                
            case 3        
                if (zoomedin == 1)
                    currentaxis = [1 numprofiles_fine 1 Profile_Size]; 
                    axis(currentaxis);
                    zoomedin = 0;
                else
                    yselzoom = yselpick;
                    ymin = max(round(yselzoom - Profile_Size/10),0);
                    ymax = min(round(yselzoom + Profile_Size/10),Profile_Size);
                    currentaxis=[1 numprofiles_fine ymin ymax ];
                    axis(currentaxis);
                    zoomedin = 1;
                end
        end
        switch char(button)
            case '1'
                maxprof = maxprof/ContrastScale;
                setcolormap(grayscaleon, maxprof);    
            case '2'
                maxprof = maxprof*ContrastScale;
                setcolormap(grayscaleon, maxprof);    
            case {'c','C'}
                grayscaleon = 1 - grayscaleon
                setcolormap(grayscaleon, maxprof);    
            case {'e','E'}
                % There are two kinds of lines that may need to be erased:
                %   already straightebed anchorlines and
                %    newly drawn anchorlines ("alignmentpoints"). 
                % They need to be treated separately.
                    
		
		closeoldanchor = Profile_Size^2 + 1;
		if length( anchorlines ) > 0 
		  [closeoldanchor, anchornum] = min( (anchorlines(:,round(xselpick)) - yselpick).^2);
		end
		closenewanchor = Profile_Size^2;
                if (count>0) [closenewanchor, anchornum_new] = min( (alignmentpoints(:,round(xselpick)) - yselpick).^2); end
                if (closeoldanchor < closenewanchor) %old anchorline closer
%                        set(h_old(anchornum),'visible','off');
                        anchorlines = anchorlines([1:anchornum-1  anchornum+1:numanchorlines],:);                        
                        stopinput = 1; %Better show user what happens when you getrid of the old one.    
                    else %new anchorline closer.
                     set(laneline(anchornum_new)   ,'visible','off');
                     set(firstsymbol(anchornum_new),'visible','off');
                     set(lastsymbol(anchornum_new) ,'visible','off');
                     newnumbers = [1: anchornum_new-1   anchornum_new+1:count];    
                     alignmentpoints = alignmentpoints(newnumbers,:);
                     laneline    = laneline(newnumbers);
                     firstsymbol = firstsymbol(newnumbers);
                     lastsymbol  = lastsymbol(newnumbers);
                     count= count-1;
		end
            case {'q','Q'}
                stopinput = 1;
                if (pickedpeak == 0)
                    stopalign = 1;
                end;
            case {'z','Z'}
                stopinput = 1;
                stopalign = 1;
        end
    end
    if (count>0)
        anchorlines = newalignmentpoints_to_anchorlines(alignmentpoints, anchorlines, t1lane,Profile_Size);
    end
    alignmentpoints = [];numanchorlines = size(anchorlines,1);
    
    profiles_align = calculate_alignedprofiles(profiles,anchorlines,t1lane);        
    hold off;
    if(renderSqrt == 1)
        image(sqrt(abs(profiles_align)));
    else
        image(abs(profiles_align));
    end
    setcolormap(grayscaleon,maxprof);
    hold on
    axis(currentaxis);
    for k=1: numanchorlines
        hold on ; h_old(k) = plot([1 numprofiles_fine], [anchorlines(k,t1lane) anchorlines(k,t1lane)],'r'); hold off
    end
    
end

profiles_combine=combinelanes(profiles_align,numfinebins);
title('Done aligning gel.');

