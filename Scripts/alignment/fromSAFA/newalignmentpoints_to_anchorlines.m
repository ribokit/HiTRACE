function [anchorlines, alignmentpoints_backconvert] = newalignmentpoints_to_anchorlines(alignmentpoints, anchorlines, t1lane,Profile_Size);
%
% We need to keep track of alignment lines "in the original frame".
%
% R. Das and S. Pearlman
%

numfinelanes   = size(anchorlines,2);
numanchorlines = size(anchorlines,1);

numalignmentpoints = size(alignmentpoints,1);
alignmentpoints_backconvert = 0*alignmentpoints;

if (numanchorlines>0)
    
    for j=1:numalignmentpoints
        %Find for each of the alignmentpoints, what pre-existing anchorlines
        %are close to that one.
        loweranchornum = min(find(alignmentpoints(j,t1lane) < anchorlines(:,t1lane)));     
        currentframe_low  = anchorlines(loweranchornum,t1lane);
        originalframe_low = anchorlines(loweranchornum,:);
        
        if (loweranchornum>1) 
            upperanchornum = loweranchornum - 1;
            currentframe_up  = anchorlines(upperanchornum,t1lane);
            originalframe_up = anchorlines(upperanchornum,:);
        else
            currentframe_up   = 0;
            originalframe_up = originalframe_low - currentframe_low;
        end
        
        alignmentpoints_backconvert(j,:) = alignmentpoints(j,:) .* (originalframe_up - originalframe_low) / (currentframe_up - currentframe_low) + ...
            (currentframe_up*originalframe_low - currentframe_low*originalframe_up)/ ...    
            (currentframe_up - currentframe_low);
        
    end    
    
    anchorlines( numanchorlines: numanchorlines + numalignmentpoints - 1,:)  = alignmentpoints_backconvert;
    
    lastline = numanchorlines+numalignmentpoints-1;
    
else
    anchorlines = alignmentpoints;
    lastline = numalignmentpoints;
end

anchorlines = sortrows(anchorlines,1);
anchorlines(lastline+1,:) = anchorlines(lastline,:) + Profile_Size-1 - anchorlines(lastline,t1lane);


