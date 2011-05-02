%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FindPeak.m
%
%   2008 bshim
%   
%   input vector: 
%                 Dyin         : Differential data
%                 yin          : Original data
%                 Th_width     : Th_width(1) <- if use turn-ON(1) else turn-OFF(0)
%                                Th_width(2) <- threhold value
%                 Th_intensity : Th_intensity(1) <- if use turn-ON(1) else turn-OFF(0)
%                                Th_intensity(2) <- threhold value
%                 Th_adjacent  : if peak is found in the adjacent region,
%                                choose only one
%                 debug_on     : debugging plots are turned-on
%
%   output vector : 
%                 peak_out : peak information (location/intensity/width)
%         (ex) peak = [0 0 0 72 0 0 7 0 0 0 14 0 0 0 5 ; % intensity info.  
%                      0 0 0 29 0 0 4 0 0 0 18 0 0 0 7]; % width info.                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peak_out] = FindPeak_Ver3(Dyin, yin, debug_on)
                                                                
% ----------------------
% 1) Assign +/- sign
% ----------------------
sgn = zeros(size(Dyin));
sgn(find(Dyin > 0)) = 1;
sgn(find(Dyin < 0)) = -1;


% ----------------------
% 2) Find zero-crossing
% ----------------------
peak = zeros(size(Dyin));

for v = 3:length(Dyin)-1

    %%% 
    %%% peak condition (integer zero crossing)
    %%%
    if ((sgn(v-1) > 0) & (sgn(v+1) < 0)) 
        [maxval maxind] = max([yin(v-2) yin(v-1) yin(v) yin(v+1)]); % <--- adjust the mistake of zero-crossing 
        peak((v-2) + maxind - 1) = 1;    
        %pkval1 = (v-2) + maxind - 1
    end                                  
    
end

for v = 2:length(Dyin)-1

    %%% 
    %%% peak condition (fractional zero crossing)
    %%%
    if ((sgn(v) > 0) & (sgn(v+1) < 0)) 
        
if 1 % ##########################
    [maxval maxind] = max([yin(v-1) yin(v) yin(v+1)]); 
    peak((v-1) + maxind - 1) = 1; 
    %pkval2 = (v-1) + maxind - 1
end            
% if 0 % ##########################     
%         if (yin(v) >= yin(v+1)) %
%             peak(v) = 1;        % <--- adjust the mistake of zero-crossing
%         else                    % 
%             peak(v+1) = 1;      %
%         end
% end
    
    end                               
    
end
peak(1) = 0; % false peak


if debug_on
    x = 1:length(sgn);
    figure;plot(x,peak,'-');
    grid;legend('peak');title('debugging');
    axis([0 1000 -2 2]);

    figure;plot(x,sgn,'-');
    grid;legend('sign');title('debugging');
    axis([0 1000 -2 2]);
end

% ----------------------
% 3) Compute peak-width
% ----------------------
ind = find(peak == 1);
for u = 1:length(ind)    
    
    [pos_length, neg_length] = compute_length(sgn, ind(u));
    peak(ind(u)) =  (pos_length + neg_length);
            
end


peak_out = zeros(2, length(yin));
peak_out(1,:) = peak;
peak_out(2,ind) = yin(ind);

peak_out = flipud(peak_out); % 1st-row is intensity, 2nd-row is width





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   compute_length
%
%   2008 bshim
%   
%   input vector: sin (sign-input), indx
%   output: pos_length, neg_length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos_length, neg_length] = compute_length(sgn, indx)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos_length = 0;pos_indx = indx-1;
while sgn(pos_indx) > 0
    
    pos_length = pos_length + 1;
    pos_indx = pos_indx - 1;
    
    if pos_indx == 0 %%%
        break;       %%% in case indx touches the 0                   
    end              %%%
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neg_length = 0;neg_indx = indx+1;
while sgn(neg_indx) < 0
    
    neg_length = neg_length + 1;    
    neg_indx = neg_indx + 1;
    
    if neg_indx > length(sgn) %%%
        break;                %%% in case indx is bigger than the vector length
    end                       %%%  
    
end




