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
function [peak_out] = FindPeakNew_Ver3(Dyin, yin, Th_intensity, Th_width, Th_adjacent, debug_on)
                                                                
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
    
%     [pos_length, neg_length] = compute_length(sgn, ind(u)); %jk
    pos_length = 1;neg_length = 1;
    peak(ind(u)) =  (pos_length + neg_length);
            
end


peak_out = zeros(2, length(yin));
peak_out(1,:) = peak;
peak_out(2,ind) = yin(ind);

peak_out = flipud(peak_out); % 1st-row is intensity, 2nd-row is width

if debug_on
 ind
end



% ----------------------
% 4) If threshold is set
%   then consider it
% ----------------------    
if Th_width(1) ~= 0
    
    for u = 1:length(ind)
        if peak_out(2,ind(u)) < Th_width(2)
            peak_out(1,ind(u)) = 0;   
            peak_out(2,ind(u)) = 0;              
        end
    end
    
end

if Th_intensity(1) ~= 0
    
    for u = 1:length(ind)
        if peak_out(1,ind(u)) < Th_intensity(2)
            peak_out(1,ind(u)) = 0;   
            peak_out(2,ind(u)) = 0;  
        end
    end
    
end


if Th_adjacent(1) ~= 0
        
    while(1) % clean-up operation should be repeated
        x = 1:length(peak_out(1,:));
        ind_nonzero = x(peak_out(1,:) ~= 0);       
        
        diff_ind = diff(ind_nonzero);
        if (min(diff_ind) >= Th_adjacent(2))
            break;
        else
            [ind_nonzero, peak_out] = remove_adjacent_peak(ind_nonzero, diff_ind, peak_out, Th_adjacent(2));                    
        end
    end
    
end

% if Th_adjacent(1) ~= 0
%     x = 1:length(peak_out(1,:));
%     ind_nonzero = x(peak_out(1,:) ~= 0);
%     %-------------------------------
%     while(1)
%         diff_ind = diff(ind_nonzero);
%     
%         if (min(diff_ind) > Th_adjacent(2)) 
%             break;
% %         elseif length(diff_ind) == 0
% %             break;
%         else
%             debug_on = 1;
%             [ind_nonzero, peak_out] = remove_adjacent_peak(ind_nonzero, diff_ind, peak_out, Th_adjacent(2), debug_on);
%         end
%         
%         if debug_on
%             ind_nonzero
%         end
%     end      
%     %-------------------------------
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   remove_adjacent_peak
%   
%
%   2008 bshim
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ind, peak_out] = remove_adjacent_peak(ind, diff_ind, peak_out, th)

x = 1:length(diff_ind);
y = x(diff_ind < th);
for v = 1:length(y)
    
    k2 = ind(y(v));
    k1 = ind(y(v)+1);
    
    if peak_out(1,k1) >= peak_out(1,k2)
        peak_out(1,k2) = 0;
        peak_out(2,k2) = 0;
    else
        peak_out(1,k1) = 0;
        peak_out(2,k1) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   remove_adjacent_peak
%   
%
%   2008 bshim
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ind, peak_out] = remove_adjacent_peak(ind, diff_ind, peak_out, th, debug_on)
%          
% 
%     for u = 1:length(diff_ind)
%        
%         if diff_ind(u) <= th            
%             
%             if peak_out(1,ind(u)) <= peak_out(1,ind(u+1))       
%                 ind = [ind(1:u-1) ind(u+1:end)];
%                 peak_out(1, ind(u)) = 0;
%                 peak_out(2, ind(u)) = 0;            
%                 
%                 if debug_on
%                     case1 = 1
%                 fprintf('1) %2.1f(@%d) > %2.1f(@%d)\n',peak_out(1,ind(u+1)),ind(u+1),peak_out(1,ind(u)),ind(u)); 
%                 end
%             else    
%                 if debug_on
%                     case2 = 2
%                 end
%                 if (u+2) >= length(diff_ind)    
%                     ind = ind(1:u);
%                     if (u+1) <= length(diff_ind)
%                         peak_out(1, ind(u+1)) = 0;
%                         peak_out(2, ind(u+1)) = 0;  
%                     end
%                 else
%                     ind = [ind(1:u) ind(u+2:end)];     
%                     peak_out(1, ind(u+1)) = 0;
%                     peak_out(2, ind(u+1)) = 0;
%                     
%                     if debug_on
%                     fprintf('2) %2.1f(@%d) > %2.1f(@%d)\n',peak_out(1,ind(u)),ind(u),peak_out(1,ind(u+1)),ind(u+1)); 
%                     end
%                 end
%             end
%             break;
%         end
%     end
%    
% 
%     for u = 1:length(diff_ind)
%        
%         if diff_ind(u) <= th
%             
%             if peak_out(1,ind(u)) <= peak_out(1,ind(u+1))       
%                 ind = [ind(1:u-1) ind(u+1:end)];
%                 peak_out(1,ind(u)) = 0;
%                 peak_out(2,ind(u)) = 0;             
%             else    
%                 if (u+2) >= length(diff_ind)      
%                     peak_out(1, ind(u+1)) = 0;
%                     peak_out(2, ind(u+1)) = 0;
%                     ind = ind(1:u);
%                 else                           
%                     peak_out(1, ind(u+1)) = 0;
%                     peak_out(2, ind(u+1)) = 0;
%                     ind = [ind(1:u) ind(u+2:end)];                
%                 end
%             end
%             break;
%         end
%     end


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
    elseif neg_indx > 14950
        break;
    end                       %%%  
    
end




