function [auto_H, auto_W] = auto_size(user_H, seqpos, user_W, mutpos, flg)

auto_H = user_H; auto_W = user_W;

is_auto_size_H = flg || (user_H == 0);
is_auto_size_W = flg || (user_W == 0);

if is_auto_size_H;
    H_temp = length(seqpos); 
    H_step = [75 175 250 300 350 400];
    auto_H = length(find(H_temp > H_step)) + 1;
end;
if is_auto_size_W;
    W_temp = length(mutpos);
    W_step = [55 130 185 225 260 300];
    auto_W = length(find(W_temp > W_step)) + 1;
end;

