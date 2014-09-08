function [sens, ppv] = calc_sens_ppv_bp(test_str, native_str)

if ~exist('native_str','var') || isempty(native_str); return; end;
if length(test_str) ~= length(native_str);
    fprintf('WARNING: secondary structure length mismatch.\n');
    return;
end;

native_bps = convert_bps_to_str(convert_structure_to_bps(native_str));
str_bps = convert_bps_to_str(convert_structure_to_bps(test_str));
FP_count = length(setdiff(str_bps, native_bps));
FDR = FP_count/length(str_bps);

% str_bps = reshape(convert_structure_to_bps(test_str),length(convert_structure_to_bps(test_str))*2,1);
% native_bps = convert_structure_to_bps(native_str);
% count = [];
% for i = 1:size(native_bps,1);
%     temp = setdiff(native_bps(i,:),str_bps);
%     if isempty(setdiff(native_bps(i,:), temp));
%         count = [count; native_bps(i,:)];
%     end;
% end;
% FNR = size(count,1)/size(native_bps,1);
FNR = 1- (length(str_bps)-FP_count)/length(native_bps);

sens = 1-FNR;
ppv = 1-FDR;


