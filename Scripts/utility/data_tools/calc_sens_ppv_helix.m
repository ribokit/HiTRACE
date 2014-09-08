function [sens, ppv] = calc_sens_ppv_helix(test_str, native_str)

if ~exist('native_str','var') || isempty(native_str); return; end;
if length(test_str) ~= length(native_str);
    fprintf('WARNING: secondary structure length mismatch.\n');
    return;
end;

native_stems = parse_stems(native_str);
% for i = 1:length(native_stems)
%     if size(native_stems{i},1) <3;
%         native_stems{i} =[];
%     end;
% end;
test_stems = parse_stems(test_str);
% fprintf(['t',test_str,'\n']);
% fprintf(['n',native_str,'\n']);

TP_count = 0;
for i = 1:length(native_stems);
    if isempty(native_stems{i}); continue; end;
    native_pairs = convert_bps_to_str(native_stems{i});
    for j = 1:length(test_stems);
        test_pairs = convert_bps_to_str(test_stems{j});
        overlap = setdiff(native_pairs, test_pairs);
        if length(overlap) < length(native_pairs);
            %fprintf([num2str(length(overlap)),'\t',num2str(length(native_pairs)),'\n']);
            if length(overlap) <= length(native_pairs)/2;
                TP_count = TP_count + 1;
            end;
            break;
        end;
    end;
end;
FDR = 1 - TP_count /length(test_stems);

% FNR_count = 0; FNR_flag = 0;
% for i = 1:length(native_stems);
%     if isempty(native_stems{i}); continue; end;
%     native_pairs = native_stems{i};
%     for j = 1:length(test_stems);
%         test_pairs = reshape(test_stems{j},size(test_stems{j},1)*2,1);
%         overlap = setdiff(native_pairs, test_pairs);
%         if ~isempty(setdiff(native_pairs, overlap));
%             FNR_flag = 1; break;
%         end;
%     end;
%     if ~FNR_flag; FNR_count = FNR_count + 1; end;
%     FNR_flag = 0;
% end;
% FNR = FNR_count/length(native_stems);
FNR = 1- TP_count/ length(native_stems);

sens = 1-FNR;
ppv = 1-FDR;
