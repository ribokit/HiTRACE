function plot_area_pred (area_pred, clr, mk_sz, ln_wd, mk_sb, is_con)

if ~exist('mk_sz','var') || isempty(mk_sz); mk_sz = 3; end;
if ~exist('ln_wd','var') || isempty(ln_wd); ln_wd = 0.5; end;
if ~exist('mk_sb','var') || isempty(mk_sb); mk_sb = 's'; end;
if ~exist('is_con','var') || isempty(is_con); is_con = 0; end;

area_pred = triu(area_pred,1) + tril(area_pred,-1);
if ~is_con;
    for i = 1:size(area_pred,2)
        for j = 1:size(area_pred,1)
            if area_pred(j,i) == 1;
                plot(j,i,mk_sb,'color',clr,'markersize',mk_sz,'linewidth',ln_wd); hold on;
            end;
        end;
    end;
else
    expand = 3;
    for k = 1:expand;
        idx = find(area_pred==1);
        idxs = [idx+1,idx-1, idx+length(area_pred)*1,idx-length(area_pred)*1];
        idxs = idxs(idxs>0 & idxs <=length(area_pred)^2);
        area_pred(idxs) = 1;
    end;
    area_pred = triu(area_pred,expand) + tril(area_pred,-expand);
    contour(1:size(area_pred,1), 1:size(area_pred,2), area_pred, 1, 'color',clr,'linewidth',ln_wd); hold on;
end;