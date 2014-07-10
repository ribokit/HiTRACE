function plot_area_pred (area_pred, clr, mk_sz, ln_wd, mk_sb)

if ~exist('mk_sz','var') || isempty(mk_sz); mk_sz = 3; end;
if ~exist('ln_wd','var') || isempty(ln_wd); ln_wd = 0.5; end;
if ~exist('mk_sb','var') || isempty(mk_sb); mk_sb = 'o'; end;

for i = 1:size(area_pred,2)
    for j = 1:size(area_pred,1)
        if area_pred(j,i) == 1 && i ~= j;
            plot(j,i,mk_sb,'color',clr,'markersize',mk_sz,'linewidth',ln_wd); hold on;
        end;
    end;
end;
