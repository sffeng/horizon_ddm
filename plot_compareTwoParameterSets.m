function [l, r, p] = plot_compareTwoParameterSets(ax, fit1, fit2, names, f1name, f2name)

global AZred AZblue

clear X1 X2
for sn = 1:length(fit1)
    for i = 1:2 % for each horizon condition
        X1{i}(:,sn) = fit1(sn).XXfit(:,i);
        X2{i}(:,sn) = fit2(sn).XXfit(:,i);
    end
end

for i = 1:length(X1)
    [r{i},p{i}] = nancorr(X1{i}', X2{i}', 'spearman');
end


clear l
for i = 1:length(X1)
    for j = 1:size(X1{i},1)
        axes(ax(j)); hold on;
        l(i,j) = plot(X1{i}(j,:), X2{i}(j,:),'.');
        
    end
end

for i = 1:size(X1{1},1)
    axes(ax(i));
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--', 'linewidth', 1)
    title(names{i}, 'fontweight', 'normal', 'interpreter', 'tex');
    xlabel(f1name)
    ylabel(f2name)
end

for i = size(X1{1},1)+1:length(ax)
    set(ax(i), 'visible', 'off')
end
% set(ax, 'xtick', [], 'ytick', [])
set(l(1,:), 'color', AZblue)
set(l(2,:), 'color', AZred)
set(l, 'linewidth', 1, 'marker', 'o', 'markersize', 5)

