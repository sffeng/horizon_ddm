function plot_RTregression_v1(ax, sub, RTmin, RTmax)

global AZred AZblue

for sn = 1:length(sub)
    % say 1 is left
    dR = sub(sn).o1(:,4) - sub(sn).o2(:,4);
    dI = -(sub(sn).n1(:,4) - sub(sn).n2(:,4))/2;
    c = sub(sn).a(:,5);
    c(c==2) = -1;
    RT = sub(sn).RT(:,5);
    RT(RT>RTmax) = nan;
    RT(RT<RTmin) = nan;
    i1 = sub(sn).gameLength == 5;
    i6 = sub(sn).gameLength == 10;
    
    B1(:, sn) = glmfit( [c(i1).*dR(i1) c(i1).*dI(i1)], (RT(i1)));
    B6(:, sn) = glmfit( [c(i6).*dR(i6) c(i6).*dI(i6)], (RT(i6)));
end


R = (rand(1,size(B1,2))-0.5)*0.25;
clear l1 l6
for i = 1:3
    axes(ax(i)); hold on;
    plot([1+R; 2+R], [B1(i,:); B6(i,:)], 'color', [1 1 1]*0.75, 'linewidth', 1)
    l1(i) = plot(1+R, B1(i,:),'o');
    l6(i) = plot(2+R, B6(i,:),'o');
    if i > 1
        plot([0.5 2.5], [0 0], 'k--', 'linewidth', 1)
    end
end
set(l1, 'color', AZblue)
set(l6, 'color', AZred)
set([l1 l6], 'markersize', 5, 'linewidth', 1)
set(ax(1), 'ylim', [0 2]);
set(ax(2), 'ylim', [-0.04 0.01]);
set(ax(3), 'ylim', [-0.2 0.2])
set(ax, 'xlim', [0.5 2.5], 'tickdir', 'out', 'fontsize', 16, ...
    'xtick', [1 2], 'xticklabel', [1 6])

axes(ax(1)); ylabel({'\beta_0 [seconds]'}); xlabel('horizon')
axes(ax(2)); ylabel({'\beta_R [seconds per point]'}); xlabel('horizon')
axes(ax(3)); ylabel({'\beta_I [seconds]'}); xlabel('horizon')

axes(ax(1));
plot([1 1 2 2], [1.8 1.9 1.9 1.8], 'k-', 'linewidth', 1)
text(1.5, 1.91, '*', 'fontsize', 30, 'horizontalalignment', 'center')

axes(ax(2));
plot([1 1 2 2], [0.005 0.008 0.008 0.005], 'k-', 'linewidth', 1)
text(1.5, 0.0082, '***', 'fontsize', 30, 'horizontalalignment', 'center')
