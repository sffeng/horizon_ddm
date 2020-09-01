function e = plot_rtCurves_sub_v1(ax, sub, RTmin, RTmax)

global AZred AZblue


binEdges = [-25:10:25];
for sn = 1:length(sub)
    RT = sub(sn).RT(:,5);
    RT(RT<RTmin) = nan;
    RT(RT>RTmax) = nan;
    dR = sub(sn).o2(:,4) - sub(sn).o1(:,4);
    i22 = sub(sn).n2(:,4) == 2;
    i13 = ~i22;
    i1 = sub(sn).gameLength == 5;
    i6 = sub(sn).gameLength == 10;
    dI = (sub(sn).n2(:,4) - sub(sn).n1(:,4))/2;
    
    [M_13_1(:,sn), ~, X] = binIt(dI(i13&i1).*dR(i13&i1), RT(i13&i1), binEdges, 'std');
    [M_13_6(:,sn), ~, X] = binIt(dI(i13&i6).*dR(i13&i6), RT(i13&i6), binEdges, 'std');
    [M_22_1(:,sn), ~, X] = binIt(dR(i22&i1), RT(i22&i1), binEdges, 'std');
    [M_22_6(:,sn), ~, X] = binIt(dR(i22&i6), RT(i22&i6), binEdges, 'std');
end


m_13_1 = nanmean(M_13_1,2);
s_13_1 = nanstd(M_13_1,[],2)/sqrt(length(sub));
m_13_6 = nanmean(M_13_6,2);
s_13_6 = nanstd(M_13_6,[],2)/sqrt(length(sub));
m_22_1 = nanmean(M_22_1,2);
s_22_1 = nanstd(M_22_1,[],2)/sqrt(length(sub));
m_22_6 = nanmean(M_22_6,2);
s_22_6 = nanstd(M_22_6,[],2)/sqrt(length(sub));

%figure(1); clf;
%set(gcf, 'position', [811   575   600   300])
% ax = easy_gridOfEqualFigures([0.25 0.08], [0.15 0.15 0.03]);
%ax = easy_gridOfEqualFigures([0.3 0.12], [0.12 0.15 0.05]);
axes(ax(1)); hold on;
e = errorbar(X, m_13_1, s_13_1);
e(2) = errorbar(X, m_13_6, s_13_6);
xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
ylabel({'RT [z-score]'})
title('unequal [1 3]', 'fontweight', 'normal')
%legend({'horizon 1' 'horizon 6'}, 'location', 'northeast')

axes(ax(2)); hold on;
e(3) = errorbar(X, m_22_1, s_22_1);
e(4) = errorbar(X, m_22_6, s_22_6);
xlabel({'difference in mean reward' 'R(left) - R(right)'});
ylabel({'RT [z-score]'})
title('equal [2 2]', 'fontweight', 'normal')


%set(ax, 'xlim', [-35 35], 'ylim', [-0.15 1.65], 'tickdir', 'out', 'fontsize', 20)
set(e, 'linewidth', 3, 'marker', '.', 'markersize', 50)
set(e([1 3]), 'color', AZblue)
set(e([2 4]), 'color', AZred)