function e = plot_choiceCurvesFak_v2(ax, sub, binEdges, RTmin, RTmax)

global AZred AZblue

% binEdges = [-25:10:25];
for sn = 1:length(sub)
    RT = sub(sn).rt;
    ind = (RT>RTmin) & (RT<RTmax);
    dR = sub(sn).dR(ind);
    
    A = sub(sn).choice(ind);
    i1 = sub(sn).gameLength(ind) == 5;
    i6 = sub(sn).gameLength(ind) == 10;
    dI = -sub(sn).dI(ind);
    i22 = sub(sn).dI(ind) == 0;
    i13 = ~i22;
    
    uID(dI>0) = 1;
    uID(dI<0) = 2;
    uID(dI==0) = nan;
    
    ind = i13&i1;
    [M_13_1(:,sn), ~, X] = binIt(-dI(ind).*dR(ind), A(ind)==uID(ind)', binEdges, 'std');
    
    ind = i13&i6;
    [M_13_6(:,sn), ~, X] = binIt(-dI(ind).*dR(ind), A(ind)==uID(ind)', binEdges, 'std');
    [M_22_1(:,sn), ~, X] = binIt(dR(i22&i1), A(i22&i1)==2, binEdges, 'std');
    [M_22_6(:,sn), ~, X] = binIt(dR(i22&i6), A(i22&i6)==2, binEdges, 'std');
end


m_13_1 = nanmean(M_13_1,2);
s_13_1 = nanstd(M_13_1,[],2)/sqrt(length(sub));
m_13_6 = nanmean(M_13_6,2);
s_13_6 = nanstd(M_13_6,[],2)/sqrt(length(sub));
m_22_1 = nanmean(M_22_1,2);
s_22_1 = nanstd(M_22_1,[],2)/sqrt(length(sub));
m_22_6 = nanmean(M_22_6,2);
s_22_6 = nanstd(M_22_6,[],2)/sqrt(length(sub));

% figure(1); clf;
% set(gcf, 'position', [811   575   600   300])
% ax = easy_gridOfEqualFigures([0.3 0.12], [0.12 0.12 0.05]);
axes(ax(1)); hold on;
e = errorbar(X, m_13_1, s_13_1);
e(2) = errorbar(X, m_13_6, s_13_6);
xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
ylabel({'p(high info)'})
title('unequal [1 3]', 'fontweight', 'normal')
% leg = legend(e([2 1]), {'horizon 6' 'horizon 1'}, 'location', 'northwest');
% set(leg, 'position', [ 0.3083    0.3150    0.1683    0.1317])
axes(ax(2)); hold on;
e(3) = errorbar(X, m_22_1, s_22_1);
e(4) = errorbar(X, m_22_6, s_22_6);
xlabel({'difference in mean reward' 'R(left) - R(right)'});
ylabel({'p(left)'})
title('equal [2 2]', 'fontweight', 'normal')


set(ax, 'xlim', [-35 35], 'ylim', [0 1], 'tickdir', 'out');%, 'fontsize', 20)
set(e, 'linewidth', 3, 'marker', '.', 'markersize', 50)
set(e([1 3]), 'color', AZblue)
set(e([2 4]), 'color', AZred)


