function e = plot_RTvsdR_v1(ax, sub, binEdges, RTmin, RTmax)

global AZred AZblue

for sn = 1:length(sub)
    RT = sub(sn).RT(:,5);
    RT( ~((RT>RTmin) & (RT<RTmax)) ) = nan;
    dR = sub(sn).o2(:,4) - sub(sn).o1(:,4);
    dI = (sub(sn).n2(:,4) - sub(sn).n1(:,4))/2;
    i22 = dI == 0;
    
    i13 = ~i22;
    i1 = sub(sn).gameLength == 5;
    i6 = sub(sn).gameLength == 10;
    
    
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

axes(ax(1)); hold on;
e = errorbar(X, m_13_1, s_13_1);
e(2) = errorbar(X, m_13_6, s_13_6);
xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
ylabel({'response time' '[seconds]'})
% title('unequal [1 3]', 'fontweight', 'normal')
% legend({'horizon 1' 'horizon 6'}, 'location', 'south')

axes(ax(2)); hold on;
e(3) = errorbar(X, m_22_1, s_22_1);
e(4) = errorbar(X, m_22_6, s_22_6);
xlabel({'difference in mean reward' 'R(left) - R(right)'});
ylabel({'response time' '[seconds]'})
% title('equal [2 2]', 'fontweight', 'normal')


set(ax(1:2), 'xlim', [-35 35], 'ylim', [0.35 1.05])
set(e, 'linewidth', 3, 'marker', '.', 'markersize', 30)
set(e([1 3]), 'color', AZblue)
set(e([2 4]), 'color', AZred)

