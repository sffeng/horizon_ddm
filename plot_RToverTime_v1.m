function [l, leg] = plot_RToverTime_v1(ax, sub, RTmin, RTmax)

global AZsand AZred AZblue

for sn = 1:length(sub)
    
    rtX = sub(sn).RT;
    rtX( ~((rtX>RTmin) & (rtX<RTmax))) = nan;
    
    %rt5 = sub(sn).RT(:,5);
    %idx = ((rt5>0.1) & (rt5<3));
    i22 = (sub(sn).n1(:,4) == 2);
    i13 = ~i22;
    
    i1 = sub(sn).gameLength == 5;
    i6 = sub(sn).gameLength == 10;
    
    ind = i13 & i1;% & idx;
    RT(:,1,sn) = nanmean(rtX(ind, :));
    ind = i13 & i6;% & idx;
    RT(:,2,sn) = nanmean(rtX(ind, :));
    ind = i22 & i1;% & idx;
    RT(:,3,sn) = nanmean(rtX(ind, :));
    ind = i22 & i6;% & idx;
    RT(:,4,sn) = nanmean(rtX(ind, :));
end

rt = nanmean(RT,3);
ss = nanstd(RT, [], 3)/ sqrt(length(sub));

axes(ax(1)); hold on;
f = fill([0.5 4.5 4.5 0.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', [1 1 1]*0.9, 'linestyle', 'none')
mix = 0.7;

f = fill([4.5 5.5 5.5 4.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', 0.1*AZsand+0.9, 'linestyle', 'none')

l = errorbar(rt(:,[2:-1:1]), ss(:,1:2));
ylabel({'response time' '[seconds]'})
xlabel('trial number')
leg = legend(l(2:-1:1), {'horizon 1' 'horizon 6'}, 'location', 'northeast');
t = title('unequal information')


axes(ax(2)); hold on;
f = fill([0.5 4.5 4.5 0.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', [1 1 1]*0.9, 'linestyle', 'none')

f = fill([4.5 5.5 5.5 4.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', 0.1*AZsand+0.9, 'linestyle', 'none')
% [1 1 1]*mix+(1-mix)*[1 1 0]
l(3:4) = errorbar(rt(:,[4:-1:3]), ss(:,1:2));
ylabel({'response time' '[seconds]'})
xlabel('trial number')
t(2) = title('equal information')

set(l, 'linewidth', 3, 'marker', '.', 'markersize', 30)

set(l([1 3]), 'color', AZred)
set(l([2 4]), 'color', AZblue)

set(ax(1:2), 'ylim', [0.3 0.9], 'xtick', [1:10], ...
    'xticklabel', {'i1' 'i2' 'i3' 'i4' 1:6}, ...
    'xlim', [0.5 10.5])