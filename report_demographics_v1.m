function report_demographics_v1(sub_all)

global AZred AZblue 


e1 = strcmp({sub_all.expt_name}, 'pilot-v1');
e2 = strcmp({sub_all.expt_name}, 'repeater-v1');

age = [sub_all.age];
gen = strcmp({sub_all.gender},'m');
nex = [sub_all.nExcluded];
ex = [sub_all.excluded];

figure(1); clf;
set(gcf, 'position', [811   556   400   270])
hold on;
i1 = (gen==1) & (ex == 0);
i2 = (gen==0) & (ex == 0);
i3 = (gen==1) & (ex == 1);
i4 = (gen==0) & (ex == 1);
l = plot(age(i1), nex(i1),'o');
l(2) = plot(age(i2), nex(i2),'o');
l(3) = plot(age(i3), nex(i3),'o');
l(4) = plot(age(i4), nex(i4),'o');

set(l, 'markersize', 5, 'linewidth', 1)
set(l([1]), 'color', AZblue)
set(l([3]), 'color', AZblue*0.5+0.5, 'marker', 'x')
set(l([2]), 'color', AZred)
set(l([4]), 'color', AZred*0.5+0.5, 'marker', 'x')
% set(gca, 'xscale', 'log', 'xlim', [15 50])
ylabel('number of trials')
xlabel('age')
legend({['male included (n = ' num2str(sum(i1)) ')']
    ['female included (n = ' num2str(sum(i2)) ')']
    ['male excluded (n = ' num2str(sum(i3)) ')']
    ['female excluded (n = ' num2str(sum(i4)) ')']
    }, ...
    'location', 'southeast')
xlim([15 50])
set(gca, 'xtick', [18 25:5:50], 'tickdir', 'out', 'fontsize', 16)


male1 = sum(gen(e1)==1);
female1 = sum(gen(e1)==0);
age1a  = min(age(e1));
age1b  = max(age(e1));
meanAge1 = mean(age(e1));

male2 = sum(gen(e2)==1);
female2 = sum(gen(e2)==0);
age2a  = min(age(e2));
age2b  = max(age(e2));
meanAge2 = nanmean(age(e2));

male3 = sum((gen==1) & (ex == 0));
female3 = sum((gen==0) & (ex == 0));
age3a = min(age((ex == 0)));
age3b = max(age((ex == 0)));
meanAge3 = nanmean(age((ex == 0)));
% str = ['Data used in this paper come from two previous published data sets: ' ...  
%     '30 participants (%d male, %d female, ages %d-%d, mean %.1f) from ' ...
%     'the original Horizon Task paper \cite{wilson2014humans} and an ' ...
%     'additional 30 participants (XXX male, XXX female, ages XXX-XXX, ' ...
%     'mean XXX) who made up additional young adults in ' ...
%     '\cite{somerville2017charting}.  Both data sets were acquired at ' ...
%     'Princeton University. In both cases participants gave informed ' ... 
%     'consent and the studies were approved by the Institutional Review ' ...
%     'Board at Princeton.']

disp(sprintf('Data used in this paper come from two previous published data sets: '))
disp(sprintf('30 participants (%d male, %d female, ages %d-%d, mean %.1f) from ', male1, female1, age1a, age1b, meanAge1))
disp('the original Horizon Task paper \cite{wilson2014humans} and an additional ')
disp(sprintf('30 participants (%d male, %d female, ages %d-%d, mean %.1f) who ', male2, female2, age2a, age2b, meanAge2))


% after exclusion
disp(sprintf('This left %d participants (%d male, %d female, ages %d-%d, mean %.1f) for the main analysis.', male3+female3, male3, female3, age3a, age3b, meanAge3));
