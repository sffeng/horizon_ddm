%% ========================================================================
%% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %%
%% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %%
%% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %%
%% ========================================================================
%%
clear
maindir     = '~/Dropbox/TEMP/2019_HorizonDDM/';
fundir      = [maindir 'matlab/'];
datadir     = [maindir 'matlab/'];
savedir     = [fundir 'modelFits/'];
addpath(fundir);
cd(fundir);
defaultPlotParameters

%% ========================================================================
%% LOAD DATA %% LOAD DATA %% LOAD DATA %% LOAD DATA %% LOAD DATA %%
%% LOAD DATA %% LOAD DATA %% LOAD DATA %% LOAD DATA %% LOAD DATA %%
%% LOAD DATA %% LOAD DATA %% LOAD DATA %% LOAD DATA %% LOAD DATA %%
%% ========================================================================

%% load human data
sub = load_humanData_v1(datadir, 'allHorizonData_v2.csv', 'DDM_demographics.csv');

%% load Sam's MCMC parameters
fname1 = 'fittedparams.csv';
fname2 = 'fittedT0.csv';
fit_MCMC = load_MCMCParameters_v1(fname1, fname2);

%% align sub sIDs to sim sIDs
[sub, sub_all] = align_subToSim_v1(sub, fit_MCMC);

%% add fields to fit_MCMC for easier comparison with data
fit_MCMC = addFrom_subToSim_v1(sub, fit_MCMC);

%% report demographics
report_demographics_v1(sub_all)
% saveFigurePdf(gcf, '~/Desktop/DDM_demographics')

%% load fake data I sent to Sam
fid = fopen('DDM_fake_20200423.csv');
hdr = textscan(fid, '%s%s%s%s%s%s', 1, 'delimiter', ',');
dat = textscan(fid, '%f%f%f%f%f%f', 'delimiter', ',');
fclose(fid)
U = unique(dat{1});
for sn = 1:length(U)
    idx = find(dat{1} == U(sn));
    
    sub_fake(sn).sID = dat{1}(idx(1));
    sub_fake(sn).dR = dat{2}(idx);
    sub_fake(sn).dI = dat{3}(idx);
    sub_fake(sn).gameLength = dat{4}(idx);
    sub_fake(sn).choice = dat{5}(idx);
    sub_fake(sn).rt = dat{6}(idx);
end
% hdr_str = '%s,%s,%s,%s,%s,%s\n';
% num_str = '%d,%f,%d,%d,%d,%f\n';
% sub_fake = 
%% fit the full model to fake data I sent Sam
% takes about 36 seconds on my machine
names  = { 'c^\mu_0' 'c^\mu_R' 'c^\mu_I' 'c^\beta_0' 'c^\beta_R' 'c^\beta_I' 'c^\alpha_0' 'c^\alpha_R' 'c^\alpha_I'  'T_0'    };
% names  = { 'c^A_0' 'c^A_R' 'c^A_I' 'c^Z_0' 'c^Z_R' 'c^Z_I' 'c^X_0' 'c^X_R' 'c^X_I'  'T_0'    };
%           cA_0 | cA_R | cA_I | cZ_0 | cZ_R | cZ_I | cX_0 | cX_R | cX_I |  T0
M      = [    1      1      1      1      1      1      1      1      1      1     ];

X0_all = [    0      0.01   0      1      0      0      0      0      0      0.05  ];
LB_all = [   -1    -10    -10      0     -3     -3     -1     -1     -1      0     ];
UB_all = [    1     10     10     10      3      3      1      1      1      3     ];

RTmin = 0.1;
RTmax = 3;
tic 
fit_fake_MLE = fit_MLE_DDM_v2(sub_fake, M, RTmin, RTmax, X0_all, LB_all, UB_all);
toc
%%
figure(1); clf
set(gcf, 'position', [811   176   600   650])
hg = ones(5,1)*0.1;
wg = ones(4,1)*0.1;
wg(1) = 0.27;
wg(end) = 0.03;
hg(1) = 0.07;
hg(end) = 0.09;
[ax, hb, wb] = easy_gridOfEqualFigures(hg, wg);
[l, r, p] = plot_compareTwoParameterSets(ax, fit_fake_MLE, fit_MCMC, names, 'sim', 'fit');


%% FIGURE XXX - illustration of different drivers of random exploration
%     0.0327    0.0466
%     0.0629    0.0395
%    -0.1054    0.3381
%     0.9381    0.8682
%    -0.0002   -0.0002
%    -0.0149   -0.0059
%    -0.0383    0.0413
%     0.0172    0.0061
%     0.0304   -0.0679
%     0.0618    0.0683
r_vals = [-30:0.01:30];

sigma1 = 1/ (2 * sqrt(2) * 0.06 * 0.9);
sigma6 = 1/ (2 * sqrt(2) * 0.04 * 0.8);


% mean parameter values
z1 = 0.9;
a1 = 0.06;

rho = z1 / a1;

A1 = -0.1 / a1;
A6 = 0.4 / a1;

% just change drift, c^\mu_R - SNR
ft1.XXfit = theory_qualitativeRT_type1_v1(A1, A6, sigma1, sigma6, rho, 1);

% just change baseline threshold, c^\beta_0
ft2.XXfit = theory_qualitativeRT_type1_v1(A1, A6, sigma1, sigma6, rho, 2);

% just change baseline drift
% ratio of c^\beta_R to c^\mu_0
% rho = 0.0002 / 0.0327; %= z1/a1
rho = 0.01;
% just change baseline threshold, c^\beta_0
% a1 = sqrt((1/ (2 * sqrt(2) * sigma1 * rr )));
% z1 = a1 * rr;
% 
% a6 = 1/ (2 * sqrt(2) * sigma6 * z1);
% 
% ft3.XXfit = [
%     a1        a6
%     0         0
%    -0.1       0.4
%     0         0
%     z1        z1
%     0         0
%     0         0
%     0         0
%     0         0
%     0.05      0.05];
ft3.XXfit = theory_qualitativeRT_type2_v1(A1, A6, sigma1, sigma6, rho, 1)
ft4.XXfit = theory_qualitativeRT_type2_v1(A1, A6, sigma1, sigma6, rho, 2)


ft3.XXfit(end,:) = 0.5;
ft4.XXfit(end,:) = 0.5;

ft1 = theory_ERDT_v1(ft1, r_vals);
ft2 = theory_ERDT_v1(ft2, r_vals);
ft3 = theory_ERDT_v1(ft3, r_vals);
ft4 = theory_ERDT_v1(ft4, r_vals);


% start figure ------------------------------------
figure(1); clf;
set(gcf, 'position', [695   365   700   450])
hg = [0.12 0.15 0.3];
wg = [0.18 0.08 0.04 0.08 0.04 0.01];
[~,hb,wb,ax] = easy_gridOfEqualFigures(hg, wg);
for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        axes(ax(i,j)); hold on;
    end
end
clear l
axes(ax(1,1)); l(:,1) = plot(ft1.xx, ft1.cc13_ana); 
axes(ax(1,2)); l(:,2) = plot(ft1.xx, ft1.rt13_ana);
axes(ax(2,1)); l(:,3) = plot(ft1.xx, ft1.cc22_ana); 
axes(ax(2,2)); l(:,4) = plot(ft1.xx, ft1.rt22_ana);

axes(ax(1,3)); l(:,5) = plot(ft1.xx, ft2.rt13_ana);
axes(ax(2,3)); l(:,6) = plot(ft1.xx, ft2.rt22_ana);

axes(ax(1,4)); l(:,7) = plot(ft1.xx, ft3.rt13_ana);
axes(ax(2,4)); l(:,8) = plot(ft1.xx, ft3.rt22_ana);

axes(ax(1,5)); l(:,9) = plot(ft1.xx, ft4.rt13_ana);
axes(ax(2,5)); l(:,10) = plot(ft1.xx, ft4.rt22_ana);

set(l(1,:), 'color', AZblue)
set(l(2,:), 'color', AZred)

set(ax(:,2:end), 'ylim', [0 1])
set(ax, 'xlim', [-35 35], 'tickdir', 'out')

delta = 0.07;
[a] = add_externalYlabels_v1(wg, hg, wb, hb, delta, ...
    { 'equal condition [2 2]' 'unequal condition [1 3]'})

delta2 = 0.25;
wg3 = [wg(1) wg(2) wg(end)];
wb3 = [wb(1) sum(wb(2:5))+sum(wg(3:5))];
[b] = add_externalXlabels(wg3, hg, wb3, hb, delta2, ...
    {'' 'response times' })


% delta2 = 0.25;
% wg2 = [wg(1) wg(2) wg(4) wg(end)];
% wb2 = [wb(1) wb(2)+wb(3)+wg(3) wb(4)+wb(5)+wg(5)];
% [b] = add_externalXlabels(wg2, hg, wb2, hb, delta2, ...
%     {'' 
%     'response times' 
%     'response times'})


delta2 = 0.13;
wg2 = [wg(1) wg(2) wg(4) wg(end)];
wb2 = [wb(1) wb(2)+wb(3)+wg(3) wb(4)+wb(5)+wg(5)];
[b] = add_externalXlabels(wg2, hg, wb2, hb, delta2, ...
    {'' 
    'c^\beta_R = c^\beta_I = 0' 
    'c^\mu_R = c^\mu_I = 0'})


delta3 = -0.12;
[b] = add_externalXlabels(wg, hg, wb, hb, delta3, ...
    {'choice curves' 
    'drift change, c^\mu_R' 
    'threshold change, c^\beta_0'
    'drift change, c^\mu_0'
    'threshold change, c^\beta_R'})
% set(b, 'linestyle','-')

delta4 = 0;
annotation('textbox', 'position', [wg(1) hb(1)+hg(1) sum(wg(2:5))+sum(wb) hg(2)-delta4] ,...
    'string', 'difference in mean reward, R(high info) - R(low info)', ...
    'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
    'fontsize', 16, ...
    'linestyle', 'none')

delta4 = 0.02;
annotation('textbox', 'position', [wg(1) 0 sum(wg(2:5))+sum(wb) hg(1)-delta4] ,...
    'string', 'difference in mean reward, R(left) - R(right)', ...
    'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
    'fontsize', 16, ...
    'linestyle', 'none')

axes(ax(1,1)); ylabel('p(high info)')
axes(ax(2,1)); ylabel('p(left)')

axes(ax(1,2)); ylabel('RT [seconds]')
axes(ax(2,2)); ylabel('RT [seconds]')

axes(ax(1,4)); ylabel('RT [seconds]')
axes(ax(2,4)); ylabel('RT [seconds]')

set(ax(:,[3 5]), 'yticklabel', [])
% 
addABCs(ax(1, :), [-0.02 0.18], 22)
saveFigurePdf(gcf, '~/Desktop/DDM_qualitative')
% addABCs(ax(1, [1 2 4]), [-0.04 0.25], 24, 'ABC')
% wg4 = [wg(1) wg(end)];
% wb4 = [sum(wb)+sum(wg(2:end-1))];
% hg4 = [hg(1) sum(hg(2:end))];
% hb4 = [hb(1)];
% delta4 = -0.2;
% [b] = add_externalXlabels(wg4, hg4, wb4, hb4, delta4, ...
%     {'11111' 
%     })

% set(b, 'linestyle', '-')

%% ========================================================================
%% theta-sigma model %% theta-sigma model %% theta-sigma model %%
%% theta-sigma model %% theta-sigma model %% theta-sigma model %%
%% theta-sigma model %% theta-sigma model %% theta-sigma model %%
%% ========================================================================

%% MLE fits
priorFlag = 1;
for sn = 1:length(sub)
    fit_logistic(sn) = fit_biasNoiseBonus_v2(sub(sn), priorFlag)
end

%% compute model choice curves
x_vals = [-30:0.01:30];
for sn = 1:length(sub)
    for i = 1:2
        cc13(:,i,sn) = fit(sn).cc13{i}(x_vals);
        cc22(:,i,sn) = fit(sn).cc22{i}(x_vals);
    end
end

CC13 = nanmean(cc13,3);
CC22 = nanmean(cc22,3);

%% basic choice curves

clear M_13_1 M_13_6 M_22_1 M_22_6
binEdges = [-25:10:25];
for sn = 1:length(sub)
    RT = sub(sn).RTz(:,5);
    dR = sub(sn).o2(:,4) - sub(sn).o1(:,4);
    
    A = sub(sn).a(:,5);
    i22 = sub(sn).n2(:,4) == 2;
    i13 = ~i22;
    i1 = sub(sn).gameLength == 5;
    i6 = sub(sn).gameLength == 10;
    dI = (sub(sn).n2(:,4) - sub(sn).n1(:,4))/2;
    uID(dI>0) = 1;
    uID(dI<0) = 2;
    uID(dI==0) = nan;
    
    ind = i13&i1&(~isnan(A));
    [M_13_1(:,sn), ~, X] = binIt(-dI(ind).*dR(ind), A(ind)==uID(ind)', binEdges, 'std');
    
    ind = i13&i6&(~isnan(A));
    [M_13_6(:,sn), ~, X] = binIt(-dI(ind).*dR(ind), A(ind)==uID(ind)', binEdges, 'std');
    
    ind = i22&i1&(~isnan(A));
    [M_22_1(:,sn), ~, X] = binIt(dR(ind), A(ind)==2, binEdges, 'std');
    ind = i22&i6&(~isnan(A));
    [M_22_6(:,sn), ~, X] = binIt(dR(ind), A(ind)==2, binEdges, 'std');
end


m_13_1 = nanmean(M_13_1,2);
s_13_1 = nanstd(M_13_1,[],2)/sqrt(length(sub));
m_13_6 = nanmean(M_13_6,2);
s_13_6 = nanstd(M_13_6,[],2)/sqrt(length(sub));
m_22_1 = nanmean(M_22_1,2);
s_22_1 = nanstd(M_22_1,[],2)/sqrt(length(sub));
m_22_6 = nanmean(M_22_6,2);
s_22_6 = nanstd(M_22_6,[],2)/sqrt(length(sub));

figure(1); clf;
set(gcf, 'position', [811   575   600   300])
ax = easy_gridOfEqualFigures([0.3 0.12], [0.12 0.12 0.05]);
% HERE
% e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
% set(e1, 'color', 'k')
% TO HERE
axes(ax(1)); hold on;
e = errorbar(X, m_13_1, s_13_1);
e(2) = errorbar(X, m_13_6, s_13_6);

l1 = plot(x_vals, CC13);

xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
ylabel({'p(high info)'})
t = title('unequal [1 3]', 'fontweight', 'normal');
leg = legend(e([2 1]), {'horizon 6' 'horizon 1'}, 'location', 'northwest');
set(leg, 'position', [ 0.3083    0.3150    0.1683    0.1317])
axes(ax(2)); hold on;
e(3) = errorbar(X, m_22_1, s_22_1);
e(4) = errorbar(X, m_22_6, s_22_6);

l2 = plot(x_vals, CC22)

set([l1(1) l2(1)], 'color', AZblue)
set([l1(2) l2(2)], 'color', AZred)

set(e, 'linestyle', 'none')

xlabel({'difference in mean reward' 'R(left) - R(right)'});
ylabel({'p(left)'})
t(2) = title('equal [2 2]', 'fontweight', 'normal');



set(ax, 'xlim', [-35 35], 'ylim', [0 1], 'tickdir', 'out', 'fontsize', 18)
set(t, 'fontsize', 20)
set(e, 'linewidth', 3, 'marker', '.', 'markersize', 40)
set(e([1 3]), 'color', AZblue)
set(e([2 4]), 'color', AZred)

saveFigurePdf(gcf, '~/Desktop/DDM_choices')

%% FIGURE 3 - basic choice curves + model parameters
XX = cat(3, fit.x);

bonus = squeeze(XX(3,:,:));
noise13 = squeeze(XX(2,:,:));
noise22 = squeeze(XX(5,:,:));

binEdges = [-25:10:25];

figure(1); clf;
set(gcf, 'position', [811   575   550   500])
ax = easy_gridOfEqualFigures([0.66  0.07], [0.14 0.14 0.07]);
ax(3:5) = easy_gridOfEqualFigures([0.1  0.6], [0.12 0.12 0.12 0.05]);
% choice curve bit --------------------------------------------------------
e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
e2 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
set([e1([2 4])], 'visible', 'off')
set([e2([1 3])], 'visible', 'off')

axes(ax(1)); hold on;

% e = errorbar(X, m_13_1, s_13_1);
% e(2) = errorbar(X, m_13_6, s_13_6);

l1 = plot(x_vals, CC13);

xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
ylabel({'p(high info)'})
t = title('unequal [1 3]', 'fontweight', 'normal');
leg = legend([e2(2) e1(1)], {'horizon 6' 'horizon 1'}, 'location', 'southeast');
set(leg, 'position', [ 0.3000    0.6706    0.1827    0.0790])

axes(ax(2)); hold on;
% e(3) = errorbar(X, m_22_1, s_22_1);
% e(4) = errorbar(X, m_22_6, s_22_6);

l2 = plot(x_vals, CC22)

set([l1(1) l2(1)], 'color', AZblue)
set([l1(2) l2(2)], 'color', AZred)

set([e1 e2], 'linestyle', 'none')

xlabel({'difference in mean reward' 'R(left) - R(right)'});
ylabel({'p(left)'})
t(2) = title('equal [2 2]', 'fontweight', 'normal');

set([e1 e2], 'linewidth', 3, 'marker', '.', 'markersize', 25)
set(e1([1 3]), 'color', AZblue)
set(e2([2 4]), 'color', AZred)
set(ax(1:2), 'xlim', [-35 35], 'ylim', [0 1], 'fontsize', 16)


% parameter values bit ----------------------------------------------------

axes(ax(3)); hold on;
R = repmat([1; 2], [1 size(bonus,2)])+0.25*(rand(size(bonus))-0.5);
plot(R, bonus, 'linewidth', 1, 'marker', 'none', 'color', [1 1 1]*0.75, ...
    'markersize', 10)
plot(R(1,:), bonus(1,:), 'o', 'markersize', 5, 'color', AZblue, 'linewidth', 1)
plot(R(2,:), bonus(2,:), 'o', 'markersize', 5, 'color', AZred, 'linewidth', 1)
t(3) = title('directed [1 3]');
ylabel('information bonus')
xlabel('horizon')
ylim([-40 60])

plot([1 1 2 2], [15 50 50 45], 'k-', 'linewidth', 1)
[~, p, ~, stats] = ttest(diff(bonus));
str = '';
if p < 0.01
    str = '*';
    if p < 0.005
        str = '**';
        if p < 0.001
            str = '***';
        end
    end
end
text(1.5, 52, str, 'fontsize', 30, 'horizontalalignment', 'center')

disp(sprintf('Info bonus: t(%d) = %.2f, p = %.2e', stats.df, stats.tstat, p));


axes(ax(4)); hold on;
plot(R, noise13, 'linewidth', 1, 'marker', 'none', 'color', [1 1 1]*0.75, ...
    'markersize', 10)
plot(R(1,:), noise13(1,:), 'o', 'markersize', 5, 'color', AZblue, 'linewidth', 1)
plot(R(2,:), noise13(2,:), 'o', 'markersize', 5, 'color', AZred, 'linewidth', 1)
plot([1 1 2 2], [33 45 45 35], 'k-', 'linewidth', 1)
[~, p, ~, stats] = ttest(diff(noise13));
str = '';
if p < 0.01
    str = '*';
    if p < 0.005
        str = '**';
        if p < 0.001
            str = '***';
        end
    end
end
text(1.5, 46, str, 'fontsize', 30, 'horizontalalignment', 'center')
disp(sprintf('Decision noise [1 3]: t(%d) = %.2f, p = %.2e', stats.df, stats.tstat, p));



t(4) = title('random [1 3]');
ylabel('decision noise')
xlabel('horizon')

axes(ax(5)); hold on;
plot(R, noise22, 'linewidth', 1, 'marker', 'none', 'color', [1 1 1]*0.75, ...
    'markersize', 10)
plot(R(1,:), noise22(1,:), 'o', 'markersize', 5, 'color', AZblue, 'linewidth', 1)
plot(R(2,:), noise22(2,:), 'o', 'markersize', 5, 'color', AZred, 'linewidth', 1)
plot([1 1 2 2], [32 45 45 40], 'k-', 'linewidth', 1)
[~, p, ~, stats] = ttest(diff(noise22));
str = '';
if p < 0.01
    str = '*';
    if p < 0.005
        str = '**';
        if p < 0.001
            str = '***';
        end
    end
end
text(1.5, 46, str, 'fontsize', 30, 'horizontalalignment', 'center')
disp(sprintf('Decision noise [2 2]: t(%d) = %.2f, p = %.2e', stats.df, stats.tstat, p));

t(5) = title('random [2 2]');
ylabel('decision noise')
xlabel('horizon')

% set(ax(4), 'yscale', 'log')
set(ax(3:5), 'xlim', [0.5 2.5], 'xticklabel', [1 6])
set(ax(4:5), 'ylim', [0 50])

set(ax,  'tickdir', 'out', 'fontsize', 16)

set(t, 'fontsize', 18, 'fontweight', 'normal')
set(t(1:5), 'units', 'normalized')
set(t(1:5), 'position',  [0.5000 1.07 0])
addABCs(ax(1:2), [-0.06 0.07], 24)
addABCs(ax(3:5), [-0.06 0.09], 24, 'CDE')

% saveFigurePdf(gcf, '~/Desktop/DDM_choices')



%% ========================================================================
%% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %%
%% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %%
%% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %%
%% ========================================================================

%% fit the full model
% takes about 36 seconds on my machine
names  = { 'c^\mu_0' 'c^\mu_R' 'c^\mu_I' 'c^\beta_0' 'c^\beta_R' 'c^\beta_I' 'c^\alpha_0' 'c^\alpha_R' 'c^\alpha_I'  'T_0'    };
% names  = { 'c^A_0' 'c^A_R' 'c^A_I' 'c^Z_0' 'c^Z_R' 'c^Z_I' 'c^X_0' 'c^X_R' 'c^X_I'  'T_0'    };
%           cA_0 | cA_R | cA_I | cZ_0 | cZ_R | cZ_I | cX_0 | cX_R | cX_I |  T0
M      = [    1      1      1      1      1      1      1      1      1      1     ];

X0_all = [    0      0.01   0      1      0      0      0      0      0      0.05  ];
LB_all = [   -1    -10    -10      0     -3     -3     -1     -1     -1      0     ];
UB_all = [    1     10     10     10      3      3      1      1      1      3     ];

RTmin = 0.1;
RTmax = 3;
tic 
fit_MLE = fit_MLE_DDM_v2(sub, M, RTmin, RTmax, X0_all, LB_all, UB_all);
toc

%% fit the "best" model
names  = { 'c^\mu_0' 'c^\mu_R' 'c^\mu_I' 'c^\beta_0' 'c^\beta_R' 'c^\beta_I' 'c^\alpha_0' 'c^\alpha_R' 'c^\alpha_I'  'T_0'    };
% names  = { 'c^A_0' 'c^A_R' 'c^A_I' 'c^Z_0' 'c^Z_R' 'c^Z_I' 'c^X_0' 'c^X_R' 'c^X_I'  'T_0'    };
%           cA_0 | cA_R | cA_I | cZ_0 | cZ_R | cZ_I | cX_0 | cX_R | cX_I |  T0
M      = [    0      1      1      1      0      0      0      1      0      1     ];

X0_all = [    0      0.01   0      1      0      0      0      0      0      0.05  ];
LB_all = [   -1    -10    -10      0     -3     -3     -1     -1     -1      0     ];
UB_all = [    1     10     10     10      3      3      1      1      1      3     ];

RTmin = 0.1;
RTmax = 3;
tic 
fit_MLEbest = fit_MLE_DDM_v2(sub, M, RTmin, RTmax, X0_all, LB_all, UB_all);
toc


%% SUPPLEMENTARY FIGURE XXX - compare MLE fit to MCMC fit
% note that to make MLE line up with MCMC I have to do a few things to 
% Sam's MCMC parameters :
%   * flip signs on A1_dI, A6_dI, z1_dI, z6_dI, x1_dI, x6_dI
%   * double x1_0, x6_0

figure(1); clf
set(gcf, 'position', [811   176   600   650])
hg = ones(5,1)*0.1;
wg = ones(4,1)*0.1;
wg(1) = 0.27;
wg(end) = 0.03;
hg(1) = 0.07;
hg(end) = 0.09;
[ax, hb, wb] = easy_gridOfEqualFigures(hg, wg);
[l, r, p] = plot_compareTwoParameterSets(ax, fit_MLE, fit_MCMC, names, 'MLE', 'MCMC');

clear t
for i = 1:length(r)
    for j = 1:length(r{i})
        
        axes(ax(j));
        t(i,j) = text(1, 0.25-(i-1)*0.15, sprintf('r = %.2g', r{i}(j,j)), 'units', 'normalized', 'fontsize', 12, ...
            'horizontalalignment', 'right')
        
    end
end
        
set(t(1,:), 'color', AZblue)
set(t(2,:), 'color', AZred)
delta = 0.08;
delta2 = 0.04;
[a, b] = add_variableNames_v1(wg, hg, wb, hb ,delta, delta2)

saveFigurePdf(gcf, '~/Desktop/DDM_MLEvsMCMCparameters')

%% parameter recovery for MLE fit -----------------------------------------
%% simulate fake data with MLE parameters
% TODO: 
%   1. could scramble parameters 
%   2. could use different set of trials

clear fak
dt = 0.001;
for sn = 1:length(fit_MLE)
    fak(sn) = simulate_Mmodel_v1(fit_MLE(sn).XXfit, dt, fit_MLE(sn).dR, fit_MLE(sn).dI, fit_MLE(sn).gameLength);
end

%% fit fake data with MLE
fit_fak = fit_MLE_DDM_v2(fak, M, RTmin, RTmax, X0_all, LB_all, UB_all);

%% compare recovered parameters with original

figure(1); clf
set(gcf, 'position', [811   176   600   650])
hg = ones(5,1)*0.1;
wg = ones(4,1)*0.1;
wg(1) = 0.27;
wg(end) = 0.03;
hg(1) = 0.07;
hg(end) = 0.09;
[ax, hb, wb] = easy_gridOfEqualFigures(hg, wg);
[l, r, p] = plot_compareTwoParameterSets(ax, fit_MLE, fit_fak, names, 'sim', 'fit');

clear t
for i = 1:length(r)
    for j = 1:length(r{i})
        
        axes(ax(j));
        t(i,j) = text(1, 0.25-(i-1)*0.15, sprintf('r = %.2g', r{i}(j,j)), 'units', 'normalized', 'fontsize', 12, ...
            'horizontalalignment', 'right')
        
    end
end
        
set(t(1,:), 'color', AZblue)
set(t(2,:), 'color', AZred)
delta = 0.08;
delta2 = 0.04;
[a, b] = add_variableNames_v1(wg, hg, wb, hb ,delta, delta2)

saveFigurePdf(gcf, '~/Desktop/DDM_parameterRecoveryMLE')


%% ------------------------------------- end parameter recovery for MLE fit

%% posterior predictive checks --------------------------------------------
%% compute theoretical average behavior for fit parameter values
r_vals = [-30:0.01:30];
fit_MLE = theory_ERDT_v1(fit_MLE, r_vals);
fit_MLEbest = theory_ERDT_v1(fit_MLEbest, r_vals);

%% plot choice and RT model vs human
dum = cat(3, fit_MLE.cc13_ana);
cc13 = nanmean(dum,3);
dum = cat(3, fit_MLE.cc22_ana);
cc22 = nanmean(dum,3);

dum = cat(3, fit_MLE.rt13_ana);
rt13 = nanmean(dum,3);
dum = cat(3, fit_MLE.rt22_ana);
rt22 = nanmean(dum,3);

binEdges = [-25:10:25];

figure(1); clf;
ax = easy_gridOfEqualFigures([0.1 0.1 0.1], [0.1 0.1 0.03]);

axes(ax(1)); hold on; l = plot(r_vals, cc13);
axes(ax(2)); hold on; l(:,2) = plot(r_vals, cc22);
axes(ax(3)); hold on; l(:,3) = plot(r_vals, rt13);
axes(ax(4)); hold on; l(:,4) = plot(r_vals, rt22);

e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
set(e1, 'linestyle', 'none', 'markersize', 30)
e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1 , 3)
set(e1, 'linestyle', 'none', 'markersize', 30)

ff = 0.75;
set(l(1,:), 'color', AZblue*ff+(1-ff));%, 'visible', 'off')
set(l(2,:), 'color', AZred*ff+(1-ff));%
% set(e1([2 4]), 'visible', 'off')
% set(e1([1 3]), 'visible', 'off')

% set([e1([2 4]) ], 'visible', 'off')

%% posterior predictive for best fitting model
dum = cat(3, fit_MLEbest.cc13_ana);
cc13 = nanmean(dum,3);
dum = cat(3, fit_MLEbest.cc22_ana);
cc22 = nanmean(dum,3);

dum = cat(3, fit_MLEbest.rt13_ana);
rt13 = nanmean(dum,3);
dum = cat(3, fit_MLEbest.rt22_ana);
rt22 = nanmean(dum,3);

binEdges = [-25:10:25];

figure(1); clf;
ax = easy_gridOfEqualFigures([0.1 0.1 0.1], [0.1 0.1 0.03]);

axes(ax(1)); hold on; l = plot(r_vals, cc13);
axes(ax(2)); hold on; l(:,2) = plot(r_vals, cc22);
axes(ax(3)); hold on; l(:,3) = plot(r_vals, rt13);
axes(ax(4)); hold on; l(:,4) = plot(r_vals, rt22);

e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
set(e1, 'linestyle', 'none', 'markersize', 30)
e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1 , 3)
set(e1, 'linestyle', 'none', 'markersize', 30)

ff = 0.75;
set(l(1,:), 'color', AZblue*ff+(1-ff));%, 'visible', 'off')
set(l(2,:), 'color', AZred*ff+(1-ff));%
% set(e1([2 4]), 'visible', 'off')
% set(e1([1 3]), 'visible', 'off')

% set([e1([2 4]) ], 'visible', 'off')

%% posterior predictive for best fitting model - 4 x 2
dum = cat(3, fit_MLEbest.cc13_ana);
cc13 = nanmean(dum,3);
dum = cat(3, fit_MLEbest.cc22_ana);
cc22 = nanmean(dum,3);

dum = cat(3, fit_MLEbest.rt13_ana);
rt13 = nanmean(dum,3);
dum = cat(3, fit_MLEbest.rt22_ana);
rt22 = nanmean(dum,3);

binEdges = [-25:10:25];

figure(1); clf;
[~,~,~,ax] = easy_gridOfEqualFigures([0.1 0.1 0.1 0.1 0.1], [0.1 0.1 0.03]);
ax = ax(:);
axes(ax(1)); hold on; l = plot(r_vals, cc13(:,1));
axes(ax(2)); hold on; l(2,1) = plot(r_vals, cc13(:,2));
axes(ax(3)); hold on; l(1,2) = plot(r_vals, cc22(:,1));
axes(ax(4)); hold on; l(2,2) = plot(r_vals, cc22(:,2));
axes(ax(5)); hold on; l(1,3) = plot(r_vals, rt13(:,1));
axes(ax(6)); hold on; l(2,3) = plot(r_vals, rt13(:,2));
axes(ax(7)); hold on; l(1,4) = plot(r_vals, rt22(:,1));
axes(ax(8)); hold on; l(2,4) = plot(r_vals, rt22(:,2));

e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
set(e1([2 3 4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)

e1 = plot_choiceCurves_v2(ax(2:3), sub, binEdges, 0.1, 3);
set(e1([1  3 4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)

e1 = plot_choiceCurves_v2(ax(2:3), sub, binEdges, 0.1, 3);
set(e1([1 2  4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)

e1 = plot_choiceCurves_v2(ax(3:4), sub, binEdges, 0.1, 3);
set(e1([1 2 3 ]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)

e1 = plot_rtCurves_sub_v1(ax(5:6), sub, 0.1 , 3)
set(e1([2 3 4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)

e1 = plot_rtCurves_sub_v1(ax(6:7), sub, 0.1 , 3)
set(e1([1  3 4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)

e1 = plot_rtCurves_sub_v1(ax(6:7), sub, 0.1 , 3)
set(e1([1 2 4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)

e1 = plot_rtCurves_sub_v1(ax(7:8), sub, 0.1 , 3)
set(e1([1 2 3]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)

ff = 0.75;
set(l(1,:), 'color', AZblue*ff+(1-ff));%, 'visible', 'off')
set(l(2,:), 'color', AZred*ff+(1-ff));%

set(ax(5:end), 'ylim', [0.4 1])

% set(e1([2 4]), 'visible', 'off')
% set(e1([1 3]), 'visible', 'off')

% set([e1([2 4]) ], 'visible', 'off')

%% posterior predictive for best fitting model - 2 x 3
dum = cat(3, fit_MLEbest.cc13_ana);
cc13 = nanmean(dum,3);
dum = cat(3, fit_MLEbest.cc22_ana);
cc22 = nanmean(dum,3);

dum = cat(3, fit_MLEbest.rt13_ana);
rt13 = nanmean(dum,3);
dum = cat(3, fit_MLEbest.rt22_ana);
rt22 = nanmean(dum,3);

binEdges = [-25:10:25];

figure(1); clf;
set(gcf, 'position', [611   305   680   400])
hg = [0.12 0.16 0.11];
wg = [0.2 0.13 0.07 0.03];
[~,hb,wb,ax] = easy_gridOfEqualFigures(hg, wg);
ax = ax(:);
axes(ax(1)); hold on; l = plot(r_vals, cc13);
axes(ax(2)); hold on; l(:,2) = plot(r_vals, cc22);

axes(ax(3)); hold on; l(1,3) = plot(r_vals, rt13(:,1));
axes(ax(5)); hold on; l(2,3) = plot(r_vals, rt13(:,2));
axes(ax(4)); hold on; l(1,4) = plot(r_vals, rt22(:,1));
axes(ax(6)); hold on; l(2,4) = plot(r_vals, rt22(:,2));

e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
% set(e1([ 3 4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 20)
% %%
% e1 = plot_choiceCurves_v2(ax(2:3), sub, binEdges, 0.1, 3);
% set(e1([1  3 4]), 'visible', 'off')
% set(e1, 'linestyle', 'none', 'markersize', 30)
% 
% e1 = plot_choiceCurves_v2(ax(2:3), sub, binEdges, 0.1, 3);
% set(e1([1 2  4]), 'visible', 'off')
% set(e1, 'linestyle', 'none', 'markersize', 30)
% 
% e1 = plot_choiceCurves_v2(ax(3:4), sub, binEdges, 0.1, 3);
% set(e1([1 2 3 ]), 'visible', 'off')
% set(e1, 'linestyle', 'none', 'markersize', 30)

e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1 , 3)
set(e1([2 4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 20)

e1 = plot_rtCurves_sub_v1(ax(5:6), sub, 0.1 , 3)
set(e1([1 3]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 20)

ff = 0.75;
set(l(1,:), 'color', AZblue*ff+(1-ff));%, 'visible', 'off')
set(l(2,:), 'color', AZred*ff+(1-ff));%

set(ax(3:end), 'ylim', [0.3 1.1])

axes(ax(1)); leg = legend([l(2,1) l(1,1) ], {'horizon 6' 'horizon 1'}, 'location', 'southeast');
axes(ax(1)); title(''); xlabel('R(high info) - R(low info)')
axes(ax(2)); title({''}); xlabel('R(left) - R(right)')
axes(ax(3)); title({'horizon 1'}); ylabel('RT [seconds]'); xlabel('R(high info) - R(low info)')
axes(ax(4)); title({''}); ylabel('RT [seconds]'); xlabel('R(left) - R(right)')
axes(ax(5)); title({'horizon 6'}); ylabel(''); xlabel('R(high info) - R(low info)')
axes(ax(6)); title({''}); ylabel(''); xlabel('R(left) - R(right)')

set(leg, 'position', [0.2956    0.5974    0.1375    0.0862])

% set(ax(5:6), 'yticklabel', [])

set(ax, 'tickdir', 'out', 'xlim', [-35 35])
wg2 = [wg(1) wg(2) wg(4)];
wb2 = [wb(1) wb(2)+wb(3)+wg(3)];
delta2 = 0.06;
[b] = add_externalXlabels(wg2, hg, wb2, hb, delta2, {'choice' 'response time'})
set(b, 'linestyle', 'none')
delta = 0.075;
[a] = add_externalYlabels_v1(wg, hg, wb, hb, delta, {'equal condition [2 2]' 'unequal condition [1 3]'})
% set(e1([2 4]), 'visible', 'off')
% set(e1([1 3]), 'visible', 'off')
set([a b], 'fontsize', 18)
set(b, 'fontsize', 20)
set(ax, 'fontsize', 14)
% set([e1([2 4]) ], 'visible', 'off')
saveFigurePdf(gcf, '~/Desktop/DDM_posteriorPredictiveMLEbest')

%% compare fake data with human data
binEdges = [-25:10:25];
figure(1); clf;
set(gcf, 'position', [811   575   600   500])
ax = easy_gridOfEqualFigures([0.15 0.2 0.05], [0.12 0.15 0.05]);

% choices -----------------------------------------------------------------
e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
e2 = plot_choiceCurvesFak_v2(ax(1:2), fak, binEdges, 0.1, 3);
set(e1, 'linestyle', 'none', 'markersize', 30)
set(e2, 'marker', 'none')

set([e1([2 4]) e2([2 4])], 'visible', 'off')
e3 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
e4 = plot_choiceCurvesFak_v2(ax(1:2), fak, binEdges, 0.1, 3);
set([e3([1 3]) e4([1 3])], 'visible', 'off')
set(e3, 'linestyle', 'none', 'markersize', 30)
set(e4, 'marker', 'none')
leg = legend([e1([1]) e3(2) e2(1) e4(2)], ...
    {'human h = 1' 'human h = 6' 'model h = 1' 'model h = 6'}, ...
    'location', 'southeast');
set(leg, 'position', [0.3216    0.6614    0.1867    0.1310]);


% RTs ---------------------------------------------------------------------
e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1 , 3)
e2 = plot_rtCurves_subFak_v2(ax(3:4), fak, 0.1 , 3)
set(e1, 'linestyle', 'none', 'markersize', 30)
set(e2, 'marker', 'none')

set([e1([2 4]) e2([2 4])], 'visible', 'off')

e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1, 3)
e2 = plot_rtCurves_subFak_v2(ax(3:4), fak, 0.1 , 3)
set([e1([1 3]) e2([1 3])], 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)
set(e2, 'marker', 'none')

set(ax(3:4), 'ylim', [0.3 1.2])

axes(ax(1)); title('unequal [1 3], horizon 1')
axes(ax(2)); title('equal [2 2], horizon 1')
axes(ax(3)); title('unequal [1 3], horizon 6'); ylabel({'reaction time' '[seconds]'})
axes(ax(4)); title('equal [2 2], horizon 6'); ylabel({'reaction time' '[seconds]'})

%% -----------------------------------------end posterior predictive checks


%% ========================================================================
%% MLE MODEL COMPARISON %% MLE MODEL COMPARISON %% MLE MODEL COMPARISON %%
%% MLE MODEL COMPARISON %% MLE MODEL COMPARISON %% MLE MODEL COMPARISON %%
%% MLE MODEL COMPARISON %% MLE MODEL COMPARISON %% MLE MODEL COMPARISON %%
%% ========================================================================

%% fit and compare ALL models by BIC and AIC
% you have to have a baseline threshold (and we'll also assume a T0) to be a DDM in 
% the first place
% one run takes about 36 seconds on my machine
% 64 models should take around 38 minutes
% actually took 1237.235324 seconds = 20.6 minutes
% so 256 should take around 80 minutes
% actually took 3823.023474 seconds = 63.71 minutes
names  = { 'c^\mu_0' 'c^\mu_R' 'c^\mu_I' 'c^\beta_0' 'c^\beta_R' 'c^\beta_I' 'c^\alpha_0' 'c^\alpha_R' 'c^\alpha_I'  'T_0'    };

count = 1;
for i = 0:1
    for j = 0:1
        for k = 0:1
            for l = 0:1
                for m = 0:1
                    for n = 0:1
                        for o = 0:1
                            for p = 0:1
                                %               cA_0 | cA_R | cA_I | cZ_0 | cZ_R | cZ_I | cX_0 | cX_R | cX_I |  T0
                                M(count,:) = [    i      j      k      1      l      m      n      o      p      1     ];
                                count = count + 1;
                            end
                        end
                    end
                end
            end
        end
    end
end
X0_all     = [    0      0.01   0      1      0      0      0      0      0      0.05  ];
LB_all     = [   -1    -10    -10      0     -3     -3     -1     -1     -1      0     ];
UB_all     = [    1     10     10     10      3      3      1      1      1      3     ];

    
RTmin = 0.1;
RTmax = 3;

tic
for count = 1:size(M,1)
    count
    fit_MLE256(count,:) = fit_MLE_DDM_v2(sub, M(count,:), RTmin, RTmax, X0_all, LB_all, UB_all);
end
toc


% save ~/Desktop/DDM_fits_20200423_256
%%



%% fit and compare a bunch of models by BIC and AIC
% one run takes about 36 seconds on my machine
% 64 models should take around 38 minutes
% actually took 1237.235324 seconds = 20.6 minutes
names  = { 'c^\mu_0' 'c^\mu_R' 'c^\mu_I' 'c^\beta_0' 'c^\beta_R' 'c^\beta_I' 'c^\alpha_0' 'c^\alpha_R' 'c^\alpha_I'  'T_0'    };

count = 1;
for i = 0:1
    for j = 0:1
        for k = 0:1
            for l = 0:1
                for m = 0:1
                    for n = 0:1
                        %               cA_0 | cA_R | cA_I | cZ_0 | cZ_R | cZ_I | cX_0 | cX_R | cX_I |  T0
                        M(count,:) = [    i      1      1      1      j      k      l      m      n      1     ];
                        count = count + 1;
                    end
                end
            end
        end
    end
end
X0_all     = [    0      0.01   0      1      0      0      0      0      0      0.05  ];
LB_all     = [   -1    -10    -10      0     -3     -3     -1     -1     -1      0     ];
UB_all     = [    1     10     10     10      3      3      1      1      1      3     ];

    
RTmin = 0.1;
RTmax = 3;

tic
for count = 1:size(M,1)
    count
    fit_MLE64(count,:) = fit_MLE_DDM_v2(sub, M(count,:), RTmin, RTmax, X0_all, LB_all, UB_all);
end
toc

%%
% save ~/Desktop/DDM_fits_20200423
%%

% gonna need to compute BIC_all which is not quite BIC1 + BIC6 because
% log(n1) + log(n2) \ne log(n1 + n2)

% for now this is an approximation
for i = 1:size(fit_MLE64,1)
    for j = 1:size(fit_MLE64,2)
        BIC(i,j) = sum(fit_MLE64(i,j).BIC);
        AIC(i,j) = sum(fit_MLE64(i,j).AIC);
    end
end
%%
for i = 1:size(fit_MLE256,1)
    for j = 1:size(fit_MLE256,2)
        BIC(i,j) = sum(fit_MLE256(i,j).BIC);
        AIC(i,j) = sum(fit_MLE256(i,j).AIC);
    end
end
%%
[nBestFit_BIC, fBestFit_BIC] = compute_nBestFit(BIC');
[nBestFit_AIC, fBestFit_AIC] = compute_nBestFit(AIC');

%%
% might need to run in Matlab 2016 to get this to work
% [alpha, exp_r_bic, xp_bic] = spm_BMS(BIC);
% [alpha, exp_r_aic, xp_aic] = spm_BMS(AIC);

%%
[sBIC,ind] = sort(fBestFit_BIC, 'descend');%   mean(BIC, 2));
% [sBIC,ind] = sort(nBestFit_BIC-mean(BIC, 2)/max(mean(BIC,2)), 'descend');
idx = find(sBIC>0);
figure(1); clf;
set(gcf, 'position', [547   102   450   380])
ax = easy_gridOfEqualFigures([0.03 0.2], [0.14 0.25]);
ax(2) = easy_gridOfEqualFigures([0.03 0.2], [0.8 0.05]);
axes(ax(1));  
imagesc(M(ind(1:length(idx)),:))
hold on;
[l1,l2] = addFacetLines(M(ind(1:length(idx)),:));
set([l1, l2], 'linewidth', 1)
yl = get(gca, 'ylim');

set(ax(1), 'xaxislocation', 'top', ...
    'xtick', 1:10, ...
    'xticklabel', names, 'ytick', [1:10])
xlab = xlabel('parameter');
ylab = ylabel('model number (ranked by BIC)');
cc = [1 1 1
    (AZred+AZblue)/2
    AZblue*0.7+0.3];
colormap(cc);

axes(ax(2)); hold on;
b = bar(nBestFit_BIC(ind(1:length(idx))));
set(b, 'facecolor', (AZred), 'barwidth', 0.8)
set(ax(2), 'view', [90 90],  'xlim', yl, 'yaxislocation', 'right', ...
    'xtick', [1:10], 'xticklabel', [])
ylab(2) = ylabel({'number' 'best fit'});



set(ax, 'tickdir', 'out', 'fontsize', 14)
set([xlab ylab], 'fontsize', 18)

saveFigurePdf(gcf, '~/Desktop/DMM_modelComparison_MLE')


%% THIS ONE??? maybe
% MISSING T0!!!!
[sBIC,ind] = sort(mean(BIC, 2));


figure(1); clf;
set(gcf, 'position', [547   102   450   750])
ax = easy_gridOfEqualFigures([0.01 0.09], [0.15 0.01 0.01 0.05 0.07]);
axes(ax(1));  
imagesc(M(ind,1:3))
hold on;
[l1,l2] = addFacetLines(M(:,1:3));
set([l1, l2], 'linewidth', 1)
yl = get(gca, 'ylim');

set(ax(1), 'xaxislocation', 'top', ...
    'xtick', 1:3, ...
    'xticklabel', names(1:3))
t = title('drift, \mu');
ylab = ylabel('model number (ranked by total BIC)')
axes(ax(2));  
imagesc(M(ind,4:6))
hold on;
[l1,l2] = addFacetLines(M(:,1:3));
set([l1, l2], 'linewidth', 1)
yl = get(gca, 'ylim');

set(ax(2), 'xaxislocation', 'top', ...
    'xtick', 1:3, ...
    'xticklabel', names(4:6))
t(2) = title('threshold, \beta');


axes(ax(3));  
imagesc(M(ind,7:9))
hold on;
[l1,l2] = addFacetLines(M(:,1:3));
set([l1, l2], 'linewidth', 1)
yl = get(gca, 'ylim');

set(ax(3), 'xaxislocation', 'top', ...
    'xtick', 1:3, ...
    'xticklabel', names(7:9))

t(3) = title('bias, \alpha');
set(ax(1:3), 'box', 'off', 'ytick', [1 10:10:60 64])
set(ax(2:3), 'ytick', [])

axes(ax(end)); hold on;
b = bar(fBestFit_BIC(ind));
set(b, 'facecolor', (AZred+AZblue)/2, 'barwidth', 1)
set(gca, 'view', [90 90], ...
    'xlim', yl, ...
    'yaxislocation', 'right')
ylab(2) = ylabel({'fraction' 'best fit'});
set(ax(end), 'xticklabel', [])

set(ax, 'tickdir', 'out', 'fontsize', 14)
set([t ylab], 'fontweight', 'normal', 'fontsize', 16)

set(ax, 'tickdir', 'out')

saveFigurePdf(gcf, '~/Desktop/DDM_MLE_modelComparison')

% axes(ax(3)); 
% b = plot(sBIC);
% % set(b, 'facecolor', AZred, 'barwidth', 1)
% set(ax(3), 'view', [90 90])
% set(ax(3), 'xlim', yl)




%% ========================================================================
%% MCMC FIT %% MCMC FIT %% MCMC FIT %% MCMC FIT %% MCMC FIT %% MCMC FIT %%
%% MCMC FIT %% MCMC FIT %% MCMC FIT %% MCMC FIT %% MCMC FIT %% MCMC FIT %%
%% MCMC FIT %% MCMC FIT %% MCMC FIT %% MCMC FIT %% MCMC FIT %% MCMC FIT %%
%% ========================================================================

%% simulate fake data with MCMC parameters
clear fak
dt = 0.001;
for sn = 1:length(fit_MLE)
    fak_MCMC(sn) = simulate_Mmodel_v1(fit_MCMC(sn).XXfit, dt, fit_MLE(sn).dR, fit_MLE(sn).dI, fit_MLE(sn).gameLength);
    fak_MCMC(sn).sID = ones(size(fak_MCMC(sn).sID))*fit_MCMC(sn).sID;
end

%% write fake data out to CSV
fname = '~/Desktop/DDM_fake_20200423.csv';

n = fieldnames(fak_MCMC(1));
n = {n{[1 3:end]}};

hdr_str = '%s,%s,%s,%s,%s,%s\n';
num_str = '%d,%f,%d,%d,%d,%f\n';

fid = fopen(fname, 'w');
fprintf(fid, hdr_str, n{:});

for sn = 1:length(fak_MCMC)
    clear X
    for i = 1:length(n)
        X(i,:) = getfield(fak_MCMC(sn), n{i});
    end
    
    fprintf(fid, num_str, X);
end
fclose(fid);

%% sensitivity analysis - how much does signal vs threshold contribute to
%% random exploration

XX = cat(3,fit_MCMC.XXfit);

iZ_0 = 4;
iA_dR = 2;

x1 = squeeze(XX(iA_dR,:,:));
x2 = squeeze(XX(iZ_0,:,:));

rA = x1(1,:) ./ x1(2,:);
rZ = x2(1,:) ./ x2(2,:);

%% compare with logistic noise change
X = cat(3,fit_logistic.x);
i13 = 2;
i22 = 5;

n13 = squeeze(X(i13, :, :));
n22 = squeeze(X(i22, :, :));

r13 = n13(2,:) ./ n13(1,:);
r22 = n22(2,:) ./ n22(1,:);
% r13(r13 > 100) = nan;
% r22(r22 > 100) = nan;


figure(1); clf; 
set(gcf, 'position', [814   511   550   300])
ax = easy_gridOfEqualFigures([0.13 0.1], [0.08 0.1 0.03]);
rnd = (rand(size(rA))-0.5)/8;

axes(ax(1)); hold on
ff = 0.5;
l = plot(1+rnd, r13, 'o' ,'color', AZcactus*ff+(1-ff));
l(2) = plot(2+rnd, r22, 'o' ,'color', AZsand*ff+(1-ff));
b = boxplot([r13' r22' ], 'notch', 'on', 'widths', 0.25);
set(b(7,:), 'visible', 'off')
set(b, 'color', 'k', 'linewidth', 2)
plot([0.5 2.5], [1 1], 'k--', 'linewidth', 1)
t = title({'change in noise' 'from logistic model'});
% ylabel('ratio \sigma(horizon 6) / \sigma(horizon 1)')
ylabel({'\sigma ratio' 'horizon 6 : horizon 1'})
xlabel('uncertainty condition')

axes(ax(2)); hold on;
l(3) = plot(1 + rnd, rA, 'o', 'color', AZmesa*ff+(1-ff));
l(4) = plot(2+rnd, rZ, 'o' ,'color', AZsky*ff+(1-ff));
plot([0.5 2.5], [1 1], 'k--', 'linewidth', 1)
set(l, 'markersize', 5, 'linewidth', 1)
b = boxplot([rA' rZ' ], 'notch', 'on', 'widths', 0.25);
set(b(7,:), 'visible', 'off')
set(b, 'color', 'k', 'linewidth', 2)
t(2) = title({'change in' 'DDM parameters'});
set(t, 'fontweight', 'normal')
xlabel('parameter')
ylabel({'parameter ratio' 'horizon 1 : horizon 6 '})
plot([1 1 2 2], [4.5 6 6 2], 'k-', 'linewidth', 1)
text(1.5, 6.1, '***', 'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
    'fontsize', 28)

set(ax(1), 'xtick', [1 2], 'xticklabel', {'[1 3]' '[2 2]'})
set(ax(2), 'xtick', [1 2], 'xticklabel', {'c^\mu_R' 'c^\beta_0'})
set(ax, 'ylim', [0 10], 'box', 'off', 'tickdir', 'out')
ax(2).TickLabelInterpreter='tex';
set(ax, 'ytick', [1:2:9])

addABCs(ax, [-0.08 0.1], 28)
saveFigurePdf(gcf, '~/Desktop/DDM_sigmaVsDDMpars')
% set(gca, 'xscale', 'log', 'yscale', 'log')




%% compute parameters on every trial for sim
for sn = 1:length(sim)
    for i = 1:length(sim(sn).dR)
        if sim(sn).h(i) == 1
            
            sim(sn).z(i,1)    = sim(sn).z1_0 + sim(sn).z1_dR * sim(sn).dR(i) + sim(sn).z1_dI * sim(sn).dI(i);
            sim(sn).bias(i,1) = sim(sn).x1_0 + sim(sn).x1_dR * sim(sn).dR(i) + sim(sn).x1_dI * sim(sn).dI(i);
            sim(sn).A(i,1)    = sim(sn).A1_0 + sim(sn).A1_dR * sim(sn).dR(i) + sim(sn).A1_dI * sim(sn).dI(i);
            sim(sn).T0(i,1)   = sim(sn).T01;
            
        else
            
            sim(sn).z(i,1)    = sim(sn).z6_0 + sim(sn).z6_dR * sim(sn).dR(i) + sim(sn).z6_dI * sim(sn).dI(i);
            sim(sn).bias(i,1) = sim(sn).x6_0 + sim(sn).x6_dR * sim(sn).dR(i) + sim(sn).x6_dI * sim(sn).dI(i);
            sim(sn).A(i,1)    = sim(sn).A6_0 + sim(sn).A6_dR * sim(sn).dR(i) + sim(sn).A6_dI * sim(sn).dI(i);
            sim(sn).T0(i,1)   = sim(sn).T06;
            
        end
        
        
        sim(sn).x0   = 2 ./ (1 + exp(-sim(sn).bias)) - 1;
        
    end
end

%% compute mean RT from DDM parameters
for sn = 1:length(sim)
    z = sim(sn).z;
    x0 = sim(sn).x0.*z;
    A = sim(sn).A;
    T0 = sim(sn).T0;
    
    a = A.^2;
    z = z./A;
    x0 = x0./A;
    
    RT = T0 + z .* tanh(z .* a) + ( 2*z.*(1-exp(-2*x0.*a))./( exp(2*z.*a) - exp(-2*z.*a) ) - x0  );
    ER = 1 ./ ( 1 + exp(2*z.*a) ) - ( (1-exp(-2*x0.*a))./(exp(2*z.*a) - exp(-2*z.*a) ) );
    sim(sn).RTana = RT;
    sim(sn).ER = ER;
end

%% compute simulated RTs
% for sn = 1:length(sim)
%
%     i1 = sim(sn).h == 1;
%     i6 = sim(sn).h == 6;
%     dR = sim(sn).dR(i1);
%     dI = sim(sn).dI(i1);
%
%     [ER, RT] = simulate_ERRT_v1(dR, dI, ...
%         sim(sn).A1_0, sim(sn).A1_dR, sim(sn).A1_dI, ...
%         sim(sn).z1_0, sim(sn).z1_dR, sim(sn).z1_dI, ...
%         sim(sn).x1_0, sim(sn).x1_dR, sim(sn).x1_dI, ...
%         sim(sn).T01);
%
%     sim(sn).ER(i1) = ER;
%     sim(sn).RTana(i1) = RT;
%
%     dR = sim(sn).dR(i6);
%     dI = sim(sn).dI(i6);
%
%     [ER, RT] = simulate_ERRT_v1(dR, dI, ...
%         sim(sn).A6_0, sim(sn).A6_dR, sim(sn).A6_dI, ...
%         sim(sn).z6_0, sim(sn).z6_dR, sim(sn).z6_dI, ...
%         sim(sn).x6_0, sim(sn).x6_dR, sim(sn).x6_dI, ...
%         sim(sn).T06);
%
%     sim(sn).ER(i6') = ER;
%     sim(sn).RTana(i6') = RT;
%
%     sim(sn).ER = sim(sn).ER';
%     sim(sn).RTana = sim(sn).RTana';
% end

%


