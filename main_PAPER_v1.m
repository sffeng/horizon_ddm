%% ========================================================================
%% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %%
%% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %%
%% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %% BASIC SETUP %%
%% ========================================================================
%%

% TODO - pair down data set

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
sub = load_humanData_v1(datadir, 'allHorizonData_v2.csv', 'DDM_demographics.csv', 0);

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



%% ========================================================================
%% illustration of model and theory %% illustration of model and theory %%
%% illustration of model and theory %% illustration of model and theory %%
%% illustration of model and theory %% illustration of model and theory %%
%% ========================================================================

%% FIGURE 1 - Schematic of the drift diffusion model showing the 
%% parameterization used in this paper 
x0 = -5;
A = 0.2;
z = 10;
rng(12)
T = 100;
r = randn(T,1) + A;

y = x0+[0; cumsum(r)];
ind = min(find(abs(y)>z));
y(ind) = z*sign(y(ind));
y(ind+1:end)= nan;


figure(1); clf;
set(gcf, 'Position', [811   572   600   250])
hg = [0.25 0.15];
wg = [0.37 0.1];
[ax, hb, wb] = easy_gridOfEqualFigures(hg, wg);
hold on;
plot([0:T], y, 'color', AZred, 'linewidth', 7)
plot([0 T], [1 1]*z,'k--', 'linewidth', 3)
plot([0 T], [1 1]*-z,'k--', 'linewidth', 3)
plot(ind, y(ind), '.', 'markersize', 50,'color', AZred)
text(T, z, ' left', 'fontsize', 20)
text(T, -z, ' right', 'fontsize', 20)
set(gca, 'ytick', [-z x0 z], ...
    'yticklabel', {'lower threshold, -\beta' 'starting point, x_0 = \alpha\beta' 'upper threshold, +\beta'}, ...
    'ylim', [-1.1 1.1]*z, 'tickdir', 'out', ...
    'xtick', [0:20:100], 'xticklabel', [0:20:100]/100, ...
    'fontsize', 20, 'linewidth', 3)
xlabel('time [seconds]')


x1 = 0;
x2 = ind/100;



X1 = x1*wb(1)+wg(1);
X2 = x2*wb(1)+wg(1)
YY = hg(1)+hb(1)+0.04;%hg(2)/2;
a = annotation('doublearrow',[X1 X2],[1 1]*YY);
text(ind/2, z*1.45, 'decision time, DT', ...
    'HorizontalAlignment', 'center', 'fontsize', 20)

x1 = 0;
x2 = 40/100;
X1 = x1*wb(1)+wg(1);
X2 = x2*wb(1)+wg(1)

y1 = (x0+z*1.1)/2.2/z;
y2 = (x0+40*A+z*1.1)/2.2/z;
Y1 = hg(1)+hb(1)*y1;
Y2 = hg(1)+hb(1)*y2;

b = annotation('arrow', [X1 X2], [Y1 Y2]);

text(x2*100, x0+40*A, ' drift rate, \mu', 'fontsize', 20)

% saveFigurePdf(gcf, '~/Desktop/DDMexample')

%% FIGURE 5 - Qualitative predictions for the logistic versions of the 
%% drift-diffusion model 
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

% just change baseline drift, c^\mu_0
rho = 0.01;
ft3.XXfit = theory_qualitativeRT_type2_v1(A1, A6, sigma1, sigma6, rho, 1)

% just change effect of dR on threshold, c^\beta_R
ft4.XXfit = theory_qualitativeRT_type2_v1(A1, A6, sigma1, sigma6, rho, 2)

% parameter key
% XX = [
%     cMu_0_1   cMu_0_6
%     cMu_R_1   cMu_R_6
%     cMu_I_1   cMu_I_6
%     cBeta_0_1 cBeta_0_6
%     cBeta_R_1 cBeta_R_6
%     cBeta_I_1 cBeta_I_6
%     0         0
%     0         0
%     0         0
%     0.05      0.05];

% just change drift, c^\mu_R - SNR
% 
%          0         0
%     0.0600    0.0356
%    -0.1000    0.2370
%     0.9000    0.9000
%          0         0
%          0         0
%          0         0
%          0         0
%          0         0
%     0.0500    0.0500
%     

% just change baseline threshold, c^\beta_0
%          0         0
%     0.0600    0.0600
%    -0.1000    0.4000
%     0.9000    0.5333
%          0         0
%          0         0
%          0         0
%          0         0
%          0         0
%     0.0500    0.0500

% just change baseline drift, c^\mu_0
%     2.3238    1.3771
%          0         0
%          0         0
%          0         0
%     0.0232    0.0232
%    -0.0387    0.1549
%          0         0
%          0         0
%          0         0
%     0.5000    0.5000

% just change effect of dR on threshold, c^\beta_R
%     2.3238    2.3238
%          0         0
%          0         0
%          0         0
%     0.0232    0.0138
%    -0.0387    0.0918
%          0         0
%          0         0
%          0         0
%     0.5000    0.5000




ft3.XXfit(end,:) = 0.5;
ft4.XXfit(end,:) = 0.5;

ft1 = theory_ERDT_v1(ft1, r_vals);
ft2 = theory_ERDT_v1(ft2, r_vals);
ft3 = theory_ERDT_v1(ft3, r_vals);
ft4 = theory_ERDT_v1(ft4, r_vals);


% start figure ------------------------------------------------------------
figure(1); clf;
set(gcf, 'position', [695   365   700   450])
hg = [0.1 0.13 0.35];
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

% delta2 = 0.33;
% wg3 = [wg(1) wg(2) wg(end)];
% wb3 = [wb(1) sum(wb(2:5))+sum(wg(3:5))];
% [b] = add_externalXlabels(wg3, hg, wb3, hb, delta2, ...
%     {'' 'response times' })
% set(b, 'fontsize', 20)
delta2 = 0.24;
wg2 = [wg(1) wg(2) wg(4) wg(end)];
wb2 = [wb(1) wb(2)+wb(3)+wg(3) wb(4)+wb(5)+wg(5)];
[b] = add_externalXlabels(wg2, hg, wb2, hb, delta2, ...
    {'' 
    'threshold independent of \DeltaR and \DeltaI' 
    'drift independent          of \DeltaR and \DeltaI'})

delta2 = 0.07;
wg2 = [wg(1) wg(2) wg(4) wg(end)];
wb2 = [wb(1) wb(2)+wb(3)+wg(3) wb(4)+wb(5)+wg(5)];
[b] = add_externalXlabels(wg2, hg, wb2, hb, delta2, ...
    {'' 
    'c^\beta_R = c^\beta_I = 0' 
    'c^\mu_R = c^\mu_I = 0'})


delta3 = -0.19;
[b] = add_externalXlabels(wg, hg, wb, hb, delta3, ...
    {'choice curves  all cases' 
    'drift change, c^\mu_R' 
    'threshold change, c^\beta_0'
    'drift change, c^\mu_0'
    'threshold change, c^\beta_R'})

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

addABCs(ax(1, :), [-0.02 0.165], 22)
saveFigurePdf(gcf, '~/Desktop/DDM_qualitative')




%% ========================================================================
%% logistic model %% logistic model %% logistic model %% logistic model %%
%% logistic model %% logistic model %% logistic model %% logistic model %%
%% logistic model %% logistic model %% logistic model %% logistic model %%
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
        cc13(:,i,sn) = fit_logistic(sn).cc13{i}(x_vals);
        cc22(:,i,sn) = fit_logistic(sn).cc22{i}(x_vals);
    end
end

CC13 = nanmean(cc13,3);
CC22 = nanmean(cc22,3);

%% FIGURE 3 - Choice behavior in the Horizon Task 
XX = cat(3, fit_logistic.x);

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
%% linear model of response times %% linear model of response times %%
%% linear model of response times %% linear model of response times %%
%% linear model of response times %% linear model of response times %%
%% ========================================================================

%% FIGURE 4 - RTs and RT curves
figure(1); clf;
set(gcf, 'position', [711   575   630   700])
ax = easy_gridOfEqualFigures([0.46 0.12 0.05], [0.15 0.15 0.08]);
ax(5:7) = easy_gridOfEqualFigures([0.08 0.72], [0.1 0.13 0.13 0.03]);
[l, leg] = plot_RToverTime_v1(ax(1:2), sub, 0.1, 3);
set(leg, 'position', [0.3168    0.8719    0.1595    0.0564]);

binEdges = [-25:10:25];
e = plot_RTvsdR_v1(ax(3:4), sub, binEdges, 0.1, 3);
plot_RTregression_v1(ax(5:7), sub, 0.1, 3)

axes(ax(1)); t = title('unequal information [1 3]')
axes(ax(2)); t(2) = title('equal information [2 2]')
axes(ax(5)); t(3) = title('baseline RT'); ylabel('\beta_0')
axes(ax(6)); t(4) = title('effect of a\DeltaR'); ylabel('\beta_R')
axes(ax(7)); t(5) = title('effect of a\DeltaI'); ylabel('\beta_I')

set(t, 'units', 'normalized');
get(t(3), 'position')
set(t(3:5), 'position', [ 0.5000 1.1 0])


set(ax, 'tickdir', 'out', 'fontsize', 16)
set(t, 'fontweight', 'normal', 'fontsize', 20)
addABCs(ax(1:4), [-0.09 0.05], 28)
addABCs(ax(5:7), [-0.075 0.06], 28, 'EFG')
saveFigurePdf(gcf, '~/Desktop/DDM_RTs')





%% ========================================================================
%% ABS DR DI MODEL %% ABS DR DI MODEL %% ABS DR DI MODEL %%
%% ABS DR DI MODEL %% ABS DR DI MODEL %% ABS DR DI MODEL %%
%% ABS DR DI MODEL %% ABS DR DI MODEL %% ABS DR DI MODEL %%
%% ========================================================================

%% fit the model with absolute value of dR on threshold
% takes about 30 seconds on my machine
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
fit_MLE_abs = fit_MLE_DDM_absDR_v1(sub, M, RTmin, RTmax, X0_all, LB_all, UB_all);
toc

%% FIGURE XXX - Drift-diffusion model parameters for abs(dR) abs(dI) model
clear X1 X6
vn1 = {
    'A1_0'
    'A1_dR'
    'A1_dI'
    'z1_0'
    'z1_dR'
    'z1_dI'
    'x1_0'
    'x1_dR'
    'x1_dI'
    'T01'
    };

vn6 = {
    'A6_0'
    'A6_dR'
    'A6_dI'
    'z6_0'
    'z6_dR'
    'z6_dI'
    'x6_0'
    'x6_dR'
    'x6_dI'
    'T06'
    };

% for i = 1:length(vn1)
%     for sn = 1:length(sim)
%         X1(i,sn) = getfield(sim(sn), vn1{i});
%         X6(i,sn) = getfield(sim(sn), vn6{i});
%     end
% end
for sn = 1:length(fit_MLE_abs)
    X1(:,sn) = fit_MLE_abs(sn).XXfit(:,1);
    X6(:,sn) = fit_MLE_abs(sn).XXfit(:,2);
end

figure(1); clf;
set(gcf, 'position', [811   176   550   650])
hg = ones(5,1)*0.1;
wg = ones(4,1)*0.07;
wg(1) = 0.24;
wg(end) = 0.02;
hg(1) = 0.07;
hg(end) = 0.09;
[~,hb,wb,ax] = easy_gridOfEqualFigures(hg, wg);
ax = ax';
ax = ax(:);
thresh = 0.01;
thresh2 = 0.005;
thresh3 = 0.001;
Ncomp = 10; % for multiple comparisons
set(ax(5), 'ylim', [-0.035 0.035]);
set(ax(8), 'ylim', [-0.06 0.09]);
for i = 1:size(X1,1)
    axes(ax(i)); hold on;
    R = (rand(1,size(X1,2))-0.5)*0.25;
    plot([1+R; 2+R], [X1(i,:); X6(i,:)], 'linewidth', 1, 'color', [1 1 1]*0.75)
    plot([1+R], [X1(i,:)], 'o', 'color', AZblue, 'linewidth', 1, 'markersize', 5)
    plot([2+R], [X6(i,:)], 'o', 'color', AZred, 'linewidth', 1, 'markersize', 5)
    plot([3+R], X6(i,:)-X1(i,:), 'o', 'color', AZsand, 'linewidth', 1, 'markersize', 5)
    if i == size(X1,1)
        xlabel('horizon')
    end
    plot([0.5 3.5], [0 0], 'k--', 'linewidth', 1)
end
% set(ax(8), 'ylim', [-0.06 0.08])
set(ax(9), 'ylim', [-0.6 0.8])
for i = 1:size(X1,1)
    axes(ax(i)); hold on;
    [~,p, ~, stats] = ttest([X1(i,:)-X6(i,:)]);
    
    if p < thresh/Ncomp
        str = '*';
        if p < thresh2/Ncomp
            str = '**';
            if p < thresh3/Ncomp
                str = '***';
            end
        end
        
        mx1 = 1.05 * max(X1(i,:));
        mx6 = 1.05 * max(X6(i,:));
        plot(1, mx1);
        yl = get(ax(i), 'ylim');
        mx = yl(1)+0.95*(yl(2)-yl(1));
        mx2 = yl(1)+0.97*(yl(2)-yl(1));
        mx1 = max(X1(i,:))+0.07*(yl(2)-yl(1));
        mx6 = max(X6(i,:))+0.07*(yl(2)-yl(1));
        
        plot([1 1 2 2], [mx1 mx mx mx6], 'k-', 'linewidth', 1);
        
        text(1.5, mx2, str, 'fontsize', 24, 'horizontalalignment', 'center')
        set(ax(i), 'color', 0.1*AZsand+0.9)
    end
    
end
for i = size(X1,1)+1:length(ax)
    set(ax(i), 'visible', 'off')
end
delta = 0.05;
a = annotation('textbox', [0 hg(1) wg(1)-delta hb(1)], 'string', 'non-decision time, T_0');
a(2) = annotation('textbox', [0 hg(1)+hg(2)+hb(1) wg(1)-delta hb(1)], 'string', 'bias, \alpha');
a(3) = annotation('textbox', [0 sum(hg(1:3))+sum(hb(1:2)) wg(1)-delta hb(1)], 'string', 'threshold, \beta');
a(4) = annotation('textbox', [0 sum(hg(1:4))+sum(hb(1:3)) wg(1)-delta hb(1)], 'string', 'drift rate, \mu');
set(a, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'right', ...
    'verticalalignment', 'middle')

delta2 = 0.06;
b = annotation('textbox', [wg(1) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'baseline');
b(2) = annotation('textbox', [wg(1)+wb(1)+wg(2) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of \DeltaR');
b(3) = annotation('textbox', [wg(1)+wb(1)+wg(2)+wb(2)+wg(3) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of \DeltaI');
set(b, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'center', ...
    'verticalalignment', 'middle')

set(ax, 'xtick', [1 2 3], 'xticklabel', {1 6 'diff'}, 'xlim', [0.5 3.5], ...
    'tickdir', 'out')
set(gcf, 'InvertHardcopy', 'off')


clear t
i=0;
i=i+1;axes(ax(i)); t(i) = title('c^\mu_0');
i=i+1;axes(ax(i)); t(i) = title('c^\mu_{R}');
i=i+1;axes(ax(i)); t(i) = title('c^\mu_{I}');
i=i+1;axes(ax(i)); t(i) = title('c^\beta_0');
i=i+1;axes(ax(i)); t(i) = title('c^\beta_{R}');
i=i+1;axes(ax(i)); t(i) = title('c^\beta_{I}');
i=i+1;axes(ax(i)); t(i) = title('c^\alpha_0');
i=i+1;axes(ax(i)); t(i) = title('c^\alpha_{R}');
i=i+1;axes(ax(i)); t(i) = title('c^\alpha_{I}');
i=i+1;axes(ax(i)); t(i) = title('T_0');
set(t, 'fontweight', 'normal', 'units', 'normalized')
set(t, 'position', [0.5000    1.08 0])
saveFigurePdf(gcf, '~/Desktop/DDM_params_absDR')

%% FIGURE XXX - posterior predictive
% compute theoretical average behavior for fit parameter values
r_vals = [-30:0.1:30];
fit_MLE_abs = theory_ERDT_abs_v1(fit_MLE_abs, r_vals);
% fit_MLEbest = theory_ERDT_v1(fit_MLEbest, r_vals);

ft = fit_MLE_abs;
dum = cat(3, ft.cc13_ana);
cc13 = nanmean(dum,3);
dum = cat(3, ft.cc22_ana);
cc22 = nanmean(dum,3);

dum = cat(3, ft.rt13_ana);
rt13 = nanmean(dum,3);
dum = cat(3, ft.rt22_ana);
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
set(e1, 'linestyle', 'none', 'markersize', 20)

e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1 , 3)
set(e1([2 4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 20)

e1 = plot_rtCurves_sub_v1(ax(5:6), sub, 0.1 , 3)
set(e1([1 3]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 20)

ff = 0.75;
set(l(1,:), 'color', AZblue*ff+(1-ff));
set(l(2,:), 'color', AZred*ff+(1-ff));

set(ax(3:end), 'ylim', [0.3 1.1])

axes(ax(1)); leg = legend([l(2,1) l(1,1) ], {'horizon 6' 'horizon 1'}, 'location', 'southeast');
axes(ax(1)); title(''); xlabel('R(high info) - R(low info)')
axes(ax(2)); title({''}); xlabel('R(left) - R(right)')
axes(ax(3)); title({'horizon 1'}); ylabel('RT [seconds]'); xlabel('R(high info) - R(low info)')
axes(ax(4)); title({''}); ylabel('RT [seconds]'); xlabel('R(left) - R(right)')
axes(ax(5)); title({'horizon 6'}); ylabel(''); xlabel('R(high info) - R(low info)')
axes(ax(6)); title({''}); ylabel(''); xlabel('R(left) - R(right)')

set(leg, 'position', [0.2956    0.5974    0.1375    0.0862])

set(ax, 'tickdir', 'out', 'xlim', [-35 35])
wg2 = [wg(1) wg(2) wg(4)];
wb2 = [wb(1) wb(2)+wb(3)+wg(3)];
delta2 = 0.06;
[b] = add_externalXlabels(wg2, hg, wb2, hb, delta2, {'choice' 'response time'})
set(b, 'linestyle', 'none')
delta = 0.075;
[a] = add_externalYlabels_v1(wg, hg, wb, hb, delta, {'equal condition [2 2]' 'unequal condition [1 3]'})

set([a b], 'fontsize', 18)
set(b, 'fontsize', 20)
set(ax, 'fontsize', 14)

saveFigurePdf(gcf, '~/Desktop/DDM_posteriorPredictiveABS')

%% parameter recovery for MLE fit -----------------------------------------
%% simulate fake data with MLE parameters
% TODO: 
%   1. could scramble parameters 
%   2. could use different set of trials

clear fak
dt = 0.001;
for sn = 1:length(fit_MLE_abs)
    fak_abs(sn) = simulate_Mmodel_abs_v1(fit_MLE_abs(sn).XXfit, dt, fit_MLE_abs(sn).dR, fit_MLE_abs(sn).dI, fit_MLE_abs(sn).gameLength);
end

%% fit fake data with MLE
fit_fak_abs = fit_MLE_DDM_absDR_v1(fak_abs, M, RTmin, RTmax, X0_all, LB_all, UB_all);

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
[l, r, p] = plot_compareTwoParameterSets(ax, fit_MLE_abs, fit_fak_abs, names, 'sim', 'fit');

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

saveFigurePdf(gcf, '~/Desktop/DDM_parameterRecoveryABS_MLE')

%% model comparison between modified and original model
mean(vertcat(fit_MLE.BIC) < vertcat(fit_MLE_abs.BIC))

%% noise check ------------------------------------------------------------
%% First simulate behavior leaving either c_R^\thresh or c_R^\drift at horizon 1 values


clear fak
dt = 0.001;
for sn = 1:length(fit_MLE_abs)
    
    % parameters
    XX = fit_MLE_abs(sn).XXfit;
    
    % keep c_R^\thresh constant with horizon cZ_R = XX(5,:);
    XX_constantThreshold = XX;
    XX_constantThreshold(5,2) = XX_constantThreshold(5,1);
    
    % keep c_R^\thresh constant with horizon cA_R = XX(2,:);
    XX_constantDrift = XX;
    XX_constantDrift(2,2) = XX_constantDrift(2,1);
    
    
    fak_abs_constantThreshold(sn) = simulate_Mmodel_abs_v1(XX_constantThreshold, dt, fit_MLE_abs(sn).dR, fit_MLE_abs(sn).dI, fit_MLE_abs(sn).gameLength);
    fak_abs_constantDrift(sn)     = simulate_Mmodel_abs_v1(XX_constantDrift, dt, fit_MLE_abs(sn).dR, fit_MLE_abs(sn).dI, fit_MLE_abs(sn).gameLength);
end

%% Then fit behavior with logisitic model to estimate noise
priorFlag = 1;
for sn = 1:length(sub)
    fit_logistic_constantThreshold(sn) = fit_logistic_XXform_v1(fak_abs_constantThreshold(sn), priorFlag);
    fit_logistic_constantDrift(sn)     = fit_logistic_XXform_v1(fak_abs_constantDrift(sn), priorFlag);
end

%% Then compare fit noise to noise fit to simulations with full model and fit to human data
XX = cat(3, fit_logistic.x);
bonus   = squeeze(XX(3,:,:));
noise13 = squeeze(XX(2,:,:));
noise22 = squeeze(XX(5,:,:));

XX_constantThreshold = cat(3, fit_logistic_constantThreshold.x);
bonus_constantThreshold   = squeeze(XX_constantThreshold(3,:,:));
noise13_constantThreshold = squeeze(XX_constantThreshold(2,:,:));
noise22_constantThreshold = squeeze(XX_constantThreshold(5,:,:));

XX_constantDrift = cat(3, fit_logistic_constantDrift.x);
bonus_constantDrift   = squeeze(XX_constantDrift(3,:,:));
noise13_constantDrift = squeeze(XX_constantDrift(2,:,:));
noise22_constantDrift = squeeze(XX_constantDrift(5,:,:));
% fit_fak_abs_constantThreshold = fit_MLE_DDM_absDR_v1(fak_abs_constantThreshold, M, RTmin, RTmax, X0_all, LB_all, UB_all);
% fit_fak_abs_constantDrift = fit_MLE_DDM_absDR_v1(fak_abs_constantDrift, M, RTmin, RTmax, X0_all, LB_all, UB_all);



M13 = [nanmean(noise13,2) nanmean(noise13_constantThreshold,2) nanmean(noise13_constantDrift,2)];
S13 = [nanstd(noise13,[],2) nanstd(noise13_constantThreshold,[],2) nanstd(noise13_constantDrift,[],2)]/sqrt(length(sub));

M22 = [nanmean(noise22,2) nanmean(noise22_constantThreshold,2) nanmean(noise22_constantDrift,2)];
S22 = [nanstd(noise22,[],2) nanstd(noise22_constantThreshold,[],2) nanstd(noise22_constantDrift,[],2)]/sqrt(length(sub));

figure(1); clf; hold on
errorbar(M22, S22)

%% compare with sigma ratio horizon 1 / horizon 6
R13 = [noise13(1,:)./noise13(2,:); noise13_constantThreshold(1,:)./noise13_constantThreshold(2,:); noise13_constantDrift(1,:)./noise13_constantDrift(2,:)];
R22 = [noise22(1,:)./noise22(2,:); noise22_constantThreshold(1,:)./noise22_constantThreshold(2,:); noise22_constantDrift(1,:)./noise22_constantDrift(2,:)];

figure(1); clf;
set(gcf, 'position', [1     1   588   307])
hb = [0.1 0.1];
wb = [0.1 0.1 0.03];
[ax, wg, hg] = easy_gridOfEqualFigures(hb, wb);
axes(ax(1)); hold on;
rang = 0.125;
l1 = plot(repmat([1 2 3]', [1 size(R13,2)])'+(rand(size(R13'))-0.5)*rang, R13','o');
% l2 = plot(repmat([4 5 6]', [1 size(R13,2)])'+(rand(size(R13'))-0.5)*rang, R13','o');
plot([0 4], [1 1], 'k--', 'linewidth', 1)

set([l1], 'markersize', 5, 'linewidth', 1)
b = boxplot([R13' ], 'notch', 'on', 'widths', 0.25);
set(b(7,:), 'visible', 'off')
set(b, 'color', 'k', 'linewidth', 2)
ylabel({'\sigma ratio' 'horizon 1 / horizon 6'})

t = title('[1 3] condition');


axes(ax(2)); hold on;
rang = 0.125;
% l1 = plot(repmat([1 2 3]', [1 size(R13,2)])'+(rand(size(R13'))-0.5)*rang, R13','o');
l2 = plot(repmat([1 2 3]', [1 size(R13,2)])'+(rand(size(R13'))-0.5)*rang, R13','o');

set([l2], 'markersize', 5, 'linewidth', 1)
b = boxplot([R22' ], 'notch', 'on', 'widths', 0.25);
set(b(7,:), 'visible', 'off')
set(b, 'color', 'k', 'linewidth', 2)
ylabel({'\sigma ratio' 'horizon 1 / horizon 6'})
plot([0 4], [1 1], 'k--', 'linewidth', 1)

t(2) = title('[2 2] condition');

set(ax, 'tickdir', 'out', 'xtick', [1:3], ...
    'xticklabel', {'humans' 'c_R^\beta constant' 'c_R^\mu constant'}, ...
    'view', [90 90], 'yaxislocation', 'left')
% set(ax(2), 'yticklabel')
set(ax, 'ylim', [-0.2 2.2], 'box', 'off')
ax(1).TickLabelInterpreter='tex';
ax(2).TickLabelInterpreter='tex';
set(t, 'fontweight', 'normal')

set([l1(1) l2(1)], 'color', AZred)
set([l1(2) l2(2)], 'color', AZblue)
set([l1(3) l2(3)], 'color', AZsand)

addABCs(ax, [-0.14 0.08], 28)
saveFigurePdf(gcf, '~/Desktop/ABS_ratios')

% add_externalYlabels_v1(hg, wg, hb, wb, 0, {'a' 'b'})

%% ========================================================================
%% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %%
%% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %%
%% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %% MLE FIT %%
%% ========================================================================



%% fit the full model
% takes about 30 seconds on my machine
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

%% FIGURE 6 - Drift-diffusion model parameters in the Horizon Task
clear X1 X6
vn1 = {
    'A1_0'
    'A1_dR'
    'A1_dI'
    'z1_0'
    'z1_dR'
    'z1_dI'
    'x1_0'
    'x1_dR'
    'x1_dI'
    'T01'
    };

vn6 = {
    'A6_0'
    'A6_dR'
    'A6_dI'
    'z6_0'
    'z6_dR'
    'z6_dI'
    'x6_0'
    'x6_dR'
    'x6_dI'
    'T06'
    };

% for i = 1:length(vn1)
%     for sn = 1:length(sim)
%         X1(i,sn) = getfield(sim(sn), vn1{i});
%         X6(i,sn) = getfield(sim(sn), vn6{i});
%     end
% end
for sn = 1:length(fit_MLE)
    X1(:,sn) = fit_MLE(sn).XXfit(:,1);
    X6(:,sn) = fit_MLE(sn).XXfit(:,2);
end

figure(1); clf;
set(gcf, 'position', [811   176   550   650])
hg = ones(5,1)*0.1;
wg = ones(4,1)*0.07;
wg(1) = 0.24;
wg(end) = 0.02;
hg(1) = 0.07;
hg(end) = 0.09;
[~,hb,wb,ax] = easy_gridOfEqualFigures(hg, wg);
ax = ax';
ax = ax(:);
thresh = 0.01;
thresh2 = 0.005;
thresh3 = 0.001;
Ncomp = 10; % for multiple comparisons
for i = 1:size(X1,1)
    axes(ax(i)); hold on;
    R = (rand(1,size(X1,2))-0.5)*0.25;
    plot([1+R; 2+R], [X1(i,:); X6(i,:)], 'linewidth', 1, 'color', [1 1 1]*0.75)
    plot([1+R], [X1(i,:)], 'o', 'color', AZblue, 'linewidth', 1, 'markersize', 5)
    plot([2+R], [X6(i,:)], 'o', 'color', AZred, 'linewidth', 1, 'markersize', 5)
    plot([3+R], X6(i,:)-X1(i,:), 'o', 'color', AZsand, 'linewidth', 1, 'markersize', 5)
    if i == size(X1,1)
        xlabel('horizon')
    end
    plot([0.5 3.5], [0 0], 'k--', 'linewidth', 1)
end
set(ax(8), 'ylim', [-0.06 0.08])
set(ax(9), 'ylim', [-0.6 0.8])
for i = 1:size(X1,1)
    axes(ax(i)); hold on;
    [~,p, ~, stats] = ttest([X1(i,:)-X6(i,:)]);
    
    if p < thresh/Ncomp
        str = '*';
        if p < thresh2/Ncomp
            str = '**';
            if p < thresh3/Ncomp
                str = '***';
            end
        end
        
        mx1 = 1.05 * max(X1(i,:));
        mx6 = 1.05 * max(X6(i,:));
        plot(1, mx1);
        yl = get(ax(i), 'ylim');
        mx = yl(1)+0.95*(yl(2)-yl(1));
        mx2 = yl(1)+0.97*(yl(2)-yl(1));
        mx1 = max(X1(i,:))+0.07*(yl(2)-yl(1));
        mx6 = max(X6(i,:))+0.07*(yl(2)-yl(1));
        
        plot([1 1 2 2], [mx1 mx mx mx6], 'k-', 'linewidth', 1);
        
        text(1.5, mx2, str, 'fontsize', 24, 'horizontalalignment', 'center')
        set(ax(i), 'color', 0.1*AZsand+0.9)
    end
    
end
for i = size(X1,1)+1:length(ax)
    set(ax(i), 'visible', 'off')
end
delta = 0.05;
a = annotation('textbox', [0 hg(1) wg(1)-delta hb(1)], 'string', 'non-decision time, T_0');
a(2) = annotation('textbox', [0 hg(1)+hg(2)+hb(1) wg(1)-delta hb(1)], 'string', 'bias, \alpha');
a(3) = annotation('textbox', [0 sum(hg(1:3))+sum(hb(1:2)) wg(1)-delta hb(1)], 'string', 'threshold, \beta');
a(4) = annotation('textbox', [0 sum(hg(1:4))+sum(hb(1:3)) wg(1)-delta hb(1)], 'string', 'drift rate, \mu');
set(a, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'right', ...
    'verticalalignment', 'middle')

delta2 = 0.06;
b = annotation('textbox', [wg(1) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'baseline');
b(2) = annotation('textbox', [wg(1)+wb(1)+wg(2) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of \DeltaR');
b(3) = annotation('textbox', [wg(1)+wb(1)+wg(2)+wb(2)+wg(3) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of \DeltaI');
set(b, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'center', ...
    'verticalalignment', 'middle')

set(ax, 'xtick', [1 2 3], 'xticklabel', {1 6 'diff'}, 'xlim', [0.5 3.5], ...
    'tickdir', 'out')
set(gcf, 'InvertHardcopy', 'off')


clear t
i=0;
i=i+1;axes(ax(i)); t(i) = title('c^\mu_0');
i=i+1;axes(ax(i)); t(i) = title('c^\mu_{R}');
i=i+1;axes(ax(i)); t(i) = title('c^\mu_{I}');
i=i+1;axes(ax(i)); t(i) = title('c^\beta_0');
i=i+1;axes(ax(i)); t(i) = title('c^\beta_{R}');
i=i+1;axes(ax(i)); t(i) = title('c^\beta_{I}');
i=i+1;axes(ax(i)); t(i) = title('c^\alpha_0');
i=i+1;axes(ax(i)); t(i) = title('c^\alpha_{R}');
i=i+1;axes(ax(i)); t(i) = title('c^\alpha_{I}');
i=i+1;axes(ax(i)); t(i) = title('T_0');
set(t, 'fontweight', 'normal', 'units', 'normalized')
set(t, 'position', [0.5000    1.08 0])
saveFigurePdf(gcf, '~/Desktop/DDM_params2')

%% STATS
% both $\rb^\drift_R$ and $\rb^\thresh_0$ decrease with horizon (STATS)
% $\rb^\drift_R$
i = 2;
[~,p, ~, stats] = ttest([X1(i,:)-X6(i,:)]);
disp(sprintf('$t(%d) = %.2f$, $p < 0.001$', stats.df, stats.tstat))
i = 3;
[~,p, ~, stats] = ttest([X1(i,:)-X6(i,:)]);
disp(sprintf('$t(%d) = %.2f$, $p < 0.001$', stats.df, -stats.tstat))
% $\rb^\thresh_0$ 
i = 4;
[~,p, ~, stats] = ttest([X1(i,:)-X6(i,:)]);
disp(sprintf('$t(%d) = %.2f$, $p < 0.001$', stats.df, stats.tstat))

%% posterior predictive checks --------------------------------------------
%% compute theoretical average behavior for fit parameter values
r_vals = [-30:0.1:30];
fit_MLE = theory_ERDT_v1(fit_MLE, r_vals);
% fit_MLEbest = theory_ERDT_v1(fit_MLEbest, r_vals);

%% FIGURE 7 - The drift-diffusion model (solid lines) captures the main 
%% qualitative effects of human choice and response time data (dots) in all 
%% conditions of the experiment.
ft = fit_MLE;
dum = cat(3, ft.cc13_ana);
cc13 = nanmean(dum,3);
dum = cat(3, ft.cc22_ana);
cc22 = nanmean(dum,3);

dum = cat(3, ft.rt13_ana);
rt13 = nanmean(dum,3);
dum = cat(3, ft.rt22_ana);
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
set(e1, 'linestyle', 'none', 'markersize', 20)

e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1 , 3)
set(e1([2 4]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 20)

e1 = plot_rtCurves_sub_v1(ax(5:6), sub, 0.1 , 3)
set(e1([1 3]), 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 20)

ff = 0.75;
set(l(1,:), 'color', AZblue*ff+(1-ff));
set(l(2,:), 'color', AZred*ff+(1-ff));

set(ax(3:end), 'ylim', [0.3 1.1])

axes(ax(1)); leg = legend([l(2,1) l(1,1) ], {'horizon 6' 'horizon 1'}, 'location', 'southeast');
axes(ax(1)); title(''); xlabel('R(high info) - R(low info)')
axes(ax(2)); title({''}); xlabel('R(left) - R(right)')
axes(ax(3)); title({'horizon 1'}); ylabel('RT [seconds]'); xlabel('R(high info) - R(low info)')
axes(ax(4)); title({''}); ylabel('RT [seconds]'); xlabel('R(left) - R(right)')
axes(ax(5)); title({'horizon 6'}); ylabel(''); xlabel('R(high info) - R(low info)')
axes(ax(6)); title({''}); ylabel(''); xlabel('R(left) - R(right)')

set(leg, 'position', [0.2956    0.5974    0.1375    0.0862])

set(ax, 'tickdir', 'out', 'xlim', [-35 35])
wg2 = [wg(1) wg(2) wg(4)];
wb2 = [wb(1) wb(2)+wb(3)+wg(3)];
delta2 = 0.06;
[b] = add_externalXlabels(wg2, hg, wb2, hb, delta2, {'choice' 'response time'})
set(b, 'linestyle', 'none')
delta = 0.075;
[a] = add_externalYlabels_v1(wg, hg, wb, hb, delta, {'equal condition [2 2]' 'unequal condition [1 3]'})

set([a b], 'fontsize', 18)
set(b, 'fontsize', 20)
set(ax, 'fontsize', 14)

saveFigurePdf(gcf, '~/Desktop/DDM_posteriorPredictiveMLEbest')
%% -----------------------------------------end posterior predictive checks

%% Sensitivity analysis ---------------------------------------------------
%% FIGURE 8 - Sensitivity analysis
XX = cat(3,fit_MLE.XXfit);

iZ_0 = 4;
iA_dR = 2;

x1 = squeeze(XX(iA_dR,:,:));
x2 = squeeze(XX(iZ_0,:,:));

rA = x1(1,:) ./ x1(2,:);
rZ = x2(1,:) ./ x2(2,:);

X = cat(3,fit_logistic.x);
i13 = 2;
i22 = 5;

n13 = squeeze(X(i13, :, :));
n22 = squeeze(X(i22, :, :));

r13 = n13(2,:) ./ n13(1,:);
r22 = n22(2,:) ./ n22(1,:);


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

%% FIGURE 8 - Sensitivity analysis
XX = cat(3,fit_MLE.XXfit);

iZ_0 = 4;
iA_dR = 2;

x1 = squeeze(XX(iA_dR,:,:));
x2 = squeeze(XX(iZ_0,:,:));

rA = x1(2,:) ./ x1(1,:);
rZ = x2(2,:) ./ x2(1,:);

X = cat(3,fit_logistic.x);
i13 = 2;
i22 = 5;

n13 = squeeze(X(i13, :, :));
n22 = squeeze(X(i22, :, :));

r13 = n13(1,:) ./ n13(2,:);
r22 = n22(1,:) ./ n22(2,:);


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
ylabel({'\sigma ratio' 'horizon 1 / horizon 6'})
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
ylabel({'parameter ratio' 'horizon 6 / horizon 1 '})
plot([1 1 2 2], [1.45 1.6 1.6 1.45], 'k-', 'linewidth', 1)
text(1.5, 1.6, '***', 'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
    'fontsize', 28)

set(ax(1), 'xtick', [1 2], 'xticklabel', {'[1 3]' '[2 2]'})
set(ax(2), 'xtick', [1 2], 'xticklabel', {'c^\mu_R' 'c^\beta_0'})
set(ax, 'ylim', [0 1.8], 'box', 'off', 'tickdir', 'out')
ax(2).TickLabelInterpreter='tex';
set(ax, 'ytick', [0:.5:9])

addABCs(ax, [-0.08 0.1], 28)
saveFigurePdf(gcf, '~/Desktop/DDM_sigmaVsDDMpars')
% set(gca, 'xscale', 'log', 'yscale', 'log')

%% simulate fake data with only c^\mu_R changing
for sn = 1:length(fit_MLE)
    sens1(sn) = fit_MLE(sn);
    xx = sens1(sn).XXfit;
    sens1(sn).XXfit = [xx(:,1) xx(:,1)];
    sens1(sn).XXfit(2,2) = xx(2,2);
    sens1(sn).XXfit(3,2) = xx(3,2);
end
%%
clear fak
dt = 0.001;
for sn = 1:length(sens1)
    sens1_fak(sn) = simulate_Mmodel_v1(sens1(sn).XXfit, dt, sens1(sn).dR, sens1(sn).dI, sens1(sn).gameLength);
end
%% fit fake data with MLE
priorFlag = 1;
for sn = 1:length(sub)
    fit_sens1_fak(sn) = fit_logistic_XXform_v1(sens1_fak(sn), priorFlag)
end
%%
sensXX = cat(3,fit_sens1_fak.x);
s1_13 = squeeze(sensXX(2,:,:));
s1_22 = squeeze(sensXX(5,:,:));
logXX = cat(3,fit_logistic.x);
l_13 = squeeze(logXX(2,:,:));
l_22 = squeeze(logXX(5,:,:));
x13 = (l_13(2,:)-l_13(1,:))./(l_13(2,:)+l_13(1,:));
s13 = (s1_13(2,:)-s1_13(1,:))./(s1_13(2,:)+s1_13(1,:));

%% ----------------------------------------------- end sensitivity analysis

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


