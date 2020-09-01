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

%% load the fake data
fname = 'DDM_fake.csv';
fid = fopen(fname);
hdr = textscan(fid, '%s%s%s%s%s%s', 1, 'delimiter', ',');
dat = textscan(fid, '%f%f%f%f%f%f', 'delimiter', ',');
fclose(fid);

sID = unique(dat{1});
sn = 1;
for i = 1:length(sID)
    ind = dat{1} == sID(i);
    fak(sn).sID = sID(i);
    fak(sn).dR = dat{2}(ind);
    fak(sn).dI = dat{3}(ind);
    fak(sn).gameLength = dat{4}(ind);
    fak(sn).choice = dat{5}(ind);
    fak(sn).RT = dat{6}(ind);
    sn = sn + 1;
end


%% load the "recovered" parameters
fname1 = 'parameterRecovery_2020-04-19_DDM_fake_fittedparams.csv';
fname2 = 'parameterRecovery_2020-04-19_DDM_fake_fittedT0.csv';
rec = load_SamsParameters_v1(fname1, fname2);

%% check that sIDs are same
if sum(abs([rec.sID] - [fak.sID]) ) ~= 0
    error('subject IDs not matched')
end

%% put dR and dI onto rec
for sn = 1:length(rec)
    rec(sn).dR = fak(sn).dR;
    rec(sn).dI = fak(sn).dI;
    rec(sn).gameLength = fak(sn).gameLength;
    rec(sn).h = fak(sn).gameLength-4;
end

%% how well do the recovered parameters fit the fake data?
rec = compute_meanERRT_v1(rec)

%% plot choice and RT curves comparing fak with rec
binEdges = [-25:10:25];
figure(1); clf;
set(gcf, 'position', [811   575   600   500])
ax = easy_gridOfEqualFigures([0.15 0.2 0.05], [0.12 0.15 0.05]);

% choices -----------------------------------------------------------------
e1 = plot_choiceCurvesFak_v1(ax(1:2), fak, binEdges, 0.1, 3);
e2 = plot_choiceCurves_sim_v3(ax(1:2), rec, binEdges, 0);

set(e1, 'linestyle', 'none', 'markersize', 30)
set(e2, 'marker', 'none')

set([e1([2 4]) e2([2 4])], 'visible', 'off')
e3 = plot_choiceCurvesFak_v1(ax(1:2), fak, binEdges, 0.1, 3);
e4 = plot_choiceCurves_sim_v3(ax(1:2), rec, binEdges, 0);
set([e3([1 3]) e4([1 3])], 'visible', 'off')
set(e3, 'linestyle', 'none', 'markersize', 30)
set(e4, 'marker', 'none')
leg = legend([e1([1]) e3(2) e2(1) e4(2)], ...
    {'fake h = 1' 'fake h = 6' 'model h = 1' 'model h = 6'}, ...
    'location', 'southeast');
set(leg, 'position', [0.3216    0.6614    0.1867    0.1310]);


% RTs ---------------------------------------------------------------------
e1 = plot_rtCurves_subFak_v1(ax(3:4), fak, 0.1 , 3)
e2 = plot_rtCurves_sim_v2(ax(3:4), rec, 1, 0)
set(e1, 'linestyle', 'none', 'markersize', 30)
set(e2, 'marker', 'none')

set([e1([2 4]) e2([2 4])], 'visible', 'off')

e1 = plot_rtCurves_subFak_v1(ax(3:4), fak, 0.1, 3)
e2 = plot_rtCurves_sim_v2(ax(3:4), rec, 1, 0)
set([e1([1 3]) e2([1 3])], 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)
set(e2, 'marker', 'none')

set(ax(3:4), 'ylim', [0.3 1.2])


axes(ax(1)); title('unequal [1 3], horizon 1')
axes(ax(2)); title('equal [2 2], horizon 1')
axes(ax(3)); title('unequal [1 3], horizon 6'); ylabel({'reaction time' '[seconds]'})
axes(ax(4)); title('equal [2 2], horizon 6'); ylabel({'reaction time' '[seconds]'})

saveFigurePdf(gcf, '~/Desktop/DDM_recVsFak')


%%
sn = 9;
ind = (fak(sn).dI == 0) & (fak(sn).gameLength == 5);

dR = fak(sn).dR(ind);
a = fak(sn).choice(ind);

[M, ~, X] = binIt(dR, a== 2, binEdges, 'std');

figure(1); clf;
plot(X, M)


