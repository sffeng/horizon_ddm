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

%% load
% load and augment data
sub = load_E1_v2([datadir 'allHorizonData_v2.csv']);


% remove bad subjects
% [sub, sub_bad] = removeBadSubjects_E1_v1(sub);

% remove bad subjects
% [sub, sub_bad] = removeBadSubjects_E1_v2(sub,0.6);

% remove subjects with same random seed
% sub = removeSameRseedSubjects_E1_v1(sub);

% princeton subjects only
i1 = strcmp({sub.expt_name}, 'pilot-v1');% 'repeater-v1'})
i2 = strcmp({sub.expt_name}, 'repeater-v1');
% i3 = strcmp({sub.expt_name}, 'short-pilot')
sub = sub(i1 | i2 );
% 'repeater-v1'    'short-pilot'

%% load demographics
fname = 'DDM_demographics.csv';
fid = fopen(fname);
hdr = textscan(fid, '%s%s%s', 1, 'delimiter', ',');
data = textscan(fid, '%s%f%s', 'delimiter', ',');
fclose(fid)

for i = 1:length(data{1})
    nm{i} = data{1}{i};
    ge{i} = data{3}{i};
end
ag = data{2};

%% pair demographics to subject
for sn = 1:length(sub)
    idx = find(strcmp(nm, sub(sn).subjectID));
    sub(sn).age = ag(idx);
    sub(sn).gender = ge{idx};
end


%% remove fast or slow RTs
for sn = 1:length(sub)
    ind = (sub(sn).RT(:,5) > 0.1) & (sub(sn).RT(:,5) < 3);
    sub(sn).a(~ind,5) = nan;
    sub(sn).nExcluded = length(ind) - sum(~ind);
end

%% load Sam's new parameters
% fname1 = 'fittedparams_nt131filtered.csv';
% fname2 = 'fittedT0_nt131filtered.csv';
fname1 = 'fittedparams.csv';
fname2 = 'fittedT0.csv';

% DDM EXP val subj_idx gLen
fid = fopen(fname1);
hdr_dum = textscan(fid,...
    '%s%s%s%s%s', ...
    1, 'delimiter', ',');
data = textscan(fid,...
    '%s%s%f%f%f', ...
    'delimiter', ',');
fclose(fid);

for i = 1:length(hdr_dum)
    hdr{i} = hdr_dum{i}{1};
end
iDDM = find(strcmp(hdr, 'DDM'));
iEXP = find(strcmp(hdr, 'EXP'));
iVAL = find(strcmp(hdr, 'val'));
iSub = find(strcmp(hdr, 'subj_idx'));
iGL = find(strcmp(hdr, 'gLen'));
%
% drift	0
% drift	oDel
% drift	cInfo
% thresh	0
% thresh	oDel
% thresh	cInfo
% bias	0
% bias	oDel
% bias	cInfo
% drift	0
% drift	oDel
% drift	cInfo
% thresh	0
% thresh	oDel
% thresh	cInfo
% bias	0
% bias	oDel
% bias	cInfo
U = unique(data{iSub});
for sn = 1:length(U)
    
    idx = find(data{iSub} == U(sn));
    dum = data{iVAL}(idx);
    sim(sn).sID = U(sn);
    
    sim(sn).A1_0  = dum(10);
    sim(sn).A1_dR = dum(11);
    sim(sn).A1_dI = -dum(12);
    sim(sn).A6_0  = dum(1);
    sim(sn).A6_dR = dum(2);
    sim(sn).A6_dI = -dum(3);
    
    sim(sn).z1_0  = dum(13);
    sim(sn).z1_dR = dum(14);
    sim(sn).z1_dI = -dum(15);
    sim(sn).z6_0  = dum(4);
    sim(sn).z6_dR = dum(5);
    sim(sn).z6_dI = -dum(6);
    
    sim(sn).x1_0  = dum(16);
    sim(sn).x1_dR = dum(17);
    sim(sn).x1_dI = -dum(18);
    sim(sn).x6_0  = dum(7);
    sim(sn).x6_dR = dum(8);
    sim(sn).x6_dI = -dum(9);
end

% now load T0
% param	val	subj_idx	gLen
fid = fopen(fname2);
hdr_dum = textscan(fid,...
    '%s%s%s%s', ...
    1, 'delimiter', ',');
data = textscan(fid,...
    '%s%f%f%f', ...
    'delimiter', ',');
fclose(fid);

for i = 1:length(hdr_dum)
    hdr{i} = hdr_dum{i}{1};
end
iVAL = find(strcmp(hdr, 'val'));
iSub = find(strcmp(hdr, 'subj_idx'));
iGL = find(strcmp(hdr, 'gLen'));

U = unique(data{iSub});
for sn = 1:length(U)
    
    idx = find(data{iSub} == U(sn));
    dum = data{iVAL}(idx);
    sim(sn).T01 = dum(2);
    sim(sn).T06 = dum(1);
    
end

%% match sim to sub
sID = [sim.sID];
[a,idx,c] = intersect([sub.sNum], sID);
sub_all = sub;
for sn = 1:length(sub);%length(idx)
    if sum(idx==sn) > 0
        sub_all(sn).excluded = 0;
    else
        sub_all(sn).excluded = 1;
    end
end

sub = sub(idx);

for sn = 1:length(sim)
    sim(sn).dR = sub(sn).o2(:,4) - sub(sn).o1(:,4);
    sim(sn).dI = (sub(sn).n2(:,4) - sub(sn).n1(:,4)) / 2;
    %sim(sn).RT = sub(sn).RT(:,5);
    sim(sn).gameLength = sub(sn).gameLength;
    sim(sn).h = sub(sn).gameLength-4;
end

%% report demographics
% all
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
saveFigurePdf(gcf, '~/Desktop/DDM_demographics')

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





%% ========================================================================
%% simulate for parameter recovery %% simulate for parameter recovery %%
%% simulate for parameter recovery %% simulate for parameter recovery %%
%% simulate for parameter recovery %% simulate for parameter recovery %%
%% ========================================================================

%% check to see whether simulations generate correct mean RTs and ERs

clear RT C
dt = 0.001;

% fak(sn)
sn = 3;
dR = sim(sn).dR;
dI = sim(sn).dI;
i1 = sim(sn).gameLength == 5;
i6 = sim(sn).gameLength == 10;

for i = 1:length(dR)
    i
    for j = 1:100
        [~, ~, DT(i,j), C(i,j)] = simluate_DDM_v1(dt, sim(sn).A(i), 1, sim(sn).z(i), sim(sn).x0(i));
        RT(i,j) = DT(i,j) + sim(sn).T0(i);
    end
end

rt = nanmean(RT,2)
er = 1-nanmean(C,2);

figure(1); clf;
ax = easy_gridOfEqualFigures([0.1 0.1], [0.1 0.1 0.03]);
axes(ax(1)); hold on;
plot(sim(sn).ER(1:length(er)), er,'.')
xl = get(gca, 'xlim');
plot(xl, xl, 'k--')

axes(ax(2)); hold on;
plot(sim(sn).RTana(1:length(rt)), rt,'.')
xl = get(gca, 'xlim');
plot(xl, xl, 'k--')


figure(2); clf;
hold on;
plot(er, rt,'.')
plot(sim(sn).ER, sim(sn).RTana,'.')

%% create a fake data set based on what subjects see
clear fak
dt = 0.00001;

% fak(sn)
for sn = 1:length(sim)
    
    clear RT C
    dR = sim(sn).dR;
    dI = sim(sn).dI;
    i1 = sim(sn).gameLength == 5;
    i6 = sim(sn).gameLength == 10;
    
    for i = 1:length(dR)
        [~, ~, DT(i), C(i)] = simluate_DDM_v1(dt, sim(sn).A(i), 1, sim(sn).z(i), sim(sn).x0(i));
        RT(i) = DT(i) + sim(sn).T0(i);
    end
    
    fak(sn).sID(1:length(dR)) = sim(sn).sID;
    fak(sn).dR = dR;
    fak(sn).dI = dI;
    fak(sn).gameLength = sim(sn).gameLength;
    fak(sn).choice = C + 1;
    fak(sn).RT = RT;
    
end

%% write fake data out to CSV
fname = '~/Desktop/DDM_fake.csv';

n = fieldnames(fak(1));

hdr_str = '%s,%s,%s,%s,%s,%s\n';
num_str = '%d,%f,%d,%d,%d,%f\n';

fid = fopen(fname, 'w');
fprintf(fid, hdr_str, n{:});

for sn = 1:length(fak)
    clear X
    for i = 1:length(n)
        X(i,:) = getfield(fak(sn), n{i});
    end
    
    fprintf(fid, num_str, X);
end
fclose(fid);

%% load parameter recovery fits
% fname1 = 'fittedparams_nt131filtered.csv';
% fname2 = 'fittedT0_nt131filtered.csv';
% fname1 = 'parameterRecovery_2020-04-19_DDM_fake_fittedparams.csv';
% fname2 = 'parameterRecovery_2020-04-19_DDM_fake_fittedT0.csv';
% fname1 = 'parameterRecovery2_fittedparams.csv';
% fname2 = 'parameterRecovery2_fittedT0.csv';
fname1 = 'AAA_fittedparams.csv';
fname2 = 'AAA_fittedT0.csv';

% DDM EXP val subj_idx gLen
fid = fopen(fname1);
hdr_dum = textscan(fid,...
    '%s%s%s%s%s', ...
    1, 'delimiter', ',');
data = textscan(fid,...
    '%s%s%f%f%f', ...
    'delimiter', ',');
fclose(fid);

for i = 1:length(hdr_dum)
    hdr{i} = hdr_dum{i}{1};
end
iDDM = find(strcmp(hdr, 'DDM'));
iEXP = find(strcmp(hdr, 'EXP'));
iVAL = find(strcmp(hdr, 'val'));
iSub = find(strcmp(hdr, 'subj_idx'));
iGL = find(strcmp(hdr, 'gLen'));
%
% drift	0
% drift	oDel
% drift	cInfo
% thresh	0
% thresh	oDel
% thresh	cInfo
% bias	0
% bias	oDel
% bias	cInfo
% drift	0
% drift	oDel
% drift	cInfo
% thresh	0
% thresh	oDel
% thresh	cInfo
% bias	0
% bias	oDel
% bias	cInfo

% drift	0
% drift	oDel
% drift	cInfo
% thresh	0
% thresh	oDel
% thresh	cInfo
% bias	0
% bias	oDel
% bias	cInfo

U = unique(data{iSub});
for sn = 1:length(U)
    
    idx = find(data{iSub} == U(sn));
    dum = data{iVAL}(idx);
    par(sn).sID = U(sn);
    
    par(sn).A1_0  = dum(10);
    par(sn).A1_dR = dum(11);
    par(sn).A1_dI = -dum(12);
    par(sn).A6_0  = dum(1);
    par(sn).A6_dR = dum(2);
    par(sn).A6_dI = -dum(3);
    
    par(sn).z1_0  = dum(13);
    par(sn).z1_dR = dum(14);
    par(sn).z1_dI = -dum(15);
    par(sn).z6_0  = dum(4);
    par(sn).z6_dR = dum(5);
    par(sn).z6_dI = -dum(6);
    
    par(sn).x1_0  = dum(16);
    par(sn).x1_dR = dum(17);
    par(sn).x1_dI = -dum(18);
    par(sn).x6_0  = dum(7);
    par(sn).x6_dR = dum(8);
    par(sn).x6_dI = -dum(9);
end

% now load T0
% param	val	subj_idx	gLen
fid = fopen(fname2);
hdr_dum = textscan(fid,...
    '%s%s%s%s', ...
    1, 'delimiter', ',');
data = textscan(fid,...
    '%s%f%f%f', ...
    'delimiter', ',');
fclose(fid);

for i = 1:length(hdr_dum)
    hdr{i} = hdr_dum{i}{1};
end
iVAL = find(strcmp(hdr, 'val'));
iSub = find(strcmp(hdr, 'subj_idx'));
iGL = find(strcmp(hdr, 'gLen'));

U = unique(data{iSub});
for sn = 1:length(U)
    
    idx = find(data{iSub} == U(sn));
    dum = data{iVAL}(idx);
    par(sn).T01 = dum(2);
    par(sn).T06 = dum(1);
    
end

%%
figure(1); clf;

nm = fieldnames(par);

hg = ones(4, 1)*0.1;
wg = ones(8, 1)*0.05;
ax = easy_gridOfEqualFigures(hg, wg);

clear X Y
for i = 1:length(nm)
    for j = 1:length(par)
        Y(j) = getfield(par(j), nm{i});
        X(j) = getfield(sim(j), nm{i});
        
    end
    r = corr(X', Y');
    YY(:,i) = Y;
    XX(:,i) = X;
    axes(ax(i)); hold on;
    plot(X, Y,'.')
    text(0.1, 0.1, sprintf('r = %.2g', r), 'units', 'normalized')
    title(nm{i} ,'interpreter', 'none', 'fontweight', 'normal')
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--')
    xlabel('sim')
    ylabel('fit')
end
set(ax, 'xtick', [], 'ytick', [])
saveFigurePdf(gcf, '~/Desktop/DDM_parameterRecovery')
%
% plot([par.A1_0], [sim.A1_0],'.')
% plot([par.A1_dR], [sim.A1_dR],'.')
% plot([par.A1_dI], [sim.A1_dI],'.')
%
% plot([par.A6_0], [sim.A6_0],'.')
% plot([par.A6_dR], [sim.A6_dR],'.')
% plot([par.A6_dI], [sim.A6_dI],'.')
%
% plot([par.z1_0], [sim.z1_0],'.')
% plot([par.z1_dR], [sim.z1_dR],'.')
% plot([par.z1_dI], [sim.z1_dI],'.')
%
% plot([par.z6_0], [sim.z6_0],'.')
% plot([par.z6_dR], [sim.z6_dR],'.')
% plot([par.z6_dI], [sim.z6_dI],'.')
%
% plot([par.x1_0], [sim.x1_0],'.')
% plot([par.x1_dR], [sim.x1_dR],'.')
% plot([par.x1_dI], [sim.x1_dI],'.')
%
% plot([par.x6_0], [sim.x6_0],'.')
% plot([par.x6_dR], [sim.x6_dR],'.')
% plot([par.x6_dI], [sim.x6_dI],'.')
%
% plot([par.T01], [sim.T01],'.')
% plot([par.T06], [sim.T06],'.')


% lsline

%% ========================================================================
%% DDM example %% DDM example %% DDM example %% DDM example %%
%% DDM example %% DDM example %% DDM example %% DDM example %%
%% DDM example %% DDM example %% DDM example %% DDM example %%
%% ========================================================================

%% FIGURE 1 - DDM example
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


%% ========================================================================
%% theta-sigma model %% theta-sigma model %% theta-sigma model %%
%% theta-sigma model %% theta-sigma model %% theta-sigma model %%
%% theta-sigma model %% theta-sigma model %% theta-sigma model %%
%% ========================================================================

%% MLE fits
priorFlag = 1;
for sn = 1:length(sub)
    fit(sn) = fit_biasNoiseBonus_v2(sub(sn), priorFlag)
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

saveFigurePdf(gcf, '~/Desktop/DDM_choices')






%% ========================================================================
%% reaction times %% reaction times %% reaction times %% reaction times %%
%% reaction times %% reaction times %% reaction times %% reaction times %%
%% reaction times %% reaction times %% reaction times %% reaction times %%
%% ========================================================================

%% FIGURE 4 - RTs and RT curves
figure(1); clf;
set(gcf, 'position', [711   575   630   700])
ax = easy_gridOfEqualFigures([0.48 0.14 0.06], [0.15 0.15 0.08]);
ax(5:7) = easy_gridOfEqualFigures([0.08 0.72], [0.1 0.13 0.13 0.03]);
[l, leg] = plot_RToverTime_v1(ax(1:2), sub, 0.1, 3);
set(leg, 'position', [0.3168    0.8719    0.1595    0.0564]);

binEdges = [-25:10:25];
e = plot_RTvsdR_v1(ax(3:4), sub, binEdges, 0.1, 3);
plot_RTregression_v1(ax(5:7), sub, 0.1, 3)

axes(ax(1)); t = title('unequal information')
axes(ax(2)); t(2) = title('equal information')
axes(ax(5)); t(3) = title('baseline RT'); ylabel('\beta_0')
axes(ax(6)); t(4) = title('effect of |\DeltaR|'); ylabel('\beta_R')
axes(ax(7)); t(5) = title('effect of \DeltaI'); ylabel('\beta_I')

set(t, 'units', 'normalized');
get(t(3), 'position')
set(t(3:5), 'position', [ 0.5000 1.1 0])


set(ax, 'tickdir', 'out', 'fontsize', 16)
set(t, 'fontweight', 'normal', 'fontsize', 20)
addABCs(ax(1:4), [-0.09 0.05], 28)
addABCs(ax(5:7), [-0.075 0.06], 28, 'EFG')
saveFigurePdf(gcf, '~/Desktop/DDM_RTs')
% addABCs(ax(1), [-0.08 0.1], 28)
% addABCs(ax(2), [-0.1 0.1], 28, 'B')
% addABCs(ax(3), [-0.09 0.1], 28, 'C')


%% ========================================================================
%% DDM parameters %% DDM parameters %% DDM parameters %% DDM parameters %%
%% DDM parameters %% DDM parameters %% DDM parameters %% DDM parameters %%
%% DDM parameters %% DDM parameters %% DDM parameters %% DDM parameters %%
%% ========================================================================

%% FIGURE 5 - plot DDM parameters vs horizon + difference
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

for i = 1:length(vn1)
    for sn = 1:length(sim)
        X1(i,sn) = getfield(sim(sn), vn1{i});
        X6(i,sn) = getfield(sim(sn), vn6{i});
    end
end

figure(1); clf;
set(gcf, 'position', [811   176   600   650])
hg = ones(5,1)*0.1;
wg = ones(4,1)*0.07;
wg(1) = 0.34;
wg(end) = 0.02;
hg(1) = 0.07;
hg(end) = 0.09;
[~,hb,wb,ax] = easy_gridOfEqualFigures(hg, wg);
ax = ax';
ax = ax(:);
thresh = 0.01;
thresh2 = 0.001;
thresh3 = 0.0001;
for i = 1:size(X1,1)
    axes(ax(i)); hold on;
    R = (rand(1,size(X1,2))-0.5)*0.25;
    plot([1+R; 2+R], [X1(i,:); X6(i,:)], 'linewidth', 1, 'color', [1 1 1]*0.75)
    plot([1+R], [X1(i,:)], 'o', 'color', AZblue, 'linewidth', 1, 'markersize', 5)
    plot([2+R], [X6(i,:)], 'o', 'color', AZred, 'linewidth', 1, 'markersize', 5)
    plot([3+R], X6(i,:)-X1(i,:), 'o', 'color', AZsand, 'linewidth', 1, 'markersize', 5)
    %plot(X1(i,:), X6(i,:), 'o', 'color', AZred, 'linewidth', 1, 'markersize', 5)
    %xl = get(ax(i), 'xlim');
    %plot(xl, xl, 'k--', 'linewidth', 1)
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
    
    if p < thresh
        str = '*';
        if p < thresh2
            str = '**';
            if p < thresh3
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
        
        %mx = max(mx1, mx6)*1.1;
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



%% ========================================================================
%% compare DDM to data %% compare DDM to data %% compare DDM to data %%
%% compare DDM to data %% compare DDM to data %% compare DDM to data %%
%% compare DDM to data %% compare DDM to data %% compare DDM to data %%
%% ========================================================================

%% FIGURE 6 - choices and RTs, model and subjects
binEdges = [-25:10:25];
figure(1); clf;
set(gcf, 'position', [811   575   600   500])
ax = easy_gridOfEqualFigures([0.15 0.2 0.05], [0.12 0.15 0.05]);

% choices -----------------------------------------------------------------
e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
e2 = plot_choiceCurves_sim_v3(ax(1:2), sim, binEdges, 0);
set(e1, 'linestyle', 'none', 'markersize', 30)
set(e2, 'marker', 'none')

set([e1([2 4]) e2([2 4])], 'visible', 'off')
e3 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);
e4 = plot_choiceCurves_sim_v3(ax(1:2), sim, binEdges, 0);
set([e3([1 3]) e4([1 3])], 'visible', 'off')
set(e3, 'linestyle', 'none', 'markersize', 30)
set(e4, 'marker', 'none')
leg = legend([e1([1]) e3(2) e2(1) e4(2)], ...
    {'human h = 1' 'human h = 6' 'model h = 1' 'model h = 6'}, ...
    'location', 'southeast');
set(leg, 'position', [0.3216    0.6614    0.1867    0.1310]);


% RTs ---------------------------------------------------------------------
e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1 , 3)
e2 = plot_rtCurves_sim_v2(ax(3:4), sim, 1, 0)
set(e1, 'linestyle', 'none', 'markersize', 30)
set(e2, 'marker', 'none')

set([e1([2 4]) e2([2 4])], 'visible', 'off')

e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1, 3)
e2 = plot_rtCurves_sim_v2(ax(3:4), sim, 1, 0)
set([e1([1 3]) e2([1 3])], 'visible', 'off')
set(e1, 'linestyle', 'none', 'markersize', 30)
set(e2, 'marker', 'none')

set(ax(3:4), 'ylim', [0.3 1.2])


axes(ax(1)); title('unequal [1 3], horizon 1')
axes(ax(2)); title('equal [2 2], horizon 1')
axes(ax(3)); title('unequal [1 3], horizon 6'); ylabel({'reaction time' '[seconds]'})
axes(ax(4)); title('equal [2 2], horizon 6'); ylabel({'reaction time' '[seconds]'})
saveFigurePdf(gcf, '~/Desktop/ddm_choiceRTmodelHuman')




%% ========================================================================
%% sensitivity analysis %% sensitivity analysis %% sensitivity analysis %%
%% sensitivity analysis %% sensitivity analysis %% sensitivity analysis %%
%% sensitivity analysis %% sensitivity analysis %% sensitivity analysis %%
%% ========================================================================

%% compare DDM to logistic model noise

for sn = 1:length(sim)
    ns_1_ddm(sn) = 1./(sim(sn).A1_dR * sim(sn).z1_0 * 2*sqrt(2));
    ns_6_ddm(sn) = 1./(sim(sn).A6_dR * sim(sn).z6_0 * 2*sqrt(2));
    ns13_1_log(sn) = fit(sn).x(2,1);
    ns13_6_log(sn) = fit(sn).x(2,2);
    ns22_1_log(sn) = fit(sn).x(5,1);
    ns22_6_log(sn) = fit(sn).x(5,2);
end

%%
% exclude the subject with negative ddm noise
idx = [sim.A1_dR]>0;
figure(1); clf;
set(gcf, 'position', [663   530   600  300]);
ax = easy_gridOfEqualFigures([0.26 0.1], [0.14 0.19 0.05]);
axes(ax(1)); hold on;
l = plot(ns22_1_log(idx), ns_1_ddm(idx), 'o');
l(2) = plot(ns22_6_log(idx), ns_6_ddm(idx), 'o');

[r1,p1] = corr(ns22_1_log(idx)', ns_1_ddm(idx)', 'type', 'spearman');
[r2,p2] = corr(ns22_6_log(idx)', ns_6_ddm(idx)', 'type', 'spearman');
plot([1 100], [1 100],'k--', 'linewidth', 1)
xlabel({'logistic model' 'noise std' })
ylabel({'DDM approximation' 'noise std' })
legend({'horizon 1' 'horizon 6'}, 'location', 'northwest')
t = title('unequal condition')
% text(0.1, 0.5, sprintf('r(%d) = %.2g\np = %.2g', sum(idx)-2, r1, p1), ...
%     'units', 'normalized', 'color', AZblue)
txt = text(0.03, 0.19, sprintf('r(%d) = %.2g', sum(idx)-2, r1), ...
    'units', 'normalized', 'color', AZblue);
txt(2) = text(0.03, 0.07, sprintf('r(%d) = %.2g', sum(idx)-2, r2), ...
    'units', 'normalized', 'color', AZred);


axes(ax(2)); hold on;
l(3) = plot(ns13_1_log(idx), ns_1_ddm(idx), 'o');
l(4) = plot(ns13_6_log(idx), ns_6_ddm(idx), 'o');
[r1,p1] = corr(ns13_1_log(idx)', ns_1_ddm(idx)', 'type', 'spearman');
[r2,p2] = corr(ns13_6_log(idx)', ns_6_ddm(idx)', 'type', 'spearman');

plot([1 100], [1 100],'k--', 'linewidth', 1)
xlabel({'logistic model' 'noise std' })
ylabel({'DDM approximation' 'noise std' })
t(2) = title('equal condition')
txt(3) = text(0.03, 0.19, sprintf('r(%d) = %.2g', sum(idx)-2, r1), ...
    'units', 'normalized', 'color', AZblue);
txt(4) = text(0.03, 0.07, sprintf('r(%d) = %.2g', sum(idx)-2, r2), ...
    'units', 'normalized', 'color', AZred);

set(l, 'markersize', 5, 'linewidth', 1)
set(l([1 3]), 'color', AZblue)
set(l([2 4]), 'color', AZred)

set(txt, 'fontsize', 16)

set(ax, 'xscale', 'log', 'yscale', 'log', 'tickdir', 'out', 'fontsize', 16, ...
    'xtick', 10.^[-6:2:2])
set(t, 'fontweight', 'normal', 'fontsize', 20)

addABCs(ax, [-0.11 0.1], 28)

% saveFigurePdf(gcf, '~/Desktop/DDM_noiseLogisticVsDDM')
