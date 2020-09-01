%% maximum likelihood approach

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

% princeton subjects only
i1 = strcmp({sub.expt_name}, 'pilot-v1');% 'repeater-v1'})
i2 = strcmp({sub.expt_name}, 'repeater-v1');
% i3 = strcmp({sub.expt_name}, 'short-pilot')
sub = sub(i1 | i2 );
% 'repeater-v1'    'short-pilot'


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




%% 
sn = 1;

rt = sub(sn).RT(:,5);
ch = sub(sn).a(:,5) == 1;
ind = (rt > 0.1) & (rt < 3);
ch = ch(ind); rt = rt(ind);

A = 0.1; 
T0 = 0;
x0 = 0;
z = 1;
c = 1;

LL = compute_lik(rt, ch, A, T0, x0, z, c)

obFunc = @(x) -compute_lik(rt, ch, x(1), x(2), x(3), x(4), c);

% linear constraints to make abs(x0) < z
A = [0 0 1 -1; 0 0 -1 -1];
B = [0 0]';

X0 = [   1   0    0   1 ]';
LB = [   0   0 -inf   0 ]';
UB = [ inf min(rt)  inf 10 ]';
X = fmincon(obFunc, X0, A, B, [], [], LB, UB);



%%
rt = sub(sn).RT(:,5);
ch = sub(sn).a(:,5) == 1;
ind = (rt > 0.1) & (rt < 3);
ch = ch(ind); rt = rt(ind);
dR = sub(sn).o2(ind,4) - sub(sn).o1(ind,4);
dI = (sub(sn).n2(ind,4) - sub(sn).n1(ind,4))/2;

cA_0 = 0;
%cA_R = 0.01;
%cA_I = 0;
% cZ_0 = 1;
cZ_R = 0;
cZ_I = 0;
cX_0 = 0;
cX_R = 0; 
cX_I = 0;
%T0 = 0.01;
c = 1;

obFunc = @(x) lik_DDMregression_v1(rt, ch, dR, dI, ...
    cA_0, x(1), x(2), ...
    x(3), cZ_R, cZ_I, ...
    cX_0, cX_R, cX_I, ...
    x(4), c);

X0 = [ 0.01 0 1 0.01 ];
LB = [ -1   -1 0 0];
UB = [ 1     1 3 2];

Xfit = fmincon(obFunc, X0, [], [], [], [], LB, UB)



%% MLE fit regression DDM 
% M for Map - which variables do we want to include in the fit?
names  = { 'cA_0' 'cA_R' 'cA_I' 'cZ_0' 'cZ_R' 'cZ_I' 'cX_0' 'cX_R' 'cX_I'  'T0'    };
%           cA_0 | cA_R | cA_I | cZ_0 | cZ_R | cZ_I | cX_0 | cX_R | cX_I |  T0
M      = [    1      1      1      1      1      1      1      1      1      1     ];

X0_all = [    0      0.01   0      1      0      0      0      0      0      0.05  ];
LB_all = [   -1    -10    -10      0     -3     -3     -1     -1     -1      0     ];
UB_all = [    1     10     10     10      3      3      1      1      1      3     ];


% initial parameters and bounds
X0 = X0_all(M==1);
LB = LB_all(M==1);
UB = UB_all(M==1);

clear Xfit
gl_vals = [5 10];
for sn = 1:length(sub)
    for i = 1:length(gl_vals)
        
        rt = sub(sn).RT(:,5);
        ch = sub(sn).a(:,5) == 2;
        gl = sub(sn).gameLength;
        ind = (rt > 0.1) & (rt < 3) & (gl == gl_vals(i));
        ch = ch(ind); rt = rt(ind);
        dR = sub(sn).o2(ind,4) - sub(sn).o1(ind,4);
        dI = (sub(sn).n2(ind,4) - sub(sn).n1(ind,4))/2;
        
        obFunc = @(x) lik_DDMregX_v1(rt, ch, dR, dI, 1, M, x);
        
        Xfit(:,sn,i) = fmincon(obFunc, X0, [], [], [], [], LB, UB);
        
    end
end

%%

% i = 5;
% figure(1); clf; hold on;
% plot(Xfit(i,:,1),'.-')
% plot(Xfit(i,:,2),'.-')

%% 
RTmin = 0.1;
RTmax = 3;
fit = fit_MLE_DDM_v2(sub, M, RTmin, RTmax, X0_all, LB_all, UB_all);

%% compare MLE fits with Sam's fits
for sn = 1:length(sim)
    sim(sn).XXfit(1,1) = sim(sn).A1_0;
    sim(sn).XXfit(1,2) = sim(sn).A6_0;
    sim(sn).XXfit(2,1) = sim(sn).A1_dR;
    sim(sn).XXfit(2,2) = sim(sn).A6_dR;
    sim(sn).XXfit(3,1) = sim(sn).A1_dI;
    sim(sn).XXfit(3,2) = sim(sn).A6_dI;
    
    sim(sn).XXfit(4,1) = sim(sn).z1_0;
    sim(sn).XXfit(4,2) = sim(sn).z6_0;
    sim(sn).XXfit(5,1) = sim(sn).z1_dR;
    sim(sn).XXfit(5,2) = sim(sn).z6_dR;
    sim(sn).XXfit(6,1) = sim(sn).z1_dI;
    sim(sn).XXfit(6,2) = sim(sn).z6_dI;
    
    sim(sn).XXfit(7,1) = sim(sn).x1_0;
    sim(sn).XXfit(7,2) = sim(sn).x6_0;
    sim(sn).XXfit(8,1) = sim(sn).x1_dR;
    sim(sn).XXfit(8,2) = sim(sn).x6_dR;
    sim(sn).XXfit(9,1) = sim(sn).x1_dI;
    sim(sn).XXfit(9,2) = sim(sn).x6_dI;
    
    sim(sn).XXfit(10,1) = sim(sn).T01;
    sim(sn).XXfit(10,2) = sim(sn).T06;
end

clear Xsim Xfit
for sn = 1:length(sim)
    Xsam(:,sn) = sim(sn).XXfit(:);
    Xmle(:,sn) = fit(sn).XXfit(:);
end

%%

figure(1); clf
hg = ones(5,1)*0.03;
wg = ones(7,1)*0.03;
ax = easy_gridOfEqualFigures(hg, wg);
for i = 1:size(Xsam,1)
    axes(ax(i)); hold on;
    plot(Xsam(i,:), Xmle(i,:),'.')
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--')
end



%% simulate an M model

clear A z b x0

sn = 1;

% time step
dt = 0.01;

% simulation parameters
XX = fit(sn).XXfit;

% experiment parameters
gameLength = fit(sn).gameLength;
dR = fit(sn).dR;
dI = fit(sn).dI;


% %%
% % unpack XX
% cA_0 = XX(1,:);
% cA_R = XX(2,:);
% cA_I = XX(3,:);
% 
% cZ_0 = XX(4,:);
% cZ_R = XX(5,:);
% cZ_I = XX(6,:);
% 
% cX_0 = XX(7,:);
% cX_R = XX(8,:);
% cX_I = XX(9,:);
% 
% T0_0   = XX(10,:);
% 
% gl_vals = [5 10];
% 
% % get a drift rate, threshold and starting point for each trial
% for i = 1:length(dR)
%     j = find(gameLength(i) == gl_vals);
%     
%     A(i) = cA_0(j) + cA_R(j) * dR(i) + cA_I(j) * dI(i);
%     z(i) = cZ_0(j) + cZ_R(j) * dR(i) + cZ_I(j) * dI(i);
%     b(i) = cX_0(j) + cX_R(j) * dR(i) + cX_I(j) * dI(i);
%     x0(i) = z(i)*(2*(1 ./ ( 1 + exp(-b(i)) )) - 1);
%     T0(i) = T0_0(j);
% end
% 
% % simulate each trial
% for i = 1:length(dR)
%     [~, ~, DT(i), C(i)] = simluate_DDM_v1(dt, A(i), 1, z(i), x0(i));
%     RT(i) = DT(i) + T0(i);
% end
% 
% fak.sID(1:length(dR)) = nan;
% fak.XX = XX;
% fak.dR = dR;
% fak.dI = dI;
% fak.gameLength = gameLength;
% fak.choice = C + 1;
% fak.RT = RT;

%%
clear fak
dt = 0.001;
for sn = 1:length(fit)
    fak(sn) = simulate_Mmodel_v1(fit(sn).XXfit, dt, fit(sn).dR, fit(sn).dI, fit(sn).gameLength);
end
%% fit simulated data
fit_fak = fit_MLE_DDM_v2(fak, M, RTmin, RTmax, X0_all, LB_all, UB_all);

% fit parameters to real data
% simulate data with fit parameters - make fak
% fit fak for parameter recovery

%%
clear Xsim Xfit
for sn = 1:length(fak)
    Xsim(:,sn) = fak(sn).XX(:);
    Xfit(:,sn) = fit_fak(sn).XXfit(:);
end


figure(1); clf
hg = ones(5,1)*0.05;
wg = ones(7,1)*0.05;
ax = easy_gridOfEqualFigures(hg, wg);
for i = 1:size(Xsim,1)
    axes(ax(i)); hold on;
    plot(Xsim(i,:), Xfit(i,:),'.')
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--')
    xlabel('sim')
    ylabel('fit')
end
set(ax, 'xtick', [], 'ytick', [])
saveFigurePdf(gcf, '~/Desktop/DDM_mleParameterRecovery')



