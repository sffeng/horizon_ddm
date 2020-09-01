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


%% load Sam's new parameters
fname1 = 'fittedparams_nt131filtered.csv';
fname2 = 'fittedT0_nt131filtered.csv';

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
sub = sub(idx);

for sn = 1:length(sim)
    sim(sn).dR = sub(sn).o2(:,4) - sub(sn).o1(:,4);
    sim(sn).dI = (sub(sn).n2(:,4) - sub(sn).n1(:,4)) / 2;
    %sim(sn).RT = sub(sn).RT(:,5);
    sim(sn).gameLength = sub(sn).gameLength;
    sim(sn).h = sub(sn).gameLength-4;
end

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
%%
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

%% OLD
% 
% z    = sim(sn).z1_0 + sim(sn).z1_dR * dR + sim(sn).z1_dI * dI;
% bias = sim(sn).x1_0 + sim(sn).x1_dR * dR + sim(sn).x1_dI * dI;
% A    = sim(sn).A1_0 + sim(sn).A1_dR * dR + sim(sn).A1_dI * dI
% T0   = sim(sn).T01;
% 
% x0   = 1 ./ (1 + exp(-bias)) .* z;
% 
% a = A.^2;
% z = z./A;
% x0 = x0./A;
% 
% RT = T0 + z .* tanh(z .* a) + ( 2*z.*(1-exp(-2*x0.*a))./( exp(2*z.*a) - exp(-2*z.*a) ) - x0  );
% ER = 1 ./ ( 1 + exp(2*z.*a) ) - ( (1-exp(-2*x0.*a))./(exp(2*z.*a) - exp(-2*z.*a) ) );
% % 
% r_vals = [-30:0.01:30];
% for sn = 1:length(sim)
%     
%     [ER, RT] = compute_contERRT_v1(...
%         r_vals, ...
%         sim(sn).z1_0, sim(sn).z1_dR, sim(sn).z1_dI, ...
%         sim(sn).x1_0, sim(sn).x1_dR, sim(sn).x1_dI, ...
%         sim(sn).A1_0, sim(sn).A1_dR, sim(sn).A1_dI, ...
%         sim(sn).T01);
%     sim(sn).cont_RT1_13 = RT;
%     sim(sn).cont_ER1_13 = ER;
% 
%     [ER, RT] = compute_contERRT_v1(...
%         r_vals, ...
%         sim(sn).z1_0, sim(sn).z1_dR, 0, ...
%         sim(sn).x1_0, sim(sn).x1_dR, 0, ...
%         sim(sn).A1_0, sim(sn).A1_dR, 0, ...
%         sim(sn).T01);
%     sim(sn).cont_RT1_22 = RT;
%     sim(sn).cont_ER1_22 = ER;
% 
%     [ER, RT] = compute_contERRT_v1(...
%         r_vals, ...
%         sim(sn).z6_0, sim(sn).z6_dR, sim(sn).z6_dI, ...
%         sim(sn).x6_0, sim(sn).x6_dR, sim(sn).x6_dI, ...
%         sim(sn).A6_0, sim(sn).A6_dR, sim(sn).A6_dI, ...
%         sim(sn).T01);
%     sim(sn).cont_RT6_13 = RT;
%     sim(sn).cont_ER6_13 = ER;
%     
%     [ER, RT] = compute_contERRT_v1(...
%         r_vals, ...
%         sim(sn).z6_0, sim(sn).z6_dR, 0, ...
%         sim(sn).x6_0, sim(sn).x6_dR, 0, ...
%         sim(sn).A6_0, sim(sn).A6_dR, 0, ...
%         sim(sn).T01);
%     sim(sn).cont_RT6_22 = RT;
%     sim(sn).cont_ER6_22 = ER;
%     
% end




%% load simulations
fname = 'FUllData.csv';
fid = fopen(fname);

hdr = textscan(fid,...
    '%s%s%s%s%s%s%s%s%s', ...
    1, 'delimiter', ',');
data = textscan(fid,...
    '%f%f%f%f%f%f%f%f%f', ...
    'delimiter', ',');
fclose(fid);

%% organize into subjects
uu = unique(data{1});
for sn = 1:length(uu)
    ind = uu(sn) == (data{1});
    sim2(sn).dR = data{3}(ind);
    sim2(sn).dI = -data{4}(ind);
    sim2(sn).RT = data{5}(ind);
    sim2(sn).h = data{2}(ind);
    
    % DDM parameters x0	z	A	T0
    sim2(sn).x0 = data{6}(ind);
    sim2(sn).z = data{7}(ind);
    sim2(sn).A = data{8}(ind);
    sim2(sn).T0 = data{9}(ind);
    
end

%% compute mean RT from DDM parameters
for sn = 1:length(sim2)
    z = sim2(sn).z;
    x0 = sim2(sn).x0.*z;
    A = sim2(sn).A;
    T0 = sim2(sn).T0;
    
    a = A.^2;
    z = z./A;
    x0 = x0./A;
    
    RT = T0 + z .* tanh(z .* a) + ( 2*z.*(1-exp(-2*x0.*a))./( exp(2*z.*a) - exp(-2*z.*a) ) - x0  );
    ER = 1 ./ ( 1 + exp(2*z.*a) ) - ( (1-exp(-2*x0.*a))./(exp(2*z.*a) - exp(-2*z.*a) ) );
    sim2(sn).RTana = RT;
    sim2(sn).ER = ER;
end
% %%
% L = load('sim2Data.mat');
% X = L.data;
%
% % organize into subjects
% uu = unique(X.Subject);
% for sn = 1:length(uu)
%     ind = uu(sn) == X.Subject;
%     sim2(sn).dR = X.dR(ind);
%     sim2(sn).dI = -X.dI(ind);
%     sim2(sn).RT = X.rt(ind);
%     sim2(sn).h = X.Horizon(ind);
% end
%


%% match sim2 to sub
clear isEq
for i1 = 1:length(sub)
    for i2 = 1:length(sim2)
        hh = sub(i1).gameLength;
        hh(hh==5) = 1;
        hh(hh==10) = 6;
        isEq(i1,i2) = isequal(sim2(i2).h, hh);
    end
end
sub = sub((sum(isEq,2)~=0));

%% back out Sam's original parameters
for sn = 1:length(sim2)
    
    i1 = sim2(sn).h == 1;
    i6 = sim2(sn).h == 6;
    
    T0 = sim2(sn).T0;
    x0 = sim2(sn).x0;
    l0 = (x0 + 1)/2;
    l0 = log((l0)./(1-l0));
    
    z = sim2(sn).z;
    A = sim2(sn).A;
    dR = sim2(sn).dR;
    dI = sim2(sn).dI;
    
    B = glmfit([dR(i1) dI(i1)], A(i1)); A1_0 = B(1); A1_dR = B(2); A1_dI = B(3);
    
    B = glmfit([dR(i6) dI(i6)], A(i6)); A6_0 = B(1); A6_dR = B(2); A6_dI = B(3);
    
    B = glmfit([dR(i1) dI(i1)], z(i1)); z1_0 = B(1); z1_dR = B(2); z1_dI = B(3);
    B = glmfit([dR(i6) dI(i6)], z(i6)); z6_0 = B(1); z6_dR = B(2); z6_dI = B(3);
    
    B = glmfit([dR(i1) dI(i1)], l0(i1)); x1_0 = B(1); x1_dR = B(2); x1_dI = B(3);
    B = glmfit([dR(i6) dI(i6)], l0(i6)); x6_0 = B(1); x6_dR = B(2); x6_dI = B(3);
    
    %A2 = A_0 + A_dR * dR(i1) + A_dI * dI(i1);
    l2 = x6_0 + x6_dR * dR(i6) + x6_dI * dI(i6);
    
    sim2(sn).T01 = unique(T0(i1));
    sim2(sn).T06 = unique(T0(i6));
    
    sim2(sn).A1_0  = A1_0;
    sim2(sn).A1_dR = A1_dR;
    sim2(sn).A1_dI = A1_dI;
    sim2(sn).A6_0  = A6_0;
    sim2(sn).A6_dR = A6_dR;
    sim2(sn).A6_dI = A6_dI;
    
    sim2(sn).z1_0  = z1_0;
    sim2(sn).z1_dR = z1_dR;
    sim2(sn).z1_dI = z1_dI;
    sim2(sn).z6_0  = z6_0;
    sim2(sn).z6_dR = z6_dR;
    sim2(sn).z6_dI = z6_dI;
    
    sim2(sn).x1_0  = x1_0;
    sim2(sn).x1_dR = x1_dR;
    sim2(sn).x1_dI = x1_dI;
    sim2(sn).x6_0  = x6_0;
    sim2(sn).x6_dR = x6_dR;
    sim2(sn).x6_dI = x6_dI;
end



%% compute parameters on every trial for sim
% for sn = 1:length(sim2)
%     for i = 1:length(sim2(sn).dR)
%         if sim2(sn).h(i) == 1
%             
%             sim2(sn).z2(i,1)    = sim2(sn).z1_0 + sim2(sn).z1_dR * sim2(sn).dR(i) + sim2(sn).z1_dI * sim2(sn).dI(i);
%             sim2(sn).bias2(i,1) = sim2(sn).x1_0 + sim2(sn).x1_dR * sim2(sn).dR(i) + sim2(sn).x1_dI * sim2(sn).dI(i);
%             sim2(sn).A2(i,1)    = sim2(sn).A1_0 + sim2(sn).A1_dR * sim2(sn).dR(i) + sim2(sn).A1_dI * sim2(sn).dI(i);
%             sim2(sn).T02(i,1)   = sim2(sn).T01;
%             
%         else
%             
%             sim2(sn).z2(i,1)    = sim2(sn).z6_0 + sim2(sn).z6_dR * sim2(sn).dR(i) + sim2(sn).z6_dI * sim2(sn).dI(i);
%             sim2(sn).bias2(i,1) = sim2(sn).x6_0 + sim2(sn).x6_dR * sim2(sn).dR(i) + sim2(sn).x6_dI * sim2(sn).dI(i);
%             sim2(sn).A2(i,1)    = sim2(sn).A6_0 + sim2(sn).A6_dR * sim2(sn).dR(i) + sim2(sn).A6_dI * sim2(sn).dI(i);
%             sim2(sn).T02(i,1)   = sim2(sn).T06;
%             
%         end
%         
%         
%         sim2(sn).x02   = 2 ./ (1 + exp(-sim2(sn).bias2)) - 1;
% 
%     end
% end


%% analytic RTs on continuous x-axis
% r_vals = [-30:0.01:30];
% for sn = 1:length(sim)
%     
%     [ER, RT] = compute_contERRT_v1(...
%         r_vals, ...
%         sim(sn).z1_0, sim(sn).z1_dR, sim(sn).z1_dI, ...
%         sim(sn).x1_0, sim(sn).x1_dR, sim(sn).x1_dI, ...
%         sim(sn).A1_0, sim(sn).A1_dR, sim(sn).A1_dI, ...
%         sim(sn).T01);
%     sim(sn).cont_RT1_13 = RT;
%     sim(sn).cont_ER1_13 = ER;
% 
%     [ER, RT] = compute_contERRT_v1(...
%         r_vals, ...
%         sim(sn).z1_0, sim(sn).z1_dR, 0, ...
%         sim(sn).x1_0, sim(sn).x1_dR, 0, ...
%         sim(sn).A1_0, sim(sn).A1_dR, 0, ...
%         sim(sn).T01);
%     sim(sn).cont_RT1_22 = RT;
%     sim(sn).cont_ER1_22 = ER;
% 
%     [ER, RT] = compute_contERRT_v1(...
%         r_vals, ...
%         sim(sn).z6_0, sim(sn).z6_dR, sim(sn).z6_dI, ...
%         sim(sn).x6_0, sim(sn).x6_dR, sim(sn).x6_dI, ...
%         sim(sn).A6_0, sim(sn).A6_dR, sim(sn).A6_dI, ...
%         sim(sn).T01);
%     sim(sn).cont_RT6_13 = RT;
%     sim(sn).cont_ER6_13 = ER;
%     
%     [ER, RT] = compute_contERRT_v1(...
%         r_vals, ...
%         sim(sn).z6_0, sim(sn).z6_dR, 0, ...
%         sim(sn).x6_0, sim(sn).x6_dR, 0, ...
%         sim(sn).A6_0, sim(sn).A6_dR, 0, ...
%         sim(sn).T01);
%     sim(sn).cont_RT6_22 = RT;
%     sim(sn).cont_ER6_22 = ER;
%     
% end


%% ========================================================================
%% DDM example %% DDM example %% DDM example %% DDM example %%
%% DDM example %% DDM example %% DDM example %% DDM example %%
%% DDM example %% DDM example %% DDM example %% DDM example %%
%% ========================================================================

%% FIGURE DDM example 
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
% remove fast or slow RTs
for sn = 1:length(sub)
    ind = (sub(sn).RT(:,5) > 0.1) & (sub(sn).RT(:,5) < 3);
    sub(sn).a(~ind,5) = nan;
end
%%
for sn = 1:length(sub)
    fit(sn) = fit_biasNoiseBonus_v2(sub(sn), priorFlag)
end

%% model choice curves
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

%% FIGURE - basic choice curves + model parameters
XX = cat(3, fit.x);

bonus = squeeze(XX(3,:,:));
noise13 = squeeze(XX(2,:,:));
noise22 = squeeze(XX(5,:,:));

figure(1); clf;
set(gcf, 'position', [811   575   550   500])
ax = easy_gridOfEqualFigures([0.66  0.07], [0.14 0.14 0.07]);
ax(3:5) = easy_gridOfEqualFigures([0.1  0.6], [0.12 0.12 0.12 0.05]);
% choice curve bit --------------------------------------------------------
axes(ax(1)); hold on;
e = errorbar(X, m_13_1, s_13_1);
e(2) = errorbar(X, m_13_6, s_13_6);

l1 = plot(x_vals, CC13);

xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
ylabel({'p(high info)'})
t = title('unequal [1 3]', 'fontweight', 'normal');
leg = legend(e([2 1]), {'horizon 6' 'horizon 1'}, 'location', 'southeast');
set(leg, 'position', [ 0.3000    0.6706    0.1827    0.0790])
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

set(e, 'linewidth', 3, 'marker', '.', 'markersize', 25)
set(e([1 3]), 'color', AZblue)
set(e([2 4]), 'color', AZred)
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

%% RT vs trial in task
clear RT
for sn = 1:length(sub)
    i22 = (sub(sn).n1(:,4) == 2);
    i13 = ~i22;
    
    i1 = sub(sn).gameLength == 5;
    i6 = sub(sn).gameLength == 10;
    
    ind = i13 & i1;
    RT(:,1,sn) = nanmean(sub(sn).RT(ind, :));
    ind = i13 & i6;
    RT(:,2,sn) = nanmean(sub(sn).RT(ind, :));
    ind = i22 & i1;
    RT(:,3,sn) = nanmean(sub(sn).RT(ind, :));
    ind = i22 & i6;
    RT(:,4,sn) = nanmean(sub(sn).RT(ind, :));
end

rt = nanmean(RT,3);
ss = nanstd(RT, [], 3)/ sqrt(length(sub));

figure(1); clf;
ax = easy_gridOfEqualFigures([0.1 0.1], [0.1 0.1 0.03]);
axes(ax(1)); hold on;
f = fill([0.5 4.5 4.5 0.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', [1 1 1]*0.85, 'linestyle', 'none')
l = errorbar(rt(:,[2:-1:1]), ss(:,1:2));
ylabel('reaction time [seconds]')
xlabel('trial number')

axes(ax(2)); hold on;
f = fill([0.5 4.5 4.5 0.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', [1 1 1]*0.85, 'linestyle', 'none')
l(3:4) = errorbar(rt(:,[4:-1:3]), ss(:,1:2));
ylabel('reaction time [seconds]')
xlabel('trial number')


set(l, 'linewidth', 3, 'marker', '.', 'markersize', 30)

set(l([1 3]), 'color', AZred)
set(l([2 4]), 'color', AZblue)

set(ax(1:2), 'ylim', [0 1.2], 'xtick', [1:10], ...
    'xticklabel', {'i1' 'i2' 'i3' 'i4' 1:6}, ...
    'tickdir', 'out', 'xlim', [0.5 10.5])

%% basic RT curves



clear M_13_1 M_13_6 M_22_1 M_22_6
binEdges = [-25:10:25];
for sn = 1:length(sim)
    RT = sim(sn).RT;
    dR = sim(sn).dR;
    i22 = sim(sn).dI == 0;
    i13 = ~i22;
    i1 = sim(sn).h == 1;
    i6 = sim(sn).h == 6;
    dI = sim(sn).dI;
    
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

figure(2); clf;
set(gcf, 'position', [811   575   600   300])
% ax = easy_gridOfEqualFigures([0.25 0.08], [0.15 0.15 0.03]);
ax = easy_gridOfEqualFigures([0.3 0.12], [0.12 0.15 0.05]);
axes(ax(1)); hold on;
e = errorbar(X, m_13_1, s_13_1);
e(2) = errorbar(X, m_13_6, s_13_6);
xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
ylabel({'RT [z-score]'})
title('unequal [1 3]', 'fontweight', 'normal')
legend({'horizon 1' 'horizon 6'}, 'location', 'south')

axes(ax(2)); hold on;
e(3) = errorbar(X, m_22_1, s_22_1);
e(4) = errorbar(X, m_22_6, s_22_6);
xlabel({'difference in mean reward' 'R(left) - R(right)'});
ylabel({'RT [z-score]'})
title('equal [2 2]', 'fontweight', 'normal')


set(ax, 'xlim', [-35 35], 'ylim', [-0.15 1.65], 'tickdir', 'out', 'fontsize', 20)
set(e, 'linewidth', 3, 'marker', '.', 'markersize', 50)
set(e([1 3]), 'color', AZblue)
set(e([2 4]), 'color', AZred)

saveFigurePdf(gcf, '~/Desktop/model')


%% FIGURE RTs and RT curves
clear RT
for sn = 1:length(sub)
    i22 = (sub(sn).n1(:,4) == 2);
    i13 = ~i22;
    
    i1 = sub(sn).gameLength == 5;
    i6 = sub(sn).gameLength == 10;
    
    %idx = sub(sn).RT(:,5) < 3;
    ind = i13 & i1;% & idx;
    RT(:,1,sn) = nanmean(sub(sn).RT(ind, :));
    ind = i13 & i6;% & idx;
    RT(:,2,sn) = nanmean(sub(sn).RT(ind, :));
    ind = i22 & i1;% & idx;
    RT(:,3,sn) = nanmean(sub(sn).RT(ind, :));
    ind = i22 & i6;% & idx;
    RT(:,4,sn) = nanmean(sub(sn).RT(ind, :));
end

rt = nanmean(RT,3);
ss = nanstd(RT, [], 3)/ sqrt(length(sub));

figure(1); clf;
set(gcf, 'position', [811   575   600   500])
ax = easy_gridOfEqualFigures([0.15 0.15 0.08], [0.13 0.15 0.03]);
axes(ax(1)); hold on;
f = fill([0.5 4.5 4.5 0.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', [1 1 1]*0.9, 'linestyle', 'none')
mix = 0.7;

f = fill([4.5 5.5 5.5 4.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', 0.1*AZsand+0.9, 'linestyle', 'none')

l = errorbar(rt(:,[2:-1:1]), ss(:,1:2));
ylabel({'reaction time' '[seconds]'})
xlabel('trial number')
legend(l(2:-1:1), {'horizon 1' 'horizon 6'}, 'location', 'northeast')
t = title('unequal information')


axes(ax(2)); hold on;
f = fill([0.5 4.5 4.5 0.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', [1 1 1]*0.9, 'linestyle', 'none')

f = fill([4.5 5.5 5.5 4.5], [0 0 1.2 1.2],'r');
set(f, 'facecolor', 0.1*AZsand+0.9, 'linestyle', 'none')
% [1 1 1]*mix+(1-mix)*[1 1 0]
l(3:4) = errorbar(rt(:,[4:-1:3]), ss(:,1:2));
ylabel({'reaction time' '[seconds]'})
xlabel('trial number')
t(2) = title('equal information')

set(l, 'linewidth', 3, 'marker', '.', 'markersize', 30)

set(l([1 3]), 'color', AZred)
set(l([2 4]), 'color', AZblue)

set(ax(1:2), 'ylim', [0 1.2], 'xtick', [1:10], ...
    'xticklabel', {'i1' 'i2' 'i3' 'i4' 1:6}, ...
    'xlim', [0.5 10.5])


% basic RT curves



clear M_13_1 M_13_6 M_22_1 M_22_6
binEdges = [-25:10:25];
for sn = 1:length(sim)
    RT = sim(sn).RT;
    RT(RT>3) = nan;
    dR = sim(sn).dR;
    i22 = sim(sn).dI == 0;
    i13 = ~i22;
    i1 = sim(sn).h == 1;
    i6 = sim(sn).h == 6;
    dI = sim(sn).dI;
    
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

axes(ax(3)); hold on;
e = errorbar(X, m_13_1, s_13_1);
e(2) = errorbar(X, m_13_6, s_13_6);
xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
ylabel({'reaction time' '[seconds]'})
% title('unequal [1 3]', 'fontweight', 'normal')
% legend({'horizon 1' 'horizon 6'}, 'location', 'south')

axes(ax(4)); hold on;
e(3) = errorbar(X, m_22_1, s_22_1);
e(4) = errorbar(X, m_22_6, s_22_6);
xlabel({'difference in mean reward' 'R(left) - R(right)'});
ylabel({'reaction time' '[seconds]'})
% title('equal [2 2]', 'fontweight', 'normal')


set(ax(3:4), 'xlim', [-35 35], 'ylim', [0 1.5])
set(e, 'linewidth', 3, 'marker', '.', 'markersize', 30)
set(e([1 3]), 'color', AZblue)
set(e([2 4]), 'color', AZred)



set(ax, 'tickdir', 'out', 'fontsize', 16)
set(t, 'fontweight', 'normal', 'fontsize', 20)
addABCs(ax, [-0.09 0.05], 28)
saveFigurePdf(gcf, '~/Desktop/DDM_RTs')

%% STATS: RTs
clear B
% glm for RTs
for sn = 1:length(sub)
    % say 1 is left
    dR = sub(sn).o1(:,4) - sub(sn).o2(:,4);
    dI = -(sub(sn).n1(:,4) - sub(sn).n2(:,4))/2;
    c = sub(sn).a(:,5);
    c(c==2) = -1;
    RT = sub(sn).RT(:,5);
    RT(RT>3) = nan;
    i1 = sub(sn).gameLength == 5;
    i6 = sub(sn).gameLength == 10;
    
    B1(:, sn) = glmfit( [c(i1).*dR(i1) c(i1).*dI(i1)], (RT(i1)));
    B6(:, sn) = glmfit( [c(i6).*dR(i6) c(i6).*dI(i6)], (RT(i6)));
end

figure(1); clf;
ax = easy_gridOfEqualFigures([0.18 0.1], [0.12 0.14 0.14 0.03]);

R = (rand(1,size(B1,2))-0.5)*0.25;
clear l1 l6
for i = 1:3
    axes(ax(i)); hold on;
    plot([1+R; 2+R], [B1(i,:); B6(i,:)], 'color', [1 1 1]*0.75, 'linewidth', 1)
    l1(i) = plot(1+R, B1(i,:),'o');
    l6(i) = plot(2+R, B6(i,:),'o');
    if i > 1
        plot([0.5 2.5], [0 0], 'k--', 'linewidth', 1)
    end
end
set(l1, 'color', AZblue)
set(l6, 'color', AZred)
set([l1 l6], 'markersize', 5, 'linewidth', 1)
set(ax(1), 'ylim', [0 2]);
set(ax(2), 'ylim', [-0.04 0.01]);
set(ax(3), 'ylim', [-0.2 0.2])
set(ax, 'xlim', [0.5 2.5], 'tickdir', 'out', 'fontsize', 16, ...
    'xtick', [1 2], 'xticklabel', [1 6])

axes(ax(1)); ylabel({'\beta_0 [seconds]'}); xlabel('horizon')
axes(ax(2)); ylabel({'\beta_R [seconds per point]'}); xlabel('horizon')
axes(ax(3)); ylabel({'\beta_I [seconds]'}); xlabel('horizon')
set(gcf, 'position', [811   517   600   300])

axes(ax(1));
plot([1 1 2 2], [1.8 1.9 1.9 1.8], 'k-', 'linewidth', 1)
text(1.5, 1.91, '*', 'fontsize', 30, 'horizontalalignment', 'center')

axes(ax(2));
plot([1 1 2 2], [0.005 0.008 0.008 0.005], 'k-', 'linewidth', 1)
text(1.5, 0.0082, '***', 'fontsize', 30, 'horizontalalignment', 'center')
addABCs(ax(1), [-0.08 0.1], 28)
addABCs(ax(2), [-0.1 0.1], 28, 'B')
addABCs(ax(3), [-0.09 0.1], 28, 'C')
saveFigurePdf(gcf, '~/Desktop/DDM_RTglm')


%% compute theoretical RT curves
dr = x_vals;
for sn = 1:length(sub)
    rtL_13_1 = B1(1,sn) + B1(2,sn)*dr + B1(3,sn);
    rtR_13_1 = B1(1,sn) - B1(2,sn)*dr - B1(3,sn);
    rt_13_1(:,sn) = rtL_13_1'.*cc13(:,1,sn) + rtR_13_1'.*(1-cc13(:,1,sn));
    rtL_13_6 = B6(1,sn) + B6(2,sn)*dr + B6(3,sn);
    rtR_13_6 = B6(1,sn) - B6(2,sn)*dr - B6(3,sn);
    rt_13_6(:,sn) = rtL_13_6'.*cc13(:,2,sn) + rtR_13_6'.*(1-cc13(:,2,sn));
end
% rtc_22_1 =

%% ========================================================================
%% DDM parameters %% DDM parameters %% DDM parameters %% DDM parameters %%
%% DDM parameters %% DDM parameters %% DDM parameters %% DDM parameters %%
%% DDM parameters %% DDM parameters %% DDM parameters %% DDM parameters %%
%% ========================================================================

%% FIGURE plot DDM parameters vs horizon
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

% vn1 = {
%     'x1_0'
%     'x1_dI'
%     'x1_dR'
%     'A1_0'
%     'A1_dI'
%     'A1_dR'
%     'T01'
%     'z1_0'
%     'z1_dI'
%     'z1_dR'
%     };
% 
% vn6 = {
%     'x6_0'
%     'x6_dI'
%     'x6_dR'
%     'A6_0'
%     'A6_dI'
%     'A6_dR'
%     'T06'
%     'z6_0'
%     'z6_dI'
%     'z6_dR'
%     };

for i = 1:length(vn1)
    for sn = 1:length(sim)
        X1(i,sn) = getfield(sim(sn), vn1{i});
        X6(i,sn) = getfield(sim(sn), vn6{i});
    end
end

figure(1); clf;
set(gcf, 'position', [811   176   600   600])
hg = ones(5,1)*0.07;
wg = ones(4,1)*0.07;
wg(1) = 0.34;
wg(end) = 0.02;
hg(1) = 0.07;
hg(end) = 0.06;
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
    %plot(X1(i,:), X6(i,:), 'o', 'color', AZred, 'linewidth', 1, 'markersize', 5)
    [~,p, ~, stats] = ttest([X1(i,:)-X6(i,:)]);
    %xl = get(ax(i), 'xlim');
    %plot(xl, xl, 'k--', 'linewidth', 1)
    if i == size(X1,1)
        xlabel('horizon')
    end
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
a(2) = annotation('textbox', [0 hg(1)+hg(2)+hb(1) wg(1)-delta hb(1)], 'string', 'starting point, x_0');
a(3) = annotation('textbox', [0 sum(hg(1:3))+sum(hb(1:2)) wg(1)-delta hb(1)], 'string', 'threshold, z');
a(4) = annotation('textbox', [0 sum(hg(1:4))+sum(hb(1:3)) wg(1)-delta hb(1)], 'string', 'drift rate, A');
set(a, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'right', ...
    'verticalalignment', 'middle')

delta2 = 0.02;
b = annotation('textbox', [wg(1) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'baseline');
b(2) = annotation('textbox', [wg(1)+wb(1)+wg(2) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of dR');
b(3) = annotation('textbox', [wg(1)+wb(1)+wg(2)+wb(2)+wg(3) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of dI');
set(b, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'center', ...
    'verticalalignment', 'middle')

set(ax, 'xtick', [1 2], 'xticklabel', [1 6], 'xlim', [0.5 2.5], ...
    'tickdir', 'out')
set(gcf, 'InvertHardcopy', 'off')
saveFigurePdf(gcf, '~/Desktop/DDM_params')

%% FIGURE plot DDM parameters vs horizon + difference
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

% vn1 = {
%     'x1_0'
%     'x1_dI'
%     'x1_dR'
%     'A1_0'
%     'A1_dI'
%     'A1_dR'
%     'T01'
%     'z1_0'
%     'z1_dI'
%     'z1_dR'
%     };
% 
% vn6 = {
%     'x6_0'
%     'x6_dI'
%     'x6_dR'
%     'A6_0'
%     'A6_dI'
%     'A6_dR'
%     'T06'
%     'z6_0'
%     'z6_dI'
%     'z6_dR'
%     };

for i = 1:length(vn1)
    for sn = 1:length(sim)
        X1(i,sn) = getfield(sim(sn), vn1{i});
        X6(i,sn) = getfield(sim(sn), vn6{i});
    end
end

figure(1); clf;
set(gcf, 'position', [811   176   600   600])
hg = ones(5,1)*0.07;
wg = ones(4,1)*0.07;
wg(1) = 0.34;
wg(end) = 0.02;
hg(1) = 0.07;
hg(end) = 0.06;
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
    [~,p, ~, stats] = ttest([X1(i,:)-X6(i,:)]);
    %xl = get(ax(i), 'xlim');
    %plot(xl, xl, 'k--', 'linewidth', 1)
    if i == size(X1,1)
        xlabel('horizon')
    end
    plot([0.5 3.5], [0 0], 'k--', 'linewidth', 1)
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
a(2) = annotation('textbox', [0 hg(1)+hg(2)+hb(1) wg(1)-delta hb(1)], 'string', 'starting point, x_0');
a(3) = annotation('textbox', [0 sum(hg(1:3))+sum(hb(1:2)) wg(1)-delta hb(1)], 'string', 'threshold, z');
a(4) = annotation('textbox', [0 sum(hg(1:4))+sum(hb(1:3)) wg(1)-delta hb(1)], 'string', 'drift rate, A');
set(a, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'right', ...
    'verticalalignment', 'middle')

delta2 = 0.02;
b = annotation('textbox', [wg(1) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'baseline');
b(2) = annotation('textbox', [wg(1)+wb(1)+wg(2) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of dR');
b(3) = annotation('textbox', [wg(1)+wb(1)+wg(2)+wb(2)+wg(3) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of dI');
set(b, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'center', ...
    'verticalalignment', 'middle')

set(ax, 'xtick', [1 2 3], 'xticklabel', {1 6 'diff'}, 'xlim', [0.5 3.5], ...
    'tickdir', 'out')
set(gcf, 'InvertHardcopy', 'off')
saveFigurePdf(gcf, '~/Desktop/DDM_params2')



%% ========================================================================
%% compare DDM to data %% compare DDM to data %% compare DDM to data %%
%% compare DDM to data %% compare DDM to data %% compare DDM to data %%
%% compare DDM to data %% compare DDM to data %% compare DDM to data %%
%% ========================================================================

%% FIGURE - choices and RTs, model and subjects
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


%% FIGURE - Choices + RTs, just model 
binEdges = [-25:10:25];
figure(1); clf;
set(gcf, 'position', [811   575   600   500])
% ax = easy_gridOfEqualFigures([0.25 0.08], [0.15 0.15 0.03]);
ax = easy_gridOfEqualFigures([0.15 0.2 0.05], [0.12 0.15 0.05]);
e1 = plot_choiceCurves_sim_v2(ax(1:2), sim, binEdges); %plot_rtCurves_sim_v1(ax(3:4), sim, 1)
e2 = plot_rtCurves_sim_v1(ax(3:4), sim, 1)
set([e1 e2], 'markersize', 30)
axes(ax(3)); ylabel('RT [seconds]'); title('')
axes(ax(4)); ylabel('RT [seconds]'); title('')
set(ax, 'tickdir', 'out')
axes(ax(1)); title('unequal [1 3]', 'fontsize', 18)
axes(ax(2)); title('equal [2 2]', 'fontsize', 18)
addABCs(ax, [-0.08 0.05], 28)

saveFigurePdf(gcf, '~/Desktop/DDM_modelChoiceRT')
% set([e1([2 4]) e2([2 4])], 'visible', 'off')


%% FIGURE - RTs for subjects
% clear M_13_1 M_13_6 M_22_1 M_22_6
% binEdges = [-25:10:25];
% for sn = 1:length(sub)
%     RT = sub(sn).RT(:,5);
%     %RT(RT>3) = nan;
%     dR = sub(sn).o2(:,4) - sub(sn).o1(:,4);
%     i22 = sub(sn).n2(:,4) == 2;
%     i13 = ~i22;
%     i1 = sub(sn).gameLength == 5;
%     i6 = sub(sn).gameLength == 10;
%     dI = (sub(sn).n2(:,4) - sub(sn).n1(:,4))/2;
%
%     [M_13_1(:,sn), ~, X] = binIt(dI(i13&i1).*dR(i13&i1), RT(i13&i1), binEdges, 'std');
%     [M_13_6(:,sn), ~, X] = binIt(dI(i13&i6).*dR(i13&i6), RT(i13&i6), binEdges, 'std');
%     [M_22_1(:,sn), ~, X] = binIt(dR(i22&i1), RT(i22&i1), binEdges, 'std');
%     [M_22_6(:,sn), ~, X] = binIt(dR(i22&i6), RT(i22&i6), binEdges, 'std');
% end
%
%
% m_13_1 = nanmean(M_13_1,2);
% s_13_1 = nanstd(M_13_1,[],2)/sqrt(length(sub));
% m_13_6 = nanmean(M_13_6,2);
% s_13_6 = nanstd(M_13_6,[],2)/sqrt(length(sub));
% m_22_1 = nanmean(M_22_1,2);
% s_22_1 = nanstd(M_22_1,[],2)/sqrt(length(sub));
% m_22_6 = nanmean(M_22_6,2);
% s_22_6 = nanstd(M_22_6,[],2)/sqrt(length(sub));
%
% figure(1); clf;
% set(gcf, 'position', [811   575   600   300])
% % ax = easy_gridOfEqualFigures([0.25 0.08], [0.15 0.15 0.03]);
% ax = easy_gridOfEqualFigures([0.3 0.12], [0.12 0.15 0.05]);
% axes(ax(1)); hold on;
% e = errorbar(X, m_13_1, s_13_1);
% e(2) = errorbar(X, m_13_6, s_13_6);
% xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
% ylabel({'RT [z-score]'})
% title('unequal [1 3]', 'fontweight', 'normal')
% legend({'horizon 1' 'horizon 6'}, 'location', 'northeast')
%
% axes(ax(2)); hold on;
% e(3) = errorbar(X, m_22_1, s_22_1);
% e(4) = errorbar(X, m_22_6, s_22_6);
% xlabel({'difference in mean reward' 'R(left) - R(right)'});
% ylabel({'RT [z-score]'})
% title('equal [2 2]', 'fontweight', 'normal')
%
%
% set(ax, 'xlim', [-35 35], 'ylim', [-0.15 1.65], 'tickdir', 'out', 'fontsize', 20)
% set(e, 'linewidth', 3, 'marker', '.', 'markersize', 50)
% set(e([1 3]), 'color', AZblue)
% set(e([2 4]), 'color', AZred)
%
% saveFigurePdf(gcf, '~/Desktop/data')

%% FIGURE - RTs for model and subjects
figure(1); clf;
set(gcf, 'position', [811   575   600   500])
% ax = easy_gridOfEqualFigures([0.25 0.08], [0.15 0.15 0.03]);
ax = easy_gridOfEqualFigures([0.15 0.2 0.05], [0.12 0.15 0.05]);

e1 = plot_rtCurves_sub_v1(ax(1:2), sub, 0.1 , 3)
e2 = plot_rtCurves_sim_v1(ax(1:2), sim, 1)
set([e1([2 4]) e2([2 4])], 'visible', 'off')
% set(e1, 'linestyle', 'none')
% set(e2, 'marker', 'none')
% set(e2, 'linestyle', '--', 'color', 'k', 'marker', 'o', 'markersize', 10)

e1 = plot_rtCurves_sub_v1(ax(3:4), sub, 0.1, 3)
e2 = plot_rtCurves_sim_v1(ax(3:4), sim, 1)
set([e1([1 3]) e2([1 3])], 'visible', 'off')
% set(e2, 'linestyle', '--', 'color', 'k', 'marker', 'o', 'markersize', 10)
% set(e1, 'linestyle', 'none')
% set(e2, 'marker', 'none')
set(ax, 'ylim', [0.4 1.5])
axes(ax(1)); title('unequal [1 3], horizon 1')
axes(ax(2)); title('equal [2 2], horizon 1')
axes(ax(3)); title('unequal [1 3], horizon 6')
axes(ax(4)); title('equal [2 2], horizon 6')
saveFigurePdf(gcf, '~/Desktop/ddm')

%% choice curves
binEdges = [-25:10:25];
figure(2); clf;
set(gcf, 'position', [811   575   600   500])
ax = easy_gridOfEqualFigures([0.15 0.2 0.05], [0.12 0.15 0.05]);

e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3);; %plot_rtCurves_sub_v1(ax(1:2), sub)
e2 = plot_choiceCurves_sim_v1(ax(1:2), sim, binEdges);; %plot_rtCurves_sim_v1(ax(1:2), sim, 1)
set(e1, 'linestyle', 'none')

set(e2, 'marker', 'none')

set([e1([2 4]) e2([2 4])], 'visible', 'off')
% set(e2, 'linestyle', '--', 'color', 'k', 'marker', 'o', 'markersize', 10)

e1 = plot_choiceCurves_v2(ax(1:2), sub, binEdges, 0.1, 3); % plot_rtCurves_sub_v1(ax(3:4), sub)
e2 = plot_choiceCurves_sim_v1(ax(1:2), sim, binEdges); %plot_rtCurves_sim_v1(ax(3:4), sim, 1)
set([e1([1 3]) e2([1 3])], 'visible', 'off')
set(e1, 'linestyle', 'none')
set(e2, 'marker', 'none')
% set(e2, 'linestyle', '--', 'color', 'k', 'marker', 'o', 'markersize', 10)
% set(ax, 'ylim', [0.4 1.5])
axes(ax(1)); title('unequal [1 3], horizon 1')
axes(ax(2)); title('equal [2 2], horizon 1')
axes(ax(3)); title('unequal [1 3], horizon 6')
axes(ax(4)); title('equal [2 2], horizon 6')
saveFigurePdf(gcf, '~/Desktop/ddm_choices')
% e = plot_choiceCurves_v1(ax, sub, binEdges);
% e = plot_choiceCurves_sim_v1(ax, sim, binEdges);


%% choice curves
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
    
    
    ind = i13&i1;
    [M_13_1(:,sn), ~, X] = binIt(-dI(ind).*dR(ind), A(ind)==uID(ind), binEdges, 'std');
    
    ind = i13&i6;
    [M_13_6(:,sn), ~, X] = binIt(-dI(ind).*dR(ind), A(ind)==uID(ind), binEdges, 'std');
    
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

figure(1); clf;
set(gcf, 'position', [811   575   600   300])
ax = easy_gridOfEqualFigures([0.3 0.12], [0.12 0.12 0.05]);
axes(ax(1)); hold on;
e = errorbar(X, m_13_1, s_13_1);
e(2) = errorbar(X, m_13_6, s_13_6);
xlabel({'difference in mean reward' 'R(high info) - R(low info)'});
ylabel({'p(high info)'})
title('unequal [1 3]', 'fontweight', 'normal')
leg = legend(e([2 1]), {'horizon 6' 'horizon 1'}, 'location', 'northwest');
set(leg, 'position', [ 0.3083    0.3150    0.1683    0.1317])
axes(ax(2)); hold on;
e(3) = errorbar(X, m_22_1, s_22_1);
e(4) = errorbar(X, m_22_6, s_22_6);
xlabel({'difference in mean reward' 'R(left) - R(right)'});
ylabel({'p(left)'})
title('equal [2 2]', 'fontweight', 'normal')


set(ax, 'xlim', [-35 35], 'ylim', [0 1], 'tickdir', 'out', 'fontsize', 20)
set(e, 'linewidth', 3, 'marker', '.', 'markersize', 50)
set(e([1 3]), 'color', AZblue)
set(e([2 4]), 'color', AZred)

saveFigurePdf(gcf, '~/Desktop/choices')


%% horizon task description
%% distribution example
figure(1); clf;

hg = [0.1 0.1];
wg = [0.5 0.03];
[ax, hb, wb] = easy_gridOfEqualFigures(hg, wg);

axes(ax(1)); hold on;
x = [0:0.1:100];
y1 = exp(-(x-40).^2/2/8^2);
y2 = exp(-(x-60).^2/2/8^2);

plot(x, y1, 'linewidth', 5, 'color', AZred)
plot(x, y2, 'linewidth', 5, 'color', AZblue)
xlabel('reward')
ylabel('probability')
set(ax(1), 'yticklabel', [])


%% task example
figure(1); clf;
set(gcf, 'Position', [811   572   600   600])
hg = [0.1 0.1];
wg = [0.5 0.03];
[ax, hb, wb] = easy_gridOfEqualFigures(hg, wg);

axes(ax(1)); hold on;


R1 = {'XX' 'XX' '65' 'XX'}
R2 = {'61' '46' 'XX' '58'}
x0 = 0;
y0 = 0;
plot(x0+[0.5 0.5 1.5 1.5 0.5 0.5], y0+[0.5 10.5 10.5 0.5 0.5 1.5], '-','color', AZred)
for i = 1:9
    plot(x0+[0.5 1.5], y0+(i+0.5)*[1 1], '-','color', AZred)
end

x0 = 2;
y0 = 0;
plot(x0+[0.5 0.5 1.5 1.5 0.5 0.5], y0+[0.5 10.5 10.5 0.5 0.5 1.5], '-','color', AZblue)
for i = 1:9
    plot(x0+[0.5 1.5], y0+(i+0.5)*[1 1], '-','color', AZblue)
end

for i = 1:4
    text(1, i, R1{i}, 'horizontalalignment', 'center', 'fontsize', 20)
    text(3, i, R2{i}, 'horizontalalignment', 'center', 'fontsize', 20)
end

plot([0 4 4 0 0 0], [11 11 0 0 11 0],'k-')

xlim([0 4])
ylim([0 11])
set(ax(1), 'ydir', 'reverse', 'visible', 'off')


