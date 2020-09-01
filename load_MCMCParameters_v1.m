function fit = load_SamsParameters_v1(fname1, fname2)

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
    
    sim(sn).x1_0  = dum(16)*2;
    sim(sn).x1_dR = dum(17);
    sim(sn).x1_dI = -dum(18);
    sim(sn).x6_0  = dum(7)*2;
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
% all variables are present so map is:
%           cA_0 | cA_R | cA_I | cZ_0 | cZ_R | cZ_I | cX_0 | cX_R | cX_I |  T0
M      = [    1      1      1      1      1      1      1      1      1      1     ];

% bit of a hack, but this is putting it into a form that we can more easily
% compare with MLE fits later
for sn = 1:length(sim)
    
    %fit(sn).gameLength = sub(sn).gameLength;
    %fit(sn).dR = sub(sn).dR;
    %fit(sn).dI = sub(sn).dI;
    %fit(sn).RT = sub(sn).rt;
    %fit(sn).choice = sub(sn).choice;
    fit(sn).sID = sim(sn).sID;
    fit(sn).M = M;
    fit(sn).XXfit(1,1) = -sim(sn).A1_0;
    fit(sn).XXfit(1,2) = -sim(sn).A6_0;
    fit(sn).XXfit(2,1) = sim(sn).A1_dR;
    fit(sn).XXfit(2,2) = sim(sn).A6_dR;
    fit(sn).XXfit(3,1) = sim(sn).A1_dI;
    fit(sn).XXfit(3,2) = sim(sn).A6_dI;
    
    fit(sn).XXfit(4,1) = sim(sn).z1_0;
    fit(sn).XXfit(4,2) = sim(sn).z6_0;
    fit(sn).XXfit(5,1) = -sim(sn).z1_dR;
    fit(sn).XXfit(5,2) = -sim(sn).z6_dR;
    fit(sn).XXfit(6,1) = -sim(sn).z1_dI;
    fit(sn).XXfit(6,2) = -sim(sn).z6_dI;
    
    fit(sn).XXfit(7,1) = -sim(sn).x1_0;
    fit(sn).XXfit(7,2) = -sim(sn).x6_0;
    fit(sn).XXfit(8,1) = sim(sn).x1_dR;
    fit(sn).XXfit(8,2) = sim(sn).x6_dR;
    fit(sn).XXfit(9,1) = sim(sn).x1_dI;
    fit(sn).XXfit(9,2) = sim(sn).x6_dI;
    
    fit(sn).XXfit(10,1) = sim(sn).T01;
    fit(sn).XXfit(10,2) = sim(sn).T06;
    
    fit(sn).Xfit = fit(sn).XXfit; % redundant for full model
    
    
    fit(sn).n = nan;
    fit(sn).k = sum(M);
    fit(sn).BIC = nan;
    fit(sn).AIC = nan;
end
