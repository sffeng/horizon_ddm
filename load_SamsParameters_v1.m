function sim = load_SamsParameters_v1(fname1, fname2)

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