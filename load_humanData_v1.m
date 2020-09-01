function sub = load_humanData_v1(datadir, dataname, demoname)

% load and augment data
sub = load_E1_v2([datadir dataname]);

% princeton subjects only
i1 = strcmp({sub.expt_name}, 'pilot-v1');% 'repeater-v1'})
i2 = strcmp({sub.expt_name}, 'repeater-v1');
sub = sub(i1 | i2 );

for sn = 1:length(sub)
    sub(sn).sID = sub(sn).sNum; % align some notations for later on
end

%% load demographics
fname = demoname;
fid = fopen(fname);
hdr = textscan(fid, '%s%s%s', 1, 'delimiter', ',');
data = textscan(fid, '%s%f%s', 'delimiter', ',');
fclose(fid)

for i = 1:length(data{1})
    nm{i} = data{1}{i};
    ge{i} = data{3}{i};
end
ag = data{2};

%% pair demographics with subject
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