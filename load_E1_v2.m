function sub = load_E1_v2(fname)

% fname = 'allHorizonData_v1.csv';
% also includes gID data for repeated games

%% read data from spreadsheet
fid = fopen(fname);

hdr = textscan(fid,...
    '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s', ...
    1, 'delimiter', ',');
data = textscan(fid,...
    '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', ...
    'delimiter', ',');
fclose(fid);


%% separate data into separate subjects
R = [data{12:21}];
C = [data{22:31}];
RT = [data{32:41}];

for sn = 1:max(data{3})
    
    ind = find(data{3} == sn);
    
    sub(sn).expt_name   = data{1}{ind(1)};
    sub(sn).subjectID   = data{2}{ind(1)};
    sub(sn).sNum        = data{3}(ind(1));
    
    sub(sn).block       = data{4}(ind);
    sub(sn).game        = data{5}(ind);
    sub(sn).gameLength  = data{6}(ind);
    sub(sn).uc          = data{7}(ind);
    sub(sn).m1          = data{8}(ind);
    sub(sn).m2          = data{9}(ind);
    sub(sn).gID         = data{10}(ind);
    sub(sn).repeatNumber = data{11}(ind);
    
    sub(sn).r           = R(ind,:);
    sub(sn).a           = C(ind,:);
    sub(sn).RT          = RT(ind,:);
    
end
clear RT C R

%% augment data structure
for sn = 1:length(sub)
    
    % z-scored RT
    sub(sn).RTz = (sub(sn).RT - nanmean(sub(sn).RT(:))) / nanstd(sub(sn).RT(:));
    
    % running total of how many times each bandit is played
    sub(sn).n1 = cumsum(sub(sn).a == 1,2);
    sub(sn).n2 = cumsum(sub(sn).a == 2,2);
    
    
    
    % running total of reward from each bandit
    sub(sn).R1 = cumsum(sub(sn).r.*(sub(sn).a==1),2);
    sub(sn).R2 = cumsum(sub(sn).r.*(sub(sn).a==2),2);
    
    % running observed mean for each bandit
    sub(sn).o1 = sub(sn).R1 ./ sub(sn).n1;
    sub(sn).o2 = sub(sn).R2 ./ sub(sn).n2;
    
    % is choice objectively correct?
    sub(sn).co = repmat((sub(sn).m1 > sub(sn).m2), [1 10]) .* (sub(sn).a==1) ...
        + repmat((sub(sn).m1 < sub(sn).m2), [1 10]) .* (sub(sn).a==2);
    sub(sn).co(isnan(sub(sn).a)) = nan;
    
    % is choice a low observed mean choice? (RANDOM EXPLORATION)
    sub(sn).lm = ...
        (sub(sn).o1(:,1:end-1) < sub(sn).o2(:,1:end-1)) .* (sub(sn).a(:,2:end)==1) + ...
        (sub(sn).o1(:,1:end-1) > sub(sn).o2(:,1:end-1)) .* (sub(sn).a(:,2:end)==2);
    sub(sn).lm(sub(sn).o1(:,1:end-1)==sub(sn).o2(:,1:end-1)) = nan;
    sub(sn).lm(isnan(sub(sn).a(:,2:end))) = nan;
    % shift over so that trials line up for later
    sub(sn).lm(:,2:end+1)=sub(sn).lm(:,1:end);
    sub(sn).lm(:,1) = nan;
    
    % is choice high info choice? (DIRECTED EXPLORATION)
    sub(sn).hi = (sub(sn).n1(:,1:end-1) < sub(sn).n2(:,1:end-1)) .* (sub(sn).a(:,2:end)==1) ...
        + (sub(sn).n1(:,1:end-1) > sub(sn).n2(:,1:end-1)) .* (sub(sn).a(:,2:end)==2);
    sub(sn).hi(sub(sn).n1(:,1:end-1)==sub(sn).n2(:,1:end-1)) = nan;
    sub(sn).hi(isnan(sub(sn).a(:,2:end))) = nan;
    % shift over so that trials line up for later
    sub(sn).hi(:,2:end+1)=sub(sn).hi(:,1:end);
    sub(sn).hi(:,1) = nan;
    
    % is choice same as the last thing they did? (ALTERNATION)
    sub(sn).rep = [nan(size(sub(sn).a, 1),1) (sub(sn).a(:,1:end-1) == sub(sn).a(:,2:end))];
    
    sub(sn).dR = sub(sn).o2(:,4) - sub(sn).o1(:,4);
    sub(sn).dI = -(sub(sn).n2(:,4) - sub(sn).n1(:,4))/2;
    sub(sn).choice = sub(sn).a(:,5);
    sub(sn).rt = sub(sn).RT(:,5);
        
    
    
end
