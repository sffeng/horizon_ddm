function [sub, sub_bad] = removeBadSubjects_E1_v1(sub, thresh)

% [sub, sub_bad] = removeBadSubjects_E1_v1(sub)

for sn = 1:length(sub)
    dum = sub(sn).co(:,end);
    sub(sn).fC = nanmean(dum(:));
end

ind_good = [sub.fC]>=thresh;
sub_bad = sub(~ind_good);
sub = sub(ind_good);