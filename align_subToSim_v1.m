function [sub, sub_all] = align_subToSim_v1(sub, sim)

% also keep sub_all for later reporting of demographics (so we know who we
% excluded)

sID = [sim.sID];
[a,idx,c] = intersect([sub.sID], sID);
sub_all = sub;
for sn = 1:length(sub);%length(idx)
    if sum(idx==sn) > 0
        sub_all(sn).excluded = 0;
    else
        sub_all(sn).excluded = 1;
    end
end

sub = sub(idx);
