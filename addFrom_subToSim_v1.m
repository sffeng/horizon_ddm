function sim = addFrom_subToSim_v1(sub, sim)

for sn = 1:length(sim)
    sim(sn).dR = sub(sn).o2(:,4) - sub(sn).o1(:,4);
    sim(sn).dI = -(sub(sn).n2(:,4) - sub(sn).n1(:,4)) / 2;
    sim(sn).gameLength = sub(sn).gameLength;
    sim(sn).h = sub(sn).gameLength-4;
end
