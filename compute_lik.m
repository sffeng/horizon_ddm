function LL = compute_lik(RT, C, A, T0, x0, z, c)

for i = 1:length(RT)
    [dpdf(i)] = ddmpdf(RT(i),C(i),A,T0,x0,z,c);
end

LL = nansum(log(dpdf));