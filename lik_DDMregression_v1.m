function LL = lik_DDMregression_v1(RT, C, dR, dI, ...
    cA_0, cA_R, cA_I, ...
    cZ_0, cZ_R, cZ_I, ...
    cX_0, cX_R, cX_I, ...
    T0, c);

for i = 1:length(RT)
    A(i) = cA_0 + cA_R * dR(i) + cA_I * dI(i);
    z(i) = cZ_0 + cZ_R * dR(i) + cZ_I * dI(i);
    b(i) = cX_0 + cX_R * dR(i) + cX_I * dI(i);
    x0(i) = z(i)*(2*(1 ./ ( 1 + exp(-b(i)) )) - 1);
    
    
    [dpdf(i)] = ddmpdf(RT(i),C(i),A(i),T0,x0(i),z(i),c);
end

LL = nansum(log(dpdf));