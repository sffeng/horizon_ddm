function nLL = lik_DDMregX_absDR_v1(RT, C, dR, dI, c, M, X);

% convert input parameters into long parameter vector
XX = zeros(1, length(M));
XX(M==1) = X;

cA_0 = XX(1);
cA_R = XX(2);
cA_I = XX(3);

cZ_0 = XX(4);
cZ_R = XX(5);
cZ_I = XX(6);

cX_0 = XX(7);
cX_R = XX(8);
cX_I = XX(9);

T0   = XX(10);


for i = 1:length(RT)
    A(i) = cA_0 + cA_R * dR(i)      + cA_I * dI(i);
    z(i) = cZ_0 + cZ_R * abs(dR(i)) + cZ_I * abs(dI(i));
    b(i) = cX_0 + cX_R * dR(i)      + cX_I * dI(i);
    x0(i) = z(i)*(2*(1 ./ ( 1 + exp(-b(i)) )) - 1);
    
    
    [dpdf(i)] = ddmpdf(RT(i), C(i), A(i), T0, x0(i), z(i), c);
end

nLL = -nansum(log(dpdf));