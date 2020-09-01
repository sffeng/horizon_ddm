function sim = compute_meanERRT_v1(sim)

for sn = 1:length(sim)
    for i = 1:length(sim(sn).dR)
        if sim(sn).h(i) == 1
            
            sim(sn).z(i,1)    = sim(sn).z1_0 + sim(sn).z1_dR * sim(sn).dR(i) + sim(sn).z1_dI * sim(sn).dI(i);
            sim(sn).bias(i,1) = sim(sn).x1_0 + sim(sn).x1_dR * sim(sn).dR(i) + sim(sn).x1_dI * sim(sn).dI(i);
            sim(sn).A(i,1)    = sim(sn).A1_0 + sim(sn).A1_dR * sim(sn).dR(i) + sim(sn).A1_dI * sim(sn).dI(i);
            sim(sn).T0(i,1)   = sim(sn).T01;
            
        else
            
            sim(sn).z(i,1)    = sim(sn).z6_0 + sim(sn).z6_dR * sim(sn).dR(i) + sim(sn).z6_dI * sim(sn).dI(i);
            sim(sn).bias(i,1) = sim(sn).x6_0 + sim(sn).x6_dR * sim(sn).dR(i) + sim(sn).x6_dI * sim(sn).dI(i);
            sim(sn).A(i,1)    = sim(sn).A6_0 + sim(sn).A6_dR * sim(sn).dR(i) + sim(sn).A6_dI * sim(sn).dI(i);
            sim(sn).T0(i,1)   = sim(sn).T06;
            
        end
        
        
        sim(sn).x0   = 2 ./ (1 + exp(-sim(sn).bias)) - 1;
        
    end
end

%% compute mean RT from DDM parameters
for sn = 1:length(sim)
    z = sim(sn).z;
    x0 = sim(sn).x0.*z;
    A = sim(sn).A;
    T0 = sim(sn).T0;
    
    a = A.^2;
    z = z./A;
    x0 = x0./A;
    
    RT = T0 + z .* tanh(z .* a) + ( 2*z.*(1-exp(-2*x0.*a))./( exp(2*z.*a) - exp(-2*z.*a) ) - x0  );
    ER = 1 ./ ( 1 + exp(2*z.*a) ) - ( (1-exp(-2*x0.*a))./(exp(2*z.*a) - exp(-2*z.*a) ) );
    sim(sn).RTana = RT;
    sim(sn).ER = ER;
end