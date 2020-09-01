function fit = fit_MLE_DDM_v1(sub, M, RTmin, RTmax, X0_all, LB_all, UB_all)

% initial parameters and bounds
X0 = X0_all(M==1);
LB = LB_all(M==1);
UB = UB_all(M==1);


clear Xfit
gl_vals = [5 10];
for sn = 1:length(sub)
    for i = 1:length(gl_vals)
        
        rt = sub(sn).RT(:,5);
        ch = sub(sn).a(:,5) == 2;
        gl = sub(sn).gameLength;
        ind = (rt > RTmin) & (rt < RTmax) & (gl == gl_vals(i));
        ch = ch(ind); rt = rt(ind);
        dR = sub(sn).o2(ind,4) - sub(sn).o1(ind,4);
        dI = (sub(sn).n2(ind,4) - sub(sn).n1(ind,4))/2;
        
        obFunc = @(x) lik_DDMregX_v1(rt, ch, dR, dI, 1, M, x);
        
        [Xfit(:,sn,i), nLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB);
        
        n = length(dR);
        k = sum(M);
        BIC = k * log(n) + 2 * nLL;
        AIC = 2*k + 2 * nLL;
        
        fit(sn).gameLength = sub(sn).gameLength;
        fit(sn).dR = sub(sn).o2(:,4) - sub(sn).o1(:,4);
        fit(sn).dI = -(sub(sn).n2(:,4) - sub(sn).n1(:,4))/2;
        fit(sn).RT = sub(sn).RT(:,5);
        fit(sn).choice = sub(sn).a(:,5);
        
        fit(sn).M = M;
        fit(sn).Xfit(:,i) = Xfit(:,sn,i);
        fit(sn).XXfit(:,i) = zeros(length(M),1);
        fit(sn).XXfit(M==1, i) = Xfit(:,sn,i);
        
        fit(sn).n = n;
        fit(sn).k = k;
        fit(sn).BIC = BIC;
        fit(sn).AIC = AIC;
        
    end
end


