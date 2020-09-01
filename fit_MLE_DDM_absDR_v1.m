function fit = fit_MLE_DDM_absDR_v1(sub, M, RTmin, RTmax, X0_all, LB_all, UB_all)

% this one uses the dR, dI, choice and rt format

% initial parameters and bounds
X0 = X0_all(M==1);
LB = LB_all(M==1);
UB = UB_all(M==1);


clear Xfit
gl_vals = [5 10];
for sn = 1:length(sub)
    for i = 1:length(gl_vals)
        
        rt = sub(sn).rt;
        ch = sub(sn).choice==2;
        gl = sub(sn).gameLength;
        ind = (rt > RTmin) & (rt < RTmax) & (gl == gl_vals(i));
        ch = ch(ind); rt = rt(ind);
        dR = sub(sn).dR(ind);
        dI = sub(sn).dI(ind);
        
        obFunc = @(x) lik_DDMregX_absDR_v1(rt, ch, dR, dI, 1, M, x);
        
        [Xfit(:,sn,i), nLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB);
        
        n = length(dR);
        k = sum(M);
        BIC = k * log(n) + 2 * nLL;
        AIC = 2*k + 2 * nLL;
        
        fit(sn).gameLength = sub(sn).gameLength;
        fit(sn).dR = sub(sn).dR;
        fit(sn).dI = sub(sn).dI;
        fit(sn).RT = sub(sn).rt;
        fit(sn).choice = sub(sn).choice;
        
        fit(sn).M = M;
        fit(sn).Xfit(:,i) = Xfit(:,sn,i);
        fit(sn).XXfit(:,i) = zeros(length(M),1);
        fit(sn).XXfit(M==1, i) = Xfit(:,sn,i);
        
        fit(sn).n(i) = n;
        fit(sn).k(i) = k;
        fit(sn).BIC(i) = BIC;
        fit(sn).AIC(i) = AIC;
        
    end
end


