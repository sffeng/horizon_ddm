function fit_MLE = theory_ERDT_abs_v1(fit_MLE, r_vals)

% r_vals = [-30:0.01:30];
for sn = 1:length(fit_MLE)
    for i = 1:2
        
        [ER, RT] = compute_contERRT_abs_v2(r_vals, +1, fit_MLE(sn).XXfit(:,i));
        fit_MLE(sn).xx = r_vals;
        fit_MLE(sn).cc13_ana(:,i) = 1-ER;
        fit_MLE(sn).rt13_ana(:,i) = RT;
        
        [ER, RT] = compute_contERRT_abs_v2(r_vals, 0, fit_MLE(sn).XXfit(:,i));
        fit_MLE(sn).cc22_ana(:,i) = 1-ER;
        fit_MLE(sn).rt22_ana(:,i) = RT;
        
    end
end