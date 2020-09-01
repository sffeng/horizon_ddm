function XX = theory_qualitativeRT_type2_v1(A1, A6, sigma1, sigma6, rho, changeFlag)

% horizon 1 parameters
cMu_0_1 = sqrt(   1 / (2 *sqrt(2)*sigma1 * rho) );
cBeta_R_1 = rho * cMu_0_1;
cBeta_I_1 = A1 / (2 * sqrt(2) * sigma1 * cMu_0_1);

% horizon 6 parameters
switch changeFlag
    
    case 1
        % cMu_0 changes
        cMu_0_6 = 1 / (2 * sqrt(2) * cBeta_R_1 * sigma6);
        
        % cBeta_I changes
        cBeta_I_6 = A6 / (2 * sqrt(2) * sigma6 * cMu_0_6 );
        
        % rest stays the same
        cBeta_R_6 = cBeta_R_1;
        
    case 2
        
        % cBeta_R changes
        cBeta_R_6 = 1 / (2 * sqrt(2) * cMu_0_1 * sigma6);
        
        % cMu_I changes
        cBeta_I_6 = A6 / (2 * sqrt(2) * sigma6 * cMu_0_1);
        
        % rest stays the same
        cMu_0_6 = cMu_0_1;
        
end

% packup parameters into XX structure
XX = [
    cMu_0_1   cMu_0_6
    0         0
    0         0
    0         0
    cBeta_R_1 cBeta_R_6
    cBeta_I_1 cBeta_I_6
    0         0
    0         0
    0         0
    0.05      0.05];