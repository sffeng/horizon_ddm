function XX = theory_qualitativeRT_type1_v1(A1, A6, sigma1, sigma6, rho, changeFlag)

% horizon 1 parameters
cMu_R_1 = sqrt(   1 / (2 *sqrt(2)*sigma1 * rho) );
cBeta_0_1 = rho * cMu_R_1;
cMu_I_1 = A1 / (2 * sqrt(2) * sigma1 * cBeta_0_1);

% horizon 6 parameters
switch changeFlag
    
    case 1
        % cMu_R changes
        cMu_R_6 = 1 / (2 * sqrt(2) * cBeta_0_1 * sigma6);
        
        % cMu_I changes
        cMu_I_6 = A6 / (2 * sqrt(2) * sigma6 * cBeta_0_1 );
        
        % rest stays the same
        cBeta_0_6 = cBeta_0_1;
        
    case 2
        
        % cBeta_0 changes
        cBeta_0_6 = 1 / (2 * sqrt(2) * cMu_R_1 * sigma6);
        
        % cMu_I changes
        cMu_I_6 = A6 / (2 * sqrt(2) * sigma6 * cBeta_0_6);
        
        % rest stays the same
        cMu_R_6 = cMu_R_1;
        
end

% packup parameters into XX structure
XX = [
    0         0
    cMu_R_1   cMu_R_6
    cMu_I_1   cMu_I_6
    cBeta_0_1 cBeta_0_6
    0         0
    0         0
    0         0
    0         0
    0         0
    0.05      0.05];