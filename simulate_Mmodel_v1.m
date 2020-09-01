function fak = simulate_Mmodel_v1(XX, dt, dR, dI, gameLength)

% unpack XX
cA_0 = XX(1,:);
cA_R = XX(2,:);
cA_I = XX(3,:);

cZ_0 = XX(4,:);
cZ_R = XX(5,:);
cZ_I = XX(6,:);

cX_0 = XX(7,:);
cX_R = XX(8,:);
cX_I = XX(9,:);

T0_0   = XX(10,:);

gl_vals = [5 10];

% get a drift rate, threshold and starting point for each trial
for i = 1:length(dR)
    j = find(gameLength(i) == gl_vals);
    
    A(i) = cA_0(j) + cA_R(j) * dR(i) + cA_I(j) * dI(i);
    z(i) = cZ_0(j) + cZ_R(j) * dR(i) + cZ_I(j) * dI(i);
    b(i) = cX_0(j) + cX_R(j) * dR(i) + cX_I(j) * dI(i);
    x0(i) = z(i)*(2*(1 ./ ( 1 + exp(-b(i)) )) - 1);
    T0(i) = T0_0(j);
end

% simulate each trial
for i = 1:length(dR)
    [~, ~, DT(i), C(i)] = simluate_DDM_v1(dt, A(i), 1, z(i), x0(i));
    RT(i) = DT(i) + T0(i);
end

fak.sID(1:length(dR),1) = nan;
fak.XX = XX;
fak.dR = dR;
fak.dI = dI;
fak.gameLength = gameLength;
fak.choice = C' + 1;
fak.rt = RT';
