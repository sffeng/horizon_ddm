function [y, t, RT, C] = simluate_DDM_v1(dt, A, c, z, y0)

y(1) = y0;
flag = 1;
i = 1;
while flag
    dW = sqrt(dt) * randn(1);
    dy = A * dt + c * dW;
    y(i+1) = y(i) + dy;
    if abs(y(i+1)) > z
        flag = 0;
    end
    i = i + 1;
end

t = [0:length(y)-1]*dt;

RT = t(end);
C = (1+sign(y(end)))/2; % cross positive bound or not?