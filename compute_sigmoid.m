function P = compute_sigmoid(xx, bonus, bias, noise)

dQ = xx + bonus;
    
if bonus == 0
    P = 1 ./ ( 1 + exp( -1/noise/sqrt(2) * dQ ));
else
    P = 1 ./ ( 1 + exp( -1/sqrt(noise^2+bias^2)/sqrt(2) * dQ ));
end