function [ER, RT] = compute_contERRT_abs_v2(r_vals, dI, XX)

% unpack XX
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

A    = cA_0 + cA_R * r_vals + cA_I * dI;
z    = cZ_0 + cZ_R * abs(r_vals) + cZ_I * abs(dI);
b    = cX_0 + cX_R * r_vals + cX_I * dI;

x0   = (2 ./ (1 + exp(-b))-1) .* z;

a = A.^2;
z = z./A;
x0 = x0./A;

RT = T0 + z .* tanh(z .* a) + ( 2*z.*(1-exp(-2*x0.*a))./( exp(2*z.*a) - exp(-2*z.*a) ) - x0  );
ER = 1 ./ ( 1 + exp(2*z.*a) ) - ( (1-exp(-2*x0.*a))./(exp(2*z.*a) - exp(-2*z.*a) ) );
