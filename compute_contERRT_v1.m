function [ER, RT] = compute_contERRT_v1(r_vals, z1_0, z1_dR, z1_dI, x1_0, x1_dR, x1_dI, A1_0, A1_dR, A1_dI, T0)

z    = z1_0 + z1_dR * r_vals + z1_dI;
bias = x1_0 + x1_dR * r_vals + x1_dI;
A    = A1_0 + A1_dR * r_vals + A1_dI;

x0   = 1 ./ (1 + exp(-bias)) .* z;

a = A.^2;
z = z./A;
x0 = x0./A;

RT = T0 + z .* tanh(z .* a) + ( 2*z.*(1-exp(-2*x0.*a))./( exp(2*z.*a) - exp(-2*z.*a) ) - x0  );
ER = 1 ./ ( 1 + exp(2*z.*a) ) - ( (1-exp(-2*x0.*a))./(exp(2*z.*a) - exp(-2*z.*a) ) );
