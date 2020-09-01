% ddmpdf.m
%
% [dpdf,meanRT,phit] = ddmpdf(T,x,A,T0,x0,z,sigma)
%
% Computes the pdf for the pure ddm using Navarro and Fuss'
% numerical method
%
% INPUTS: 
% T: array of time points
% x: 0 or 1 to return pdf for lower/upper boundary
% A: Drift rate
% T0: Nondecision time (Just a shift)
% x0: Initial condition
% z: Decision thresholds are set at + and - z
% sigma: Standard deviation of noise increments
%
% Samuel Feng, 3 July 2013

function [dpdf] = ddmpdf(T,x,A,T0,x0,z,sigma)

% if abs(A)<0.00001
%     error('Drift rate is dangerously close to 0');
% end

s = sigma; errtol = 1e-4; nT = length(T);

dpdf = zeros(1,nT);
for k = 1:nT
    t = T(k)-T0;
    if t < 0
        dpdf(k) = 0; 
    else
        dpdf(k) = wfpt(t, x, A/s, 2*z/s, (x0+z)/s, errtol);
        % wfpt(t,x,drift rate,boundary,starting point,err)
    end
end
dpdf(isnan(dpdf))=0;


% flipped = 0;
% if A<0
%     x0 = -x0;
%     A = -A;
%     flipped = 1;
% end
% 
% zz = z/A;
% aa = (A/s)^2;
% xx = x0/A;
% 
% ER = 1 / (1 + exp(2*zz*aa)) - ...
%      (1 - exp(-2*xx*aa))/(exp(2*zz*aa)-exp(-2*zz*aa));
% 
% if x == 0 
%     phit = ER;
% else
%     phit = 1-ER;
% end
% 
% if flipped
%     phit = 1 - phit;
% end
% 
% if A == 0
%     phit = .5;
% end

% Compute expected value, careful to normalize
% dt = T(2)-T(1);
% I = sum(trapz(dpdf)*dt);
% cpdf = (1/I) .* dpdf;
% meanRT = sum(trapz(T.*cpdf)*dt);