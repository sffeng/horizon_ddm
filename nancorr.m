function [r,p] = nancorr(X, Y, type)


ind1 = sum(isnan(X),2) == 0;
ind2 = sum(isnan(Y),2) == 0;
ind = ind1 & ind2;
X = X(ind,:);
Y = Y(ind,:);
[r,p] = corr(X, Y, 'type', type);
