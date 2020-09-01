function [nBestFit, fBestFit] = compute_nBestFit(XXX)

[~,ind] = min(XXX');

nBestFit = hist(ind,1:size(XXX,2));
fBestFit = nBestFit/ sum(nBestFit);