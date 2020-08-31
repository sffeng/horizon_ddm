function [binMeans, binContents, binCentres, binSem] = ...
    binIt(X, Y, binEdges, type)

% [binMeans, binContents, binCentres,  binSem] = ...
%    binIt(X, Y, binEdges, type)
% Bins Y, as a function of X, with bin edges given by binEdges


% Robert Wilson
% 08/04/11

if exist('type') ~= 1
    type = 'std';
end

binSpace = diff(binEdges);

for i = 1:length(binEdges)+1
    
    
    switch i
        
        case 1
            
            ind = X <= binEdges(i);
            binCentres(i) = binEdges(i) - max(binSpace)/2;
            
        case length(binEdges)+1
            
            ind = X > binEdges(end);
            binCentres(i) = binEdges(end) + max(binSpace)/2;
            
            
        otherwise
            
            ind = (X > binEdges(i-1)) & (X <= binEdges(i) );
            binCentres(i) = ( binEdges(i-1) + binEdges(i) ) / 2;
            
    end
    
    binContents{i} = Y(ind);
    binMeans(i) = nanmean(binContents{i});
    
    switch type
    
        case 'beta'
            
            alpha = sum(binContents{i}==1);
            beta = sum(binContents{i}==0);
            binSem(i) = ...
                sqrt(alpha * beta / (alpha + beta)^2 / (alpha + beta + 1));
    
        case 'std'
    
            binSem(i) = nanstd(binContents{i}) ...
                / sqrt(sum(~isnan(binContents{i})));
        
    end
end
