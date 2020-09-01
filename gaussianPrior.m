function logP = gaussianPrior(x, mu, sigma)

logP = -(x-mu).^2/sigma^2/2 - log(sqrt(2*pi)*sigma);