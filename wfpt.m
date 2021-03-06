% Modified by Sam Feng
%
% original code from Navarro Fuss 2009

function p=wfpt(t,x,v,a,z,err)

% calculates probability that particle reaches lower bound (or upper) at
% time t

% t = time of first passage
% x = 1 upper boundary; x=0 lower boundary
% v = drift rate
% a = "boundary separation" (threshold) how much more evidence for A than for B?
% z = starting point
% err = target approximation error level

tt=t/(a^2); % use normalized time
w=z/a; % convert to relative start point

if x == 1  % if you want upper boundary, then do this
    v = -v;
    w = 1-w;
end

% calculate (minimum) number of terms needed for large t
if pi*tt*err<1 % if error threshold is set low enough
    kl=sqrt(-2*log(pi*tt*err)./(pi^2*tt)); % bound % Equation 10
    kl=max(kl,1/(pi*sqrt(tt))); % ensure boundary conditions met % line 160
else % if error threshold set too high
    kl=1/(pi*sqrt(tt)); % set to boundary condition % line 160
end

% calculate (minimum) number of terms needed for small t
if 2*sqrt(2*pi*tt)*err<1 % if error threshold is set low enough
    ks=2+sqrt(-2*tt.*log(2*sqrt(2*pi*tt)*err)); % bound % Equation 11
    ks=max(ks,sqrt(tt)+1); % ensure boundary conditions are met
else % if error threshold was set too high
    ks=2; % minimal kappa for that case
end

% compute f(tt|0,1,w)
p=0; %initialize density
if ks<kl % if small t is better (i.e., lambda<0)
	K=ceil(ks); % round to smallest integer meeting error
    for k=-floor((K-1)/2):ceil((K-1)/2) % loop over k
        p=p+(w+2*k)*exp(-((w+2*k)^2)/2/tt); % increment sum
    end
    p=p/sqrt(2*pi*tt^3); % add constant term
        
else % if large t is better... (i.e. lambda>=0)
	K=ceil(kl); % round to smallest integer meeting error
	for k=1:K
		p=p+k*exp(-(k^2)*(pi^2)*tt/2)*sin(k*pi*w); % increment sum
	end
	p=p*pi; % add constant term
end

% convert to f(t|v,a,w)
p=p*exp(-v*a*w -(v^2)*t/2)/(a^2); 

