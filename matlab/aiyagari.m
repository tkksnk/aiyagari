% Bewley-Aiyagari model
% May 2020, Takeki Sunakawa
clear all;

BETA  = 0.96;
ALPHA = 0.36; 
DELTA = 1-0.92;
SIGMA = 1.5;
mu = 0.05;
pee = 0.925;
puu = 0.5;

% grid points
ne = 2;
Ge = [1.0; mu]; 
Pe = [pee 1-pee;
    1-puu puu]; 

nb = 1001;
kmin = 0;
kmax = 20.0;
knotsb = linspace(kmin,kmax,nb)';
% knotsb = logspace(log(kbounds(1) - 1.0d0*kbounds(1) + 1.0d0)/log(10.0d0), log(kbounds(2) - 1.0d0*kbounds(1) + 1.0d0)/log(10.0d0), nb)';
% knotsb = knotsb + kbounds(1) - 1.0d0;

mue = Pe^10000;
mue = mue(1,:)';
znow = 1.0d0;
% efficiency unit of labor
lnow = Ge'*mue;

% initial distribution
vmat0 = zeros(nb,ne);
mu0 = ones(nb,ne)/(nb*ne);

% initial value of r
mnow = lnow*((1.0-BETA*(1.0-DELTA))/(ALPHA*BETA))^(1.0/(ALPHA-1.0));
mnow = 5.2074; % from My_Aiyagari.m
m0 = mnow;
r0 = (ALPHA)*znow*mnow^(ALPHA-1)*lnow^(1-ALPHA);


start = tic;
critin = 1e-4;
critmu = 1e-8;
critout = 1e-3;
diffout = 10.0;
damp = 0.01;
iter = 0;
maxiter = 1000;

while (diffout>critout && iter<maxiter)

    [r1,vmat0,mu0] = driver(r0,vmat0,mu0,Ge,Pe,knotsb,BETA,ALPHA,DELTA,SIGMA,znow,lnow,critin,critmu);

    diffout = abs(log(r1)-log(r0));    
    iter = iter+1;    
    disp(sprintf("iter = %4d, diff = %5.6f, oldr = %5.6f, newr = %5.6f",iter,diffout,r0,r1))

    % Update K
    r0 = damp*r1 + (1.0-damp)*r0;
    
end

t = toc(start);
disp(sprintf("Elapsed time is %5.10f.",t));