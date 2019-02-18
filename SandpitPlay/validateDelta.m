%{
Validates dD>log(bD)/bD rule for PIRK2 by finding the dD at
which the PIRK2 scheme reduces by a factor of 10 each step,
given the fast decay is a rate bD.  Reducing by a factor of
two is barely any different.  AJR, 5 Feb 2019
%}
clear all, close all
global betaD
betaDs=logspace(log10(exp(1)),3,41);
bTtenth=nan(size(betaDs));
for j=1:length(betaDs)
    betaD=betaDs(j);
    bT0=log(betaD)/betaD;
    bTtenth(j)=fsolve(@pirk2fac,bT0);
end
clf()
loglog(betaDs,[log(betaDs)./betaDs; bTtenth]')
grid
xlabel('\beta\Delta'),ylabel('(\delta/\Delta)_{min}')

function raterr=pirk2fac(bT)
stepFactor=0.1;
x = PIRK2(@eburst, bT, 0:2, 1);
raterr=x(end,:)./x(end-1,:) - stepFactor;
end

function [ts, xs] = eburst(ti, xi, bT) 
global betaD
    ts = linspace(ti,ti+bT,11)';
    xs = xi.*exp(-betaD.*(ts-ti));
end
