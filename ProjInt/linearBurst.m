% Used by PIRKexample.m
% AJR, Apr 2019
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\paragraph{A micro-burst simulation}
Used by \verb|PIRKexample.m|. Code the micro-burst function
using simple Euler steps. As a rule of thumb, the time-steps
\verb|dt| should satisfy $\verb|dt|  \le
1/|\verb|fastband|(1)|$ and the time to simulate with each
application of the microsolver, \verb|bT|, should be larger
than or equal to $1/|\verb|fastband|(2)|$. We set the
integration scheme to be used in the microsolver. Since the
time-steps are so small, we just use the forward Euler
scheme
\begin{matlab}
%}
function [ts, xs] = linearBurst(ti, xi, varargin) 
global dxdt
dt = 0.001;
ts = ti+(0:dt:0.05)'; 
nts = length(ts);
xs = NaN(nts,length(xi));
xs(1,:)=xi;
for k=2:nts
    xi = xi + dt*dxdt(ts(k),xi.').';
    xs(k,:)=xi;
end
end
%{
\end{matlab}
%}
