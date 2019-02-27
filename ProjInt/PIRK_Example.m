% Linear example of PIRK4(). JM and AJR, Sept 18 -- Feb 2019.
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{Example: PI using Runge--Kutta macrosolvers }
\label{sec:ExPIRK}
\localtableofcontents

This script is a demonstration of the \verb|PIRK()| schemes,
that use a Runge--Kutta macrosolver, applied to a simple
linear system with some slow and fast directions. 

\begin{devMan}
Clear workspace and set a seed.
\begin{matlab}
%}
clear
rng(1)
global dxdt
%{
\end{matlab}

The majority of this example involves setting up details for
the microsolver. We use a simple function
\verb|gen_linear_system()| that outputs a function \({f(t,x)
= {A}\vec x+\vec b}\), where matrix~\({A}\) has some
eigenvalues with large negative real part, corresponding to
fast variables and some eigenvalues with real part close to
zero, corresponding to slow variables. The function
\verb|gen_linear_system()| requires that we specify bounds
on the real part of the strongly stable eigenvalues,
\begin{matlab}
%}
fastband = [-5e2; -1e2]; 
%{
\end{matlab}
and bounds on the real part of the weakly stable/unstable
eigenvalues,
\begin{matlab}
%}
slowband = [-0.002; 0.002];
%{
\end{matlab}
We now generate a random linear system with seven fast and
three slow variables.
\begin{matlab}
%}
dxdt = gen_linear_system(7,3,fastband,slowband); 
%{
\end{matlab}



Set the macroscale times at which we request output from the
\textsc{pi} scheme and the initial conditions.
\begin{matlab}
%}
tSpan = 0: 1 : 20; 
x0 = linspace(-10,10,10)'; 
%{
\end{matlab}



We implement the \textsc{pi} scheme, saving the coarse
states in \verb|x|, the `trusted' applications of the
microsolver in \verb|tms| and~\verb|xms|, and the additional
applications of the microsolver in~\verb|rm| (the second,
third and fourth outputs are optional).
\begin{matlab}
%}
[x, tms, xms, rm] = PIRK4(@linearBurst, tSpan, x0); 
%{
\end{matlab}
To verify, we also compute the trajectories using a
standard integrator.
\begin{matlab}
%}
[tt,ode45x] = ode45(dxdt,tSpan([1,end]),x0);
%{
\end{matlab}

\begin{figure}
\caption{Demonstration of PIRK4(). From initial conditions,
the system rapidly trannsitions to an attracting invariant
manifold. The \textsc{pi} solution accurately tracks the
evolution of the variables over time while requiring only a
fraction of the computations of the standard
solver.}\label{fig:PIRK}
\includegraphics[scale=0.9]{PIRK}
\end{figure}
\cref{fig:PIRK} plots the output.
\begin{matlab}
%}
clf()
hold on
PI_sol=plot(tSpan,x,'bo');
std_sol=plot(tt,ode45x,'r');
plot(tms,xms,'k.', rm.t,rm.x,'g.');
legend([PI_sol(1),std_sol(1)],'PI Solution',...
    'Standard Solution','Location','NorthWest')
xlabel('Time'), ylabel('State')
%{
\end{matlab}
Save plot to a file.
\begin{matlab}
%}
set(gcf,'PaperPosition',[0 0 14 10]), print('-depsc2','PIRK')
%{
\end{matlab}

Code the micro-burst function using simple Euler steps. As a
rule of thumb, the time-steps \verb|dt| should satisfy
$\verb|dt|  \le 1/|\verb|fastband|(1)|$ and the time to
simulate with each application of the microsolver,
\verb|bT|, should be larger than or equal to
$1/|\verb|fastband|(2)|$. We set the integration scheme to
be used in the microsolver. Since the time-steps are so
small, we just use the forward Euler scheme
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
\end{devMan}
%}
