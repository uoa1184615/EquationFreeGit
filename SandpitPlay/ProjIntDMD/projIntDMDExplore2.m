% Explore the time integration function projIntDMD() 
% on discretisation of a nonlinear diffusion PDE.  
% Test dependence of errors on macro-time-step.
% AJR, Feb 2018
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{projIntDMDExplore2}: explore effect of varying parameters}
\label{sec:pi1eevp}
\localtableofcontents

Seek to simulate the nonlinear diffusion \pde\
\begin{equation*}
\D tu=u\DD xu\quad \text{such that }u(\pm 1)=0,
\end{equation*}
with random positive initial condition.

Set the number of interior points in the domain~\([-1,1]\), and the macroscale time-step. 
\begin{matlab}
%}
function projIntDMDExplore2
n=9
dt=2/n^2
ICNoise=0
%{
\end{matlab}

Set micro-simulation parameters.
Rank two is fine when starting on the slow manifold.
Choose middle of the road transient and analysed time.
\begin{matlab}
%}
rank=2
timeSteps=[0.2 0.2]
%{
\end{matlab}
Try integrating with macro time-steps up to this sort of magnitude.
\begin{matlab}
%}
Ttot=9
%{
\end{matlab}





Set the initial condition to parabola or some skewed random positive values.
Without noise this initial condition is already on the slow manifold so only little reason for transient time.
\begin{matlab}
%}
x=linspace(-1,1,n+2)';
u0=(0.5+ICNoise*rand(n+2,1)).*(1-x.^2);
%{
\end{matlab}

First find a reference solution of the microscale dynamics over all time, here stored in \verb|Uss|.
\begin{matlab}
%}
[Us,Uss,Tss]=projIntDMD(@dudt,u0,[0 Ttot],2,dt,[0 Ttot]);
%{
\end{matlab}


Projectively integrate two steps in time with various parameters.
But remember that \verb|projIntDMD| rounds \verb|timeSteps| etc to nearest multiple of~\verb|dt|, so some of the following is a little dodgy but should not matter for overall trend.
\begin{matlab}
%}
Dts=0.1*[1 2 4 6 10 16 26]
errs=[]; relerrs=[]; DTs=[];
for p=Dts
[~,j]=min(abs(sum(timeSteps)+p-Tss))
ts=Tss(j)*(0:2)
js=1+(j-1)*(0:2);
[us,uss,tss]=projIntDMD(@dudt,u0,ts,rank,dt,timeSteps);
%{
\end{matlab}
Plot the macroscale predictions 
\begin{matlab}
%}
if 1
  clf,plot(x,Uss(:,js),'o-',x,us,'x--')
  xlabel('space x'),ylabel('u(x,t)')
  pause(0.01)
end
%{
\end{matlab}
Accumulate errors as function of time.
\begin{matlab}
%}
err=sqrt(sum((us-Uss(:,js)).^2))
errs=[errs;err];
relerrs=[relerrs;err./sqrt(sum(Uss(:,js).^2))];
%{
\end{matlab}

End the loop over parameters.
\begin{matlab}
%}
end
%{
\end{matlab}

Plot errors
\begin{matlab}
%}
loglog(Dts,errs(:,2:3),'o:')
xlabel('projective time-step')
ylabel('steps error')
legend('one','two')
grid
matlab2tikz('pi1x2.ltx')
%{
\end{matlab}


End the main function (not needed for new enough Matlab).
\begin{matlab}
%}
end
%{
\end{matlab}





\paragraph{The nonlinear PDE discretisation}
Code the simple centred difference discretisation of the nonlinear diffusion \pde\ with constant (usually zero) boundary values.
\begin{matlab}
%}
function ut=dudt(t,u)
n=length(u);
dx=2/(n-1);
j=2:n-1;
ut=[0
    u(j).*(u(j+1)-2*u(j)+u(j-1))/dx^2
    0];
end
%{
\end{matlab}
%}
