% Michaelis--Menton example of projective integrating
% fast-slow system.  This example simply introduces basic
% usage of the PIRK2() function. AJR, 29 Sep 2018
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{egPIerrs}: Errors in projective integration
of Michaelis--Menton kinetics}
\label{sec:egPIerrs}
\localtableofcontents

The Michaelis--Menten enzyme kinetics is expressed as a
singularly perturbed system of differential equations for
\(x(t)\) and~\(y(t)\):
\begin{equation*}
\frac{dx}{dt}=-x+(x+\tfrac12)y \quad\text{and}\quad
\frac{dy}{dt}=\frac1\epsilon\big[x-(x+1)y\big].
\end{equation*}
The slow
variable~\(x(t)\) evolves on a time scale of one, whereas
the fast variable~\(y(t)\) evolves on a time scale of the
small parameter~\(\epsilon\).


\subsection{Invoke projective integration}

Clear, and set the scale separation parameter~\(\epsilon\)
to something small.  Need something like \(\epsilon=0.001\)
to get a clear error plot.
\begin{matlab}
%}
clear all
global epsilon odeopts
epsilon = 0.001
Tend=4   
x0=[1;0]
%{
\end{matlab}
Generate a reference end-value, assumed correct to this specified error.
\begin{matlab}
%}
odeopts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[ts, xs] = MMburst(0, x0, Tend);
xend=xs(end,:)
%{
\end{matlab}
Over a range of parameters: the projective step size, and the burst time length.
\begin{matlab}
%}
nlogdt=13
nburst=4
nSteps=ceil( 1.5 .^(1:nlogdt) )
meps=(4:nburst+3)
%{
\end{matlab}
Create storage for results.
\begin{matlab}
%}
xerrs=nan(nburst,nlogdt);
dts=nan(1,nlogdt);
%{
\end{matlab}

First, \cref{sec:egPIMMacc} encodes the computation of
bursts of the Michaelis--Menten system in a
function~\verb|MMburst()|. Second, here set macroscale times
of computation and interest into vector~\verb|ts|. Then,
invoke Projective Integration with \verb|PIRK2()| applied to
the burst function, say using bursts of simulations of
length~\(2\epsilon\), and starting from the initial
condition for the Michaelis--Menten system of
\((x(0),y(0))=(1,0)\) (off the slow manifold).
\begin{matlab}
%}
for j=1:nlogdt
  ts = linspace(0,Tend,nSteps(j)+1);
  dts(j)=diff(ts(1:2));
  for i=1:nburst
    xs = PIRK2(@MMburst, ts, x0, meps(i)*epsilon);
    xerrs(i,j)=norm(xs(end,:)-xend);
end, end
dts=dts
%{
\end{matlab}
Plot results.
\begin{matlab}
%}
clf(),loglog(dts,xerrs,'o:'), ylim([1e-6 1])
xlabel('macro time-step'), ylabel('error at Tend')
title(['\epsilon = ' num2str(epsilon) ',  for different burst lengths'])
legend(num2str(meps', '%i\\epsilon'))
grid off, grid
set(gcf,'PaperPosition',[0 0 14 10])
print -depsc2 egPIerrs.eps
%{
\end{matlab}


\subsection{Code an accurate burst of Michaelis--Menten enzyme kinetics}
\label{sec:egPIMMacc}
Say use \verb|ode45()| to accurately integrate a burst of the
differential equations for the Michaelis--Menten enzyme
kinetics. 
\begin{matlab}
%}
function [ts, xs] = MMburst(ti, xi, bT) 
    global epsilon odeopts
    dMMdt = @(t,x) [ -x(1)+(x(1)+0.5)*x(2)
          1/epsilon*( x(1)-(x(1)+1)*x(2) ) ];
    [ts, xs] = ode45(dMMdt, [ti ti+bT], xi, odeopts);
end
%{
\end{matlab}
%}
