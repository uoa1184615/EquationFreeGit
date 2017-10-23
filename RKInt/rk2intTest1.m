%{
% Script to test the time integration function rk2int() 
% on a simple nonlinear ODE.  
% AJR, 29 Mar 2017
%!TEX root = ../equationFreeDoc.tex
\subsection{A 1D test of RK2 integration}
\label{sec:2tpi}

Try the nonlinear, but separable, \ode
\begin{equation*}
\dot x=-2x/t
\quad\implies x=c/t^2.
\end{equation*}
\begin{matlab}
%}
fn=@(x,t) [-2*x/t]
%{
\end{matlab}
Solve \ode\ over \(1\leq t\leq 4\) with initial condition \(x(1)=1\)
\begin{matlab}
%}
ts=linspace(1,4,21);
[xs,errs]=rk2int(fn,1,ts)
plot(ts,xs,ts,100*errs,'o',ts,100*(xs-1 ./ts.^2),'x')
print -depsc2 RKInt/rk2intPlot1.eps
%{
\end{matlab}
\autoref{fig:1tpi} shows the step-errors would accumulate to a reasonably conservative estimate of the solution error.
\begin{figure}
\centering
\begin{tabular}{c@{}c}
\rotatebox{90}{\hspace{10ex}$x,\ 100\times$step-errs,\ $100\times$actual} &
\includegraphics[width=0.85\linewidth]{RKInt/rk2intPlot1}\\
& $t$
\end{tabular}
\caption{example trajectory and errors}
\label{fig:1tpi}
\end{figure}

%}
