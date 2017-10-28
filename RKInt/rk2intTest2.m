% Script to test the time integration function rk2int() 
% on a simple nonlinear system of ODEs.  
% AJR, 29 Mar 2017
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{\texttt{rk2intTest2}: A 2D test of RK2 integration}
\label{sec:2tpi}

Try the nonlinear fast-slow system
\begin{equation*}
\dot x=-xy\,,\quad \dot y=-y+x^2.
\end{equation*}
\begin{matlab}
%}
fn=@(x,t) [-x(1)*x(2);-x(2)+x(1)^2]
%{
\end{matlab}
Solve over \(0\leq t\leq5\) with initial condition way off the slow manifold of \(x(0)=1\) and \(y(0)=-1\).
\begin{matlab}
%}
ts=linspace(0,5,21);
[xs,errs]=rk2int(fn,[1;-1],ts)
plot(ts,xs,ts,100*errs,'o')
print -depsc2 RKInt/rk2intPlot2.eps
%{
\end{matlab}

\begin{figure}
\centering
\begin{tabular}{c@{}c}
\rotatebox{90}{\hspace{10ex}$x,\ y,\ 100\times$step-errs} &
\includegraphics[width=0.85\linewidth]{RKInt/rk2intPlot2}\\
& $t$
\end{tabular}
\caption{example trajectory and errors}
\label{fig:2tpi}
\end{figure}

%}
