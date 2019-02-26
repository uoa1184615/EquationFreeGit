% Basic example of PIG. JM, Sept 18.
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{Example: Projective Integration using General macrosolvers }
\label{sec:ExPIG}
In this example the Projective Integration-General scheme is
applied to a singularly perturbed ordinary differential
equation. The aim is to use a standard non-stiff numerical
integrator, such as \verb|ode45()|, on the slow, long-time
macroscale.  For this stiff system, \verb|PIG()| is an order
of magnitude faster than ordinary use of~\verb|ode45|.

\begin{devMan}
\begin{matlab}
%}
clear all, close all
%{
\end{matlab}

Set time scale separation and model.
\begin{matlab}
%}
epsilon = 1e-4;
dxdt=@(t,x) [ cos(x(1))*sin(x(2))*cos(t)
             (cos(x(1))-x(2))/epsilon ];
%{
\end{matlab}

Set the `black-box' microsolver to be an integration using a
standard solver, and set the standard time of simulation for
the microsolver.
\begin{matlab}
%}
bT = epsilon*log(1/epsilon);
microBurst = @(tb0, xb0) ode45(dxdt,[tb0 tb0+bT],xb0);
%{
\end{matlab}

Set initial conditions, and the time to be covered by the
macrosolver. 
\begin{matlab}
%}
x0 = [1 1.4];
tSpan=[0 15];
%{
\end{matlab}
Now time and integrate the above system over \verb|tspan|
using \verb|PIG()| and, for comparison, a brute force
implementation of \verb|ode45()|. Report the time taken by
each method.
\begin{matlab}
%}
tic
[ts,xs,tms,xms] = PIG('ode45',microBurst,tSpan,x0);
tPIGusingODE45asMacro = toc
tic
[t45,x45] = ode45(dxdt,tSpan,x0);
tODE45alone = toc
%{
\end{matlab}


Plot the output on two figures, showing the truth and
macrosteps on both, and all applications of the microsolver
on the first figure.
\begin{matlab}
%}
figure
h = plot(ts,xs,'o', t45,x45,'-', tms,xms,'.');
legend(h(1:2:5),'PI Solution','ode45 Solution','PI microsolver')
xlabel('Time'), ylabel('State')

figure
h = plot(ts,xs,'o', t45,x45,'-');
legend(h([1 3]),'PI Solution','ode45 Solution')
xlabel('Time'), ylabel('State')
set(gcf,'PaperPosition',[0 0 14 10]), print('-depsc2','PIGExample')
%{
\end{matlab}
\cref{fig:PIGExample} plots the output.
\begin{figure}\centering
\caption{Accurate simulation of a stiff nonautonomous system
by PIG().  The microsolver is called on-the-fly by the
macrosolver \texttt{ode45}.}\label{fig:PIGExample}
\includegraphics[scale=0.8]{../ProjInt/PIGExample}
\end{figure}

\begin{itemize}
\item The problem may be made more, or less, stiff by
changing the time-scale separation parameter
\(\epsilon=\verb|epsilon|\). The compute time of
\verb|PIG()| is almost independent of~\(\epsilon\), whereas
that of \verb|ode45()| is proportional to~\(1/\epsilon\).

But if the problem is insufficiently stiff
(larger~\(\epsilon\)), then \verb|PIG()| produces nonsense.
This nonsense is overcome by \verb|cdmc()|
(\cref{sec:Excdmc}).

\item The mildly stiff problem in \cref{sec:ExPIRK} may be
efficiently solved by a standard solver (e.g.,
\verb|ode45()|). The stiff but low dimensional problem in
this example can be solved efficiently by a standard stiff
solver (e.g., \verb|ode15s()|). The real advantage of the
Projective Integration schemes is in high dimensional stiff
problems, that cannot be efficiently solved by most standard
methods.
\end{itemize}

\end{devMan}
%}