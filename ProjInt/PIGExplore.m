%Exploration of cdmc. JM, Feb 19.
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{Explore: Projective Integration using constraint-defined manifold computing}
\label{sec:Excdmc}
\localtableofcontents

In this example the Projective Integration-General scheme is applied to a
singularly perturbed ordinary differential equation in which
the time scale separation is not large. The results demonstrate the value of
 the default \verb|cdmc()| wrapper for the
microsolver. 
\begin{devMan}
\begin{matlab}
%}
clear all, close all
%{
\end{matlab}

Set a weak time scale separation, and model.
\begin{matlab}
%}
epsilon = 0.01;
dxdt=@(t,x) [ cos(x(1))*sin(x(2))*cos(t)
             (cos(x(1))-x(2))/epsilon ];
%{
\end{matlab}

Set the microsolver to be an integration using a
standard solver, and set the standard time of simulation for
the microsolver.
\begin{matlab}
%}
bT = epsilon*log(1/epsilon);
microBurst = @(tb0,xb0) ode45(dxdt,[tb0 tb0+bT],xb0);
%{
\end{matlab}

Set initial conditions, and the time to be covered by the
macrosolver.  
\begin{matlab}
%}
x0 = [1 0];
tSpan=0:0.5:15;
%{
\end{matlab}


Simulate using \verb|PIG()|, first without the default treatment of \verb|cdmc|
for the microsolver and second with. Generate a trusted solution using standard
numerical methods.
\begin{matlab}
%}
[nt,nx] = PIG('ode45',microBurst,tSpan,x0,[],[],'no cdmc');
[ct,cx] = PIG('ode45',microBurst,tSpan,x0);
[t45,x45] = ode45(dxdt,tSpan([1 end]),x0);
%{
\end{matlab}



\cref{fig:PIGE} plots the output. 
\begin{figure}
\centering
\caption{Accurate simulation of a weakly stiff non-autonomous
system by \texttt{PIG()} using cdmc(), and an inaccurate solution
using a naive application of~\texttt{PIG()}.}\label{fig:PIGE}
\includegraphics[scale=0.9]{PIGExplore}
\end{figure}
\begin{matlab}
%}
figure
h = plot(nt,nx,'rx', ct,cx,'bo', t45,x45,'-');
legend(h(1:2:5),'Naive PIG','default: PIG + cdmc','Accurate')
xlabel('Time'), ylabel('State')
set(gcf,'PaperPosition',[0 0 14 10]), print('-depsc2','PIGExplore')
%{
\end{matlab}
The source of the error in
the standard \verb|PIG()| scheme is the burst length
\verb|bT|, that is significant on the slow time scale. Set
\verb|bT| to \verb|20*epsilon| or
\verb|50*epsilon|\footnote{this example is quite extreme: at
bT=50*epsilon, it would be computationally much cheaper to
simulate the entire length of tSpan using the microsolver
alone.} to worsen the error in both schemes. This example
reflects a general principle, that most Projective Integration schemes will
incur a global error term which is proportional to the
simulation time of the microsolver and independent of the
order of the microsolver. The \verb|PIRK()| schemes have
been written to minimise, if not eliminate entirely, this
error, but by design \verb|PIG()| works with any
user-defined macrosolver and cannot reduce this error. The
function \verb|cdmc()| reduces this error term by attempting
to mimic the microsolver without advancing time. 
\end{devMan}
%}