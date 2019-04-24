%Exploration of cdmc. JM, Feb 19.
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{Explore: Projective Integration using constraint-defined manifold computing}
\label{sec:Excdmc}
\localtableofcontents

In this example the Projective Integration-General scheme is
applied to a singularly perturbed ordinary differential
equation in which the time scale separation is not large.
The results demonstrate the value of the default
\verb|cdmc()| wrapper for the microsolver. 
\begin{devMan}
\begin{matlab}
%}
clear all, close all
%{
\end{matlab}

Set a weak time scale separation, and the underlying \ode{}s:
\begin{equation*}
\de t{x_1}=\cos x_1\sin x_2\cos t\,,
\quad
\de t{x_2}=\frac1\epsilon(-x_2+\cos x_1).
\end{equation*}
\begin{matlab}
%}
epsilon = 0.01;
dxdt=@(t,x) [ cos(x(1))*sin(x(2))*cos(t)
             (cos(x(1))-x(2))/epsilon ];
%{
\end{matlab}

Set the microsolver to be an integration using a standard
solver, and set the standard time of simulation for the
microsolver.
\begin{matlab}
%}
bT = epsilon*log(1/epsilon);
if ~exist('OCTAVE_VERSION','builtin')
    micB='ode45'; else micB='rk2Int'; end
microBurst = @(tb0, xb0) feval(micB,dxdt,[tb0 tb0+bT],xb0);
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


Simulate using \verb|PIG()|, first without the default
treatment of \verb|cdmc| for the microsolver and second
with. Generate a trusted solution using standard numerical
methods.
\begin{matlab}
%}
if ~exist('OCTAVE_VERSION','builtin')
    macInt='ode45'; else macInt='odeOct'; end
[nt,nx] = PIG(macInt,microBurst,tSpan,x0,[],[],'no cdmc');
[ct,cx] = PIG(macInt,microBurst,tSpan,x0);
[tClassic,xClassic] = feval(macInt,dxdt,tSpan,x0);
%{
\end{matlab}



\cref{fig:PIGE} plots the output. 
\begin{figure}
\centering
\caption{\label{fig:PIGE}Accurate simulation of a weakly
stiff non-autonomous system by \texttt{PIG()} using
\texttt{cdmc()}, and an inaccurate solution using a naive
application of~\texttt{PIG()}.}
\includegraphics[scale=0.9]{PIGExplore}
\end{figure}
\begin{matlab}
%}
figure
h = plot(nt,nx,'rx', ct,cx,'bo', tClassic,xClassic,'-');
legend(h(1:2:5),'Naive PIG','PIG + cdmc','Accurate')
xlabel('Time'), ylabel('State')
if ~exist('OCTAVE_VERSION','builtin')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 10])
%print('-depsc2','PIGExplore')
end
%{
\end{matlab}
A source of error in the standard \verb|PIG()| scheme is the
finite length of each burst,~\verb|bT|.  This computes a
time derivative at a time that is significantly different to
that requested by standard coded schemes.  Set \verb|bT| to
\verb|20*epsilon| or \verb|50*epsilon|\footnote{This example
is quite extreme: at \texttt{bT=50*epsilon}, it would be
computationally much cheaper to simulate the entire length
of tSpan using the microsolver alone.} to worsen the error
in both schemes. This example reflects a general principle:
most Projective Integration schemes incur a global error
term proportional to the burst time of the microsolver and
independent of the order of the microsolver. The
\verb|PIRKn()| schemes are written to eliminate this error,
but \verb|PIG()| works with any user-defined macrosolver and
cannot reduce this error, except by using the function
\verb|cdmc()|, its default. 
\end{devMan}
%}