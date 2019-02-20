%Implementation of the 'legacy codes' approach to relaxing a given
%set of coordinates near the slow manifold. This scheme introduces non-
%trivial error if the fast dynamics is insufficiently stiff.
%JM, July 2018
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{cdmc()}}
\label{sec:cdmc}
\verb|cdmc()| iteratively applies the micro-burst and then projects backwards in time to the initial conditions. The cumulative effect is to relax the variables to the attracting slow manifold, while keeping the final time for the output the same as the input time.

\begin{matlab}
%}
function [ts, xs] = cdmc(microBurst,t0,x0)
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|microBurst()|, a black box micro-burst function suitable for Projective Integration. See any of \verb|PIRK2()|, \verb|PIRK4()|, or \verb|PIG()| for a description of \verb|microBurst()|.
\item \verb|t0|, an initial time
\item \verb|x0|, an initial state
\end{itemize}
\paragraph{Output}
\begin{itemize}
\item \verb|ts|, a vector of times. \verb|tout(end)| will equal \verb|t|.
\item \verb|xs|, an array of state estimates produced by \verb|microBurst()|. 
\end{itemize}


This function is a wrapper for the micro-burst. For instance if the problem of interest is a dynamical system that is not too stiff, and which can be simulated by the microBurst \verb|sol(t,x,T)|, one would define
\begin{verbatim}
cSol = @(t,x) cdmc(sol,t,x)|
\end{verbatim}
and thereafter use \verb|csol()| in place of \verb|sol()| as the microBurst for any Projective Integration scheme.
The original microBurst \verb|sol()| could create large errors if used in a Projective Integration scheme, but the output of \verb|cdmc()| should not.

\begin{devMan}
Begin with a standard application of the micro-burst.
\begin{matlab}
%}
[ts,xs] = feval(microBurst,t0,x0);
bT = ts(end)-ts(1);
%{
\end{matlab}

Project backwards to before the initial time, then
simulate just one burst forward to obtain a simulation burst that ends at the original~\verb|t0|.
\begin{matlab}
%}
dxdt = (xs(end,:) - xs(end-1,:))/(ts(end,:) - ts(end-1,:));
x0 = xs(end,:)-2*bT*dxdt;
t0 = ts(1)-bT;
[ts,xs] = feval(microBurst,t0,x0.');
%{
\end{matlab}
\end{devMan}
%}