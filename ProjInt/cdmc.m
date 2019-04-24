% Relax a given initial condition to one onto the slow
% manifold by two steps of the 'xmas-tree' algorithm.
% JM & AJR, July 2018 -- Apr 2019
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{cdmc()}: constraint defined manifold computing}
\label{sec:cdmc}

The function \verb|cdmc()| iteratively applies the given
micro-burst and then projects backward to the initial time.
The cumulative effect is to relax the variables to the
attracting slow manifold, while keeping the `final' time for
the output the same as the input time.

\begin{matlab}
%}
function [ts, xs] = cdmc(microBurst, t0, x0)
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|microBurst()|, a black-box micro-burst function
suitable for Projective Integration. See any of
\verb|PIRK2()|, \verb|PIRK4()|, or \verb|PIG()| for a
description of \verb|microBurst()|.
\item \verb|t0|, an initial time.
\item \verb|x0|, an initial state vector.
\end{itemize}
\paragraph{Output}
\begin{itemize}
\item \verb|ts|, a vector of times.
\item \verb|xs|, an array of state estimates produced by
\verb|microBurst()|. 
\end{itemize}


This function is a wrapper for the micro-burst. For instance
if the problem of interest is a dynamical system that is not
too stiff, and which is simulated by the micro-burst function
\verb|sol(t,x)|, one would invoke \verb|cdmc()| by defining
\begin{verbatim}
cdmcSol = @(t,x) cdmc(sol,t,x)|
\end{verbatim}
and thereafter use \verb|cdmcSol()| in place of \verb|sol()|
as the microBurst in any Projective Integration scheme. The
original microBurst \verb|sol()| could create large errors
if used in the \verb|PIG()| scheme, but the output via
\verb|cdmc()| should not.

\begin{devMan}
Begin with a standard application of the micro-burst. Need
\verb|feval| as \verb|microBurst| has multiple outputs.
\begin{matlab}
%}
[t1,x1] = feval(microBurst,t0,x0);
bT = t1(end)-t1(1);
%{
\end{matlab}

Project backwards to before the initial time, then simulate
just one burst forward to obtain a simulation burst that
ends at the original~\verb|t0|.
\begin{matlab}
%}
dxdt = (x1(end,:) - x1(end-1,:))/(t1(end) - t1(end-1));
x0 = x1(end,:)-2*bT*dxdt;
t0 = t1(1)-bT;
[t2,x2] = feval(microBurst,t0,x0.');
%{
\end{matlab}
Return both sets of output(?), although only \verb|(t2,x2)|
should be used in Projective Integration---maybe safer to
return only \verb|(t2,x2)|.
\begin{matlab}
%}
ts = [t1(:); t2(:)];
xs = [x1; x2];
%{
\end{matlab}
\end{devMan}
%}