%Implementation of the 'legacy codes' approach to relaxing a given
%set of coordinates near the slow manifold. This scheme introduces non-
%trivial error if the fast dynamics is insufficiently stiff.
%JM, July 2018
%!TEX root = ../equationFreeDoc.tex
%{
\subsubsection{\texttt{cdmc()}}
\label{sec:cdmc}
\verb|cdmc()| iteratively applies the microsolver and then projects backwards in time to the initial conditions. The cumulative effect is to relax the variables to the attracting slow manifold, while keeping the final time for the output the same as the input time.
 \paragraph{Input}
\begin{itemize}
\item \verb|solver|, a black box microsolver suitable for PI. See any of \\{\verb|PIRK2(), PIRK4(), PIG()|} for a description of \verb|solver()|.
\item \verb|t|, an initial time
\item \verb|x|, an initial state
\item \verb|T|, a time period to apply \verb|solver()| for
\end{itemize}
\paragraph{Output}
\begin{itemize}
\item \verb|tout|, a vector of times. \verb|tout(end)| will equal \verb|t|.
\item \verb|xout|, an array of state estimates produced by \verb|solver()|. 
\end{itemize}
This function is a wrapper for the microsolver. For instance if the problem of interest is a dynamical system that is not too stiff, and which can be simulated by the solver \verb|sol(t,x,T)|, one would define\\ 
\verb|cSol = @(t,x,T) cdmc(sol,t,x,T)|,\\
 and thereafter use \verb|csol()| in place of \verb|sol()| as the solver for any PI scheme.
The original solver \verb|sol()| would create large errors if used in a PI scheme, but the output of \verb|cdmc()| will not.

\begin{matlab}
%}
function [tout, xout] = cdmc(solver,t,x,T)
%{
\end{matlab}
Begin with a standard application of the microsolver.
\begin{matlab}
%}
    [tt,xx]=feval(solver,t,x,T);
%{
\end{matlab}



Project backwards to before the initial time, then
simulate one burst forwards to obtain coordinates at t.
\begin{matlab}
%}
del = (xx(end,:) - xx(end-1,:))/(tt(end) - tt(end-1));
b_prop = 0.2; 
x = xx(end,:) + (1+b_prop)*(tt(1)-tt(end))*del;
tr=2*tt(1)-tt(end);
[tout,xout]=feval(solver,tr,x',b_prop*T);
%{
\end{matlab}
This concludes the function.
%}