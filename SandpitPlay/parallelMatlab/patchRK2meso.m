%patchRK2meso() is a simple example of Runge--Kutta, 2nd order, 
%integration of a given deterministic system on patches.
%AJR, 3 Sept 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchRK2meso()}}
\label{sec:patchRK2meso}
\localtableofcontents

This is a simple example of Runge--Kutta, 2nd order,
integration of a given deterministic \ode.

\begin{matlab}
%}
function [xs,errs] = patchRK2meso(fun,ts,x0,nMicro)
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|fun()| is a function such as
\verb|dxdt=fun(t,x)| that computes the right-hand side of
the \ode\ \(d\xv/dt=\fv(t,\xv)\) where \xv~is a
vector\slash array, \(t\)~is a
scalar, and the result~\fv\ is a correspondingly sized vector\slash array.
\item \verb|x0| is an vector\slash array of initial values at
the time \verb|ts(1)|.
\item \verb|ts| is a vector of meso-scale times to compute the
approximate solution, say in~\(\RR^\ell\) for
\(\ell\geq2\)\,.
\item \verb|nMicro|, optional, default~\(10\), is the number of sub-steps taken in each meso-scale step.
\end{itemize}

\paragraph{Output}
\begin{itemize}
\item  \verb|xs|, array in \(\RR^{\ell\times n}\) of
approximate solution row vector at the specified times.
\item \verb|errs|, column vector in \(\RR^{\ell}\) of error
estimate for  the step from \(t_{k-1}\) to~\(t_k\).
\end{itemize}

Set default number of micro-scale time-steps in each requested meso-scale step of~\verb|ts|.
\begin{matlab}
%}
if nargin<4, nMicro=10; end
%{
\end{matlab}

Compute the time-steps and create storage for outputs.
\begin{matlab}
%}
dt = diff(ts)/nMicro;
xs = cell(numel(ts),1);
errs = nan(numel(ts),1);
%{
\end{matlab}
Initialise first result to the given initial condition, and
evaluate the initial time derivative into~\verb|f1|.
\begin{matlab}
%}
xs{1} = x0;
errs(1) = 0;
f1 = fun(ts(1),xs{1});%??????????
%{
\end{matlab}
Compute the time-steps from~\(t_k\) to~\(t_{k+1}\), copying
the derivative~\verb|f1| at the end of the last time-step to
be the derivative at the start of this one.
\begin{matlab}
%}
for k = 1:numel(dt)
%{
\end{matlab}
Perform meso-scale time step with new interpolation of edge values.
\begin{matlab}
%}
  xs{k+1}=patchEdgeInt3(xs{k});
  for m=1:nMicro
    f0 = f1;
    assert(iscodistributed(f0),'f0 not codist')
%{
\end{matlab}
Simple second-order accurate micro-scale time-step.
\begin{matlab}
%}
    xh = xs{k+1}+f0*dt(k)/2;
    fh = patches.fun(ts(k)+dt(k)*(m-0.5),xh);
    xs{k+1} = xs{k}+fh*dt(k);
    f1 = patches.fun(ts(k)+m*dt(k),xs{k+1}); %???????
%{
\end{matlab}
End the micro-scale burst
\begin{matlab}
%}
  end
%{
\end{matlab}
Use the time derivative at~\(t_{k+1}\) to estimate an error
by storing the difference with what Simpson's rule would
estimate.
\begin{matlab}
%}
  f0=f0-2*fh+f1;
  assert(iscodistributed(f0),'f2ndDeriv not codist')
  errs(k+1) = sqrt(gather(mean(f0(:).^2)))*dt(k)/6;
end
%{
\end{matlab}
End of the function with results returned in~\verb|xs|
and~\verb|errs|.
%}
