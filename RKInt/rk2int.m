%rk2int() is a simple example of Runge--Kutta, 2nd order, 
%integration of a given deterministic ODE.
%AJR, 29 Mar 2017
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{\texttt{rk2int()}}
\label{sec:rk2int}

This is a simple example of Runge--Kutta, 2nd order, integration of a given deterministic \ode.

\begin{matlab}
%}
function [xs,errs]=rk2int(fun,x0,ts)
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|fun()| is a function such as \verb|dxdt=fun(x,t)| that computes the right-hand side of the \ode\ \(d\xv/dt=\fv(\xv,t)\) where \xv~is a column vector, say in \(\RR^n\) for \(n\geq1\)\,, \(t\)~is a scalar, and the result~\fv\ is a column vector in~\(\RR^n\).
\item \verb|x0| is an \(\RR^n\) vector of initial values at the time \verb|ts(1)|.
\item \verb|ts| is a vector of times to compute the approximate solution, say in~\(\RR^\ell\) for \(\ell\geq2\)\,.
\end{itemize}

\paragraph{Output}
\begin{itemize}
\item  \verb|xs|, array in \(\RR^{n\times\ell}\) of approximate solution vector at the specified times---but the transpose of what \Matlab\ does!
\item \verb|errs|, vector in \(\RR^{1\times\ell}\) of error estimate for  the step from \(t_{k-1}\) to~\(t_k\).
\end{itemize}

Compute the time steps and create storage for outputs.
\begin{matlab}
%}
dt=diff(ts);
xs=nan(length(x0),length(ts));
errs=nan(1,length(ts));
%{
\end{matlab}
Initialise first result to the given initial condition, and evaluate the initial time derivative into~\verb|f1|.
\begin{matlab}
%}
xs(:,1)=x0(:);
errs(1)=0;
f1=fun(xs(:,1),ts(1));
%{
\end{matlab}
Compute the time-steps from~\(t_k\) to~\(t_{k+1}\), copying the derivative~\verb|f1| at the end of the last time-step to be the derivative at the start of this one.
\begin{matlab}
%}
for k=1:length(dt)
  f0=f1;
%{
\end{matlab}
Simple second order accurate time step.
\begin{matlab}
%}
  xh=xs(:,k)+f0*dt(k)/2;
  fh=fun(xh,ts(k)+dt(k)/2);
  xs(:,k+1)=xs(:,k)+fh*dt(k);
  f1=fun(xs(:,k+1),ts(k+1));
%{
\end{matlab}
Use the time derivative at~\(t_{k+1}\) to estimate an error by storing  the difference with what Simpson's rule would estimate.
\begin{matlab}
%}
  errs(k+1)=norm(f0-2*fh+f1)*dt(k)/6;
end
%{
\end{matlab}
End of the function with results returned in~\verb|xs| and~\verb|errs|.
%}
