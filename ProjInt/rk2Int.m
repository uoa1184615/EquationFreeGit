% rk2Int() is a simple example of Runge--Kutta, 2nd order,
% integration of a given deterministic ODE.  Used by
% PIG.m, PIGExample.m, PIGExplore.m
% AJR, 4 Apr 2019
%{
This is a simple example of Runge--Kutta, 2nd order,
integration of a given deterministic \ode.
\begin{matlab}
%}
function [ts,xs,errs] = rk2Int(dxdt,ts,x0)
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|dxdt()| is a function such as
\verb|dxdt=dxdt(t,x)| that computes the right-hand side of
the \ode\ \(d\xv/dt=\fv(\xv,t)\) where \xv~is a column
vector, say in \(\RR^n\) for \(n\geq1\)\,, \(t\)~is a
scalar, and the result~\fv\ is a column vector in~\(\RR^n\).

\item \verb|x0| is an \(\RR^n\) vector of initial values at
the time \verb|ts(1)|.

\item \verb|ts| is the begin and end times of the integration.
\end{itemize}

\paragraph{Output}
\begin{itemize}
\item \verb|ts|, vector of~$\ell$ times (guess $\ell=21$).

\item  \verb|xs|, array in \(\RR^{\ell\times n}\) of
approximate solution row vector at the specified times.
\end{itemize}

Compute the time-steps and create storage for outputs. Guess
that twenty time-steps is adequate.
\begin{matlab}
%}
ts = linspace(ts(1),ts(end),21).';
dt = diff(ts);
xs = nan(numel(x0),numel(ts));
%{
\end{matlab}
Initialise first result to the given initial condition, and
evaluate the initial time derivative into~\verb|f1|.
\begin{matlab}
%}
xs(:,1) = x0(:);
errs(1) = 0;
f1 = dxdt(ts(1),xs(:,1));
%{
\end{matlab}
Compute the time-steps from~\(t_k\) to~\(t_{k+1}\), copying
the derivative~\verb|f1| at the end of the last time-step to
be the derivative at the start of this one.
\begin{matlab}
%}
for k = 1:numel(dt)
  f0 = f1;
%{
\end{matlab}
Simple second-order accurate time-step.
\begin{matlab}
%}
  xh = xs(:,k)+f0*dt(k)/2;
  fh = dxdt(ts(k)+dt(k)/2,xh);
  xs(:,k+1) = xs(:,k)+fh*dt(k);
  f1 = dxdt(ts(k+1),xs(:,k+1));
end
xs = xs.';
%{
\end{matlab}
%}
