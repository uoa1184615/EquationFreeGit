% projInt1() is a basic projective integration of a
% given system of stiff deterministic ODEs.  AJR, Oct 2017
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{\texttt{projInt1()}}
\label{sec:projInt1}
\def\dmd{\textsc{dmd}}

This is a basic example of projective integration of a given system of stiff deterministic \ode{}s via \dmd, the Dynamic Mode Decomposition \cite[]{Kutz2016}.


\begin{matlab}
%}
function [xs,xss,tss]=projInt1(fun,x0,Ts,rank,dt,nTimeSteps)
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|fun()| is a function such as \verb|dxdt=fun(t,x)| that computes the right-hand side of the \ode\ \(d\xv/dt=\fv(t,\xv)\) where \xv~is a column vector, say in \(\RR^n\) for \(n\geq1\)\,, \(t\)~is a scalar, and the result~\fv\ is a column vector in~\(\RR^n\).
\item \verb|x0| is an \(n\)-vector of initial values at the time \verb|ts(1)|.  
If any entries in~\verb|x0| are~\verb|NaN|, then \verb|fun()| must cope, and only the non-\verb|NaN| components are projected in time.
\item \verb|Ts| is a vector of times to compute the approximate solution, say in~\(\RR^\ell\) for \(\ell\geq2\)\,.
\item \verb|rank| is the rank of the \dmd\ extrapolation over macroscale time steps.  
Suspect \verb|rank| should be at least one more than the effective number of slow variables.
\item \verb|dt| is the size of the microscale time-step.  Must be small enough so that RK2 integration of the \ode{}s is stable.
\item \verb|nTimeSteps| is a two element vector: 
\begin{itemize}
\item \verb|nTimeSteps(1)| is the number of microscale time-steps~\verb|dt| thought to be needed for microscale simulation to reach the slow manifold; 
\item \verb|nTimeSteps(2)|, must be bigger than~\verb|rank|, is the number of subsequent time-steps to take which \dmd\ analyses to mode the slow manifold. 
\end{itemize}
\end{itemize}

\paragraph{Output}
\begin{itemize}
\item  \verb|xs|,  \(n\times\ell\) array of approximate solution vector at the specified times (the transpose of what \Matlab\ integrators do!)
%\item \verb|errs|, vector in \(\RR^{1\times\ell}\) of error estimate for  the step from \(t_{k-1}\) to~\(t_k\).
\item \verb|xss|, optional, \(n\times\text{big}\) array of the microscale simulation bursts---separated by NaNs for possible plotting.
\item \verb|tss|, optional, \(1\times\text{big}\) vector of times corresponding to the columns of~\verb|xss|.
\end{itemize}

Compute the time steps and create storage for outputs.
\begin{matlab}
%}
DT=diff(Ts);
n=length(x0);
xs=nan(n,length(Ts));
xss=[];tss=[];
%{
\end{matlab}
If any~\verb|x0| are \verb|NaN|, then assume the time derivative routine can cope, and here we just exclude these from \dmd\ projection and from any error estimation.
This allows a user to have space in the solutions for breaks in the data vector (that, for example, may be filled in with boundary values for a \pde\ discretisation).
\begin{matlab}
%}
j=find(~isnan(x0));
%{
\end{matlab}
Initialise first result to the given initial condition.
\begin{matlab}
%}
xs(:,1)=x0(:);
%{
\end{matlab}
Projectively integrate each of the time-steps from~\(t_k\) to~\(t_{k+1}\).
\begin{matlab}
%}
for k=1:length(DT)
%{
\end{matlab}
Microscale integration is simple, second order, Runge--Kutta method.
\begin{matlab}
%}
x=[x0(:) nan(n,sum(nTimeSteps))];
for i=1:sum(nTimeSteps)
    xh=x(:,i)+dt/2*fun(Ts(k)+(i-1)*dt,x(:,i));
    x(:,i+1)=x(:,i)+dt*fun(Ts(k)+(i-0.5)*dt,xh);
end
%{
\end{matlab}
If user requests microscale bursts, then store.
\begin{matlab}
%}
if nargout>1,xss=[xss x nan(n,1)];
if nargout>2,tss=[tss Ts(k)+(0:sum(nTimeSteps))*dt nan];
end,end
%{
\end{matlab}
Grossly check on whether the microscale integration is stable.
Is this any use??
\begin{matlab}
%}
if norm(x(j,nTimeSteps(1)+(1:nTimeSteps(2)))) ...
  > 3*norm(x(j,1:nTimeSteps(1)))
  xMicroscaleIntegration=x, macroTime=Ts(k)
  error('projInt1: microscale integration appears unstable')
end
%{
\end{matlab}

\paragraph{DMD extrapolation over the macroscale}
\dmd\ appears to work better when ones are adjoined to the data vectors for some unknown reason.
\begin{matlab}
%}
iFin=1+sum(nTimeSteps);
iStart=1+nTimeSteps(1);
x=[x;ones(1,iFin)]; j1=[j;n+1];
%{
\end{matlab}
Then the basic \dmd\ algorithm: first the fit.
\begin{matlab}
%}
[U,S,V]=svd(x(j1,iStart:iFin-1),'econ');
S=diag(S);
Sr = S(1:rank) % rx1
AUr=bsxfun(@rdivide,x(j1,iStart+1:iFin)*V(:,1:rank),Sr.');%nxr
Atilde = U(:,1:rank)'*AUr; % low-rank dynamics, rxr
[Wr, D] = eig(Atilde); % rxr
Phi = AUr*Wr; % DMD modes, nxr
%{
\end{matlab}
Second, reconstruct a prediction for the time step.
The current micro-simulation time is \verb|dt*iFin|, so step forward an amount to predict the systems state at \verb|Ts(k+1)|.
Perhaps should test~\(\omega\) and abort if 'large' and/or positive??
Answer: not necessarily as if the rank is large then the omega could contain large negative values.
\begin{matlab}
%}
omega = log(diag(D))/dt % continuous-time eigenvalues, rx1
bFin=Phi\x(j1,iFin); % rx1
x0(j)=Phi(1:end-1,:)*(bFin.*exp(omega*(DT(k)-iFin*dt))); % nx1
%{
\end{matlab}
Since some of the~\(\omega\) may be complex, if the simulation burst is real, then force the \dmd\ prediction to be real.
\begin{matlab}
%}
if isreal(x), x0=real(x0); end
xs(:,k+1)=x0;  
%{
\end{matlab}

End the macroscale time stepping.
\begin{matlab}
%}
end
%{
\end{matlab}

If requested, then add the final point to the microscale data.
\begin{matlab}
%}
if nargout>1,xss=[xss x0];
if nargout>2,tss=[tss Ts(end)];
end,end
%{
\end{matlab}
End of the function with result vectors returned in columns of~\verb|xs|, one column for each time in~\verb|Ts|.
%}
