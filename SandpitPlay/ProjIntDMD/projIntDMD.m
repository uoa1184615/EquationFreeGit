% projIntDMD() is a basic projective integration of a
% given system of stiff deterministic ODEs.  AJR, Oct 2017
%!TEX root = ../Doc/equationFreeDoc.tex
%{
\subsection{\texttt{projIntDMD()}}
\label{sec:projIntDMD}
\localtableofcontents
\def\dmd{\textsc{dmd}}

This is a basic example of projective integration of a given system of stiff deterministic \ode{}s via \dmd, the Dynamic Mode Decomposition \cite[]{Kutz2016}.


\begin{matlab}
%}
function [xs,xss,tss]=projIntDMD(fun,x0,Ts,rank,dt,timeSteps)
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

\item \verb|timeSteps| is a two element vector: 
\begin{itemize}
\item \verb|timeSteps(1)| is the time thought to be needed for microscale simulation to reach the slow manifold; 
\item \verb|timeSteps(2)| is the subsequent time which \dmd\ analyses to model the slow manifold (must be longer than \(\verb|rank|\cdot\verb|dt|\)). 
\end{itemize}
\end{itemize}

\paragraph{Output}
\begin{itemize}
\item  \verb|xs|,  \(n\times\ell\) array of approximate solution vector at the specified times (the transpose of what \Matlab\ integrators do!)

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
If either of the timeSteps are non-integer valued, then assume they are both times, instead of micro-time-steps, so set the number of time-steps accordingly (multiples of~\verb|dt|).
\begin{matlab}
%}
timeSteps=round(timeSteps/dt); 
timeSteps(2)=max(rank+1,timeSteps(2));
%{
\end{matlab}

Set an algorithmic tolerance for miscellaneous purposes.
As at Jan 2018, it is a guess.  
It might be similar to some level of microscale `noise' in the burst.
Also, in an oscillatory mode for projection, set the expected maximum number of cycles in a projection.
\begin{matlab}
%}
algTol=log(1e8);
cycMax=3;
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
Reasons: the start-up time for implicit integrators, such as ode15s, is too onerous to be worthwhile for each short burst; the microscale time step needed for stability of explicit integrators is so small that a low order method is usually accurate enough.
\begin{matlab}
%}
x=[x0(:) nan(n,sum(timeSteps))];
for i=1:sum(timeSteps)
    xmid  =x(:,i)+dt/2*fun(Ts(k)+(i-1)*dt,x(:,i));
    x(:,i+1)=x(:,i)+dt*fun(Ts(k)+(i-0.5)*dt,xmid);
end
%{
\end{matlab}
If user requests microscale bursts, then store.
\begin{matlab}
%}
if nargout>1,xss=[xss x nan(n,1)];
if nargout>2,tss=[tss Ts(k)+(0:sum(timeSteps))*dt nan];
end,end
%{
\end{matlab}
Grossly check on whether the microscale integration is stable.
Use the 1-norm, the largest column sum of the absolute values, for little reason.
Is this any use??
\begin{matlab}
%}
if norm(x(j,ceil(end/2):end),1) ...
  > 10*norm(x(j,1:floor(end/2)),1)
  xMicroscaleIntegration=x, macroTime=Ts(k)
  warning('projIntDMD: microscale integration appears unstable')
  break%out of the integration loop
end
%{
\end{matlab}
Similarly if any non-numbers generated.
\begin{matlab}
%}
if sum(~isfinite(x(:)))>0
  break%ou of integration loop
end
%{
\end{matlab}


\paragraph{DMD extrapolation over the macroscale}
But skip if the simulation has already reached the next time.
\begin{matlab}
%}
iFin=1+sum(timeSteps);
DTgap=DT(k)-iFin*dt;
if DTgap*sign(dt)<=1e-9
   i=round(DT(k)/dt); x0(j)=x(:,i+1);
   else
%{
\end{matlab}

\dmd\ appears to work better when ones are adjoined to the data vectors.
\footnote{A reason is as follows.  
Consider the one variable linear \ode\ \(\dot x=f+J(x-x_0)\) with \(x(0)=x_0\) (as from a local linearisation of nonlinear \ode).
The solution is \(x(t)=(x_0-f/J)+(f/J)e^{Jt}\) which sampled at a time-step~\(\tau\) is \(x_k=(x_0-f/J)+(f/J)G^k\) for \(G:=e^{J\tau}\).
Then \(x_{k+1}\neq ax_k\) for any~\(a\).
However, \(\begin{bmat} x_{k+1}\\1 \end{bmat}
=\begin{bmat} G&a\\0&1 \end{bmat}
\begin{bmat} x_k\\1 \end{bmat}\) for a constant \(a:=(x_0-f/J)(1-G)\).
That is, with ones adjoined, the data from the \ode\ fits the \dmd\ approach.
}
\begin{matlab}
%}
iStart=1+timeSteps(1);
x=[x;ones(1,iFin)]; j1=[j;n+1];
%{
\end{matlab}
Then the basic \dmd\ algorithm: first the fit.
However, need to test whether we need to worry about the microscale time step being too small and leading to an effect analogous to `numerical differentiation' errors: 
akin to the rule-of-thumb in fitting chaos with time-delay coordinates that a good time-step is approximately the time of the first zero of the autocorrelation.
\begin{matlab}
%}
[U,S,V]=svd(x(j1,iStart:iFin-1),'econ');
S=diag(S);
Sr = S(1:rank); % singular values, rx1
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
omega = log(diag(D))/dt; % continuous-time eigenvalues, rx1
bFin=Phi\x(j1,iFin); % rx1
%{
\end{matlab}
But we want to neglect modes that are insignificant as characterised by \autoref{tbl:negbad}, or be warned of modes that grow too rapidly.
\begin{table}
\caption{\label{tbl:negbad}criterion for deciding if some \dmd\ modes are to be neglected, and if not neglected then are they growing too badly?}
\begin{equation*}
\begin{array}{llp{0.5\linewidth}}\hline
\text{neglectness}&
\parbox[t]{4em}{\raggedright range for \(\varepsilon\approx 10^{-8}\)}&
\text{reason}\\\hline
\max(0,-\log_e|b_i|)& 0-19 &
Very small noise in the burst implies a numerical error mode.
\\
\max(0,-\Re\omega_i\,\Delta t)& 0-19 &
Rapidly decaying mode of the macro-time-step is a micro-mode that happened to be resolved in the data.
\\\hline
\text{badness}&  &
provided not already neglected
\\\hline
\max(0,+\Re\omega_i\,\Delta t)& 0-19 &
Micro-scale mode that rapidly grows, so macro-step should be smaller.
\\
\frac3C|\Im\omega_i|\Delta t& 0-19 &
An oscillatory mode with\({}\geq C\) cycles in macro-step~\(\Delta t\).
\\\hline
\end{array}
\end{equation*}
\end{table}
Assume appropriate to sum the neglect-ness, and the badness, for testing. 
Then warn if there is a mode that is too bad.
\begin{matlab}
%}
DTgap=DT(k)-iFin*dt;
negness=max(0,-log(abs(bFin)))+max(0,-real(omega*DTgap));
badness=max(0,+real(omega*DTgap))+3/cycMax*abs(imag(omega))*DTgap;
iOK=find(negness<algTol);
iBad=find(badness(iOK)>algTol);
if ~isempty(iBad) 
    warning('projIntDMD: some bad modes in projection')
    badness=badness(iOK(iBad))
    rank=rank
    burstDt=timeSteps*dt
    break
    end
%{
\end{matlab}
Scatter the prediction into the non-Nan elements of \verb|x0|.
\begin{matlab}
%}
x0(j)=Phi(1:end-1,iOK)*(bFin(iOK).*exp(omega(iOK)*DTgap)); % nx1
%{
\end{matlab}
End the omission of the projection in the case when the burst is longer than the macroscale step.
\begin{matlab}
%}
end
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
