% PIRK4 implements fourth-order Projective Integration with
% a user-specified microsolver.  The macrosolver adapts the
% explicit fourth-order Runge--Kutta scheme. JM, Oct 2018. 
%!TEX root = ../Doc/eqnFreeDevMan.tex

%{
\section{\texttt{PIRK4()}: projective integration of fourth-order accuracy}
\label{sec:PIRK4}
\localtableofcontents

\subsection{Introduction}

This Projective Integration scheme implements a macrosolver
analogous to the fourth-order Runge--Kutta method.

\begin{matlab}
%}
function [x, tms, xms, rm, svf] = PIRK4(microBurst, tSpan, x0, bT)
%{
\end{matlab}

See \cref{sec:PIRK2} as the inputs and outputs are the 
same as \verb|PIRK2()|.

\paragraph{If no arguments, then execute an example}
%\label{sec:pirk4eg}
\begin{matlab}
%}
if nargin==0
%{
\end{matlab}
\subparagraph{Example of Michaelis--Menton backwards in time} 
The Michaelis--Menten enzyme kinetics is expressed as a
singularly perturbed system of differential equations for
\(x(t)\) and~\(y(t)\) (encoded in function \verb|MMburst|):
\begin{equation*}
\frac{dx}{dt}=-x+(x+\tfrac12)y \quad\text{and}\quad
\frac{dy}{dt}=\frac1\epsilon\big[x-(x+1)y\big].
\end{equation*}
With initial conditions \(x(0)=y(0)=0.2\), the following
code uses forward time bursts in order to integrate
backwards in time to \(t=-5\). It plots the computed
solution over time \(-5\leq t\leq0\) for parameter
\(\epsilon=0.1\)\,. Since the rate of decay is
\(\beta\approx 1/\epsilon\) we choose a burst length
\(\epsilon\log(|\Delta|/\epsilon)\) as here the macroscale
time-step \(\Delta=-1\).
\begin{matlab}
%}
global MMepsilon
MMepsilon = 0.1
ts = 0:-1:-5
bT = MMepsilon*log(abs(ts(2)-ts(1))/MMepsilon)
[x,tms,xms,rm,svf] = PIRK4(@MMburst, ts, 0.2*[1;1], bT);
figure, plot(ts,x,'o:',tms,xms)
xlabel('time t'), legend('x(t)','y(t)')
title('Backwards-time projective integration of Michaelis--Menten')
%{
\end{matlab}
Upon finishing execution of the example, exit this function.
\begin{matlab}
%}
return
end%if no arguments
%{
\end{matlab}

\input{../ProjInt/MMburst.m}



\begin{devMan}


\paragraph{Input}
\begin{itemize}
\item \verb|microBurst()|, a function that produces output from
the user-specified code for microscale simulation. 
\begin{verbatim}
[tOut, xOut] = microBurst(tStart, xStart, bT)
\end{verbatim}
\begin{itemize}
\item Inputs:
  \verb|tStart|,~the start time of a burst of simulation;
  \(\verb|xStart|\),~the row \(n\)-vector of the starting
  state; \verb|bT|, optional, the total time to simulate in
  the burst---if \verb|microBurst()| determines~\verb|bT|,
  then replace~\verb|bT| in the argument list
  by~\verb|varargin|. 
\item Outputs:
  \verb|tOut|,~the column vector of solution times; and 
  \verb|xOut|,~an array in which each \emph{row} contains
  the system state at corresponding times.
\end{itemize}

\item \verb|tSpan| is an \(\ell\)-vector of times at which
the user requests output, of which the first element is
always the initial time. \verb|PIRK4()| does not use
adaptive time-stepping; the macroscale time-steps are
(nearly) the steps between elements of \verb|tSpan|.

\item \verb|x0| is an \(n\)-vector of initial values at the
initial time \verb|tSpan(1)|. Elements of~\verb|x0| may be
\verb|NaN|: they are included in the simulation and output,
and often represent boundaries in space fields.

\item \verb|bT|, optional, either missing, or
empty~(\verb|[]|), or a scalar: if a given scalar, then it
is the length of the micro-burst simulations---the minimum
amount of time needed for the microscale simulation to relax
to the slow manifold; else if missing or~\verb|[]|, then
\verb|microBurst()| must itself determine the length of a
computed burst.
\begin{matlab}
%}
if nargin<4, bT=[]; end
%{
\end{matlab}

\end{itemize}


\paragraph{Output}
If there are no output arguments specified, then a plot is 
drawn of the computed solution~\verb|x| versus \verb|tSpan|.
\begin{itemize}
\item  \verb|x|, an \(\ell \times n\) array of the
approximate solution vector. Each row is an estimated state
at the corresponding time in \verb|tSpan|.  The simplest
usage is then \verb|x = PIRK4(microBurst,tSpan,x0,bT)|.

However, microscale details of the underlying Projective
Integration computations may be helpful. \verb|PIRK4()|
provides two to four optional outputs of the microscale
bursts. 

\item \verb|tms|, optional, is an \(L\) dimensional column
vector containing microscale times of burst simulations,
each burst separated by~\verb|NaN|; 

\item \verb|xms|,
optional, is an \(L\times n\) array of the corresponding
microscale states---this data is an accurate simulation of
the state and may help visualise more details of the
solution. 

\item \verb|rm|, optional, a struct containing the
`remaining' applications of the micro-burst required by the
Projective Integration method during the calculation of the
macrostep: \begin{itemize}
\item \verb|rm.t|~is a column vector of microscale times; and 
\item \verb|rm.x|~is the array of corresponding burst states.
\end{itemize}
The states \verb|rm.x| do not have the same physical
interpretation as those in \verb|xms|; the \verb|rm.x| are
required in order to estimate the slow vector field during
the calculation of the Runge--Kutta increments, and do not
in general resemble the true dynamics.

\item  \verb|svf|, optional, a struct containing the
Projective Integration estimates of the slow vector field.
\begin{itemize}
\item \verb|svf.t| is a \(4\ell\) dimensional column vector
containing all times at which the Projective Integration
scheme is extrapolated along micro-burst data to form a
macrostep. 
\item \verb|svf.dx| is a \(4\ell\times n\) array containing
the estimated slow vector field.
\end{itemize}

\end{itemize}


\subsection{The projective integration code}

Determine the number of time-steps and preallocate storage
for macroscale estimates.
\begin{matlab}
%}
nT=length(tSpan); 
x=nan(nT,length(x0)); 
%{
\end{matlab}

Get the number of expected outputs and set logical indices
to flag what data should be saved.
\begin{matlab}
%}
nArgs=nargout(); 
saveMicro = (nArgs>1); 
saveFullMicro = (nArgs>3); 
saveSvf = (nArgs>4); 
%{
\end{matlab}


Run a preliminary application of the micro-burst on the
initial conditions to help relax to the slow manifold. This
is done in addition to the micro-burst in the main loop,
because the initial conditions are often far from the
attracting slow manifold. Require the user to input and
output rows of the system state.
\begin{matlab}
%}
x0 = reshape(x0,1,[]); 
[relax_t,relax_x0] = microBurst(tSpan(1),x0,bT);
%{
\end{matlab}

Use the end point of the micro-burst as the initial
conditions.
\begin{matlab}
%}
tSpan(1) = relax_t(end); 
x(1,:)=relax_x0(end,:); 
%{
\end{matlab}

If saving information, then record the first application of
the micro-burst. Allocate cell arrays for times and states
for outputs requested by the user, as concatenating cells is
much faster than iteratively extending arrays. 
\begin{matlab}
%}
if saveMicro 
    tms = cell(nT,1);
    xms = cell(nT,1);
    tms{1} = reshape(relax_t,[],1);
    xms{1} = relax_x0;
    if saveFullMicro 
        rm.t = cell(nT,1);
        rm.x = cell(nT,1);
        if saveSvf 
            svf.t = nan(4*nT-4,1); 
            svf.dx = nan(4*nT-4,length(x0)); 
        end
    end
end
%{
\end{matlab}


\paragraph{Loop over the macroscale time-steps}
\begin{matlab}
%}
for jT = 2:nT
    T = tSpan(jT-1);
%{
\end{matlab}
If four applications of the micro-burst would cover the
entire macroscale time-step, then do so (setting some
internal states to \verb|NaN|); else proceed to projective
step.
\begin{matlab}
%}
    if ~isempty(bT) && 4*abs(bT)>=abs(tSpan(jT)-T) && bT*(tSpan(jT)-T)>0
        [t1,xm1] = microBurst(T, x(jT-1,:), tSpan(jT)-T);
        x(jT,:) = xm1(end,:);
        t2=nan; xm2=nan(1,size(xm1,2));
        t3=nan; t4=nan; xm3=xm2; xm4 = xm2; dx1=xm2; dx2=xm2;
    else
%{
\end{matlab}

Run the first application of the micro-burst; since this
application directly follows from the initial conditions, or
from the latest macrostep, this microscale information is
physically meaningful as a simulation of the system. Extract
the size of the final time-step. 
   \begin{matlab}
%}
    [t1,xm1] = microBurst(T, x(jT-1,:), bT);
    del = t1(end)-t1(end-1);
%{
\end{matlab}
Check for round-off error.
\begin{matlab}
%}
    xt=[reshape(t1(end-1:end),[],1) xm1(end-1:end,:)];
    roundingTol=1e-8;
    if norm(diff(xt))/norm(xt,'fro') < roundingTol
    warning(['significant round-off error in 1st projection at T=' num2str(T)])
    end
%{
\end{matlab}

Find the needed time-step to reach time \verb|tSpan(n+1)|
and form a first estimate \verb|dx1| of the slow vector
field.
\begin{matlab}
%}
    Dt = tSpan(jT)-t1(end);  
    dx1 = (xm1(end,:)-xm1(end-1,:))/del; 
%{
\end{matlab}
Assume burst times are the same length for this macro-step,
or effectively so (recall that \verb|bT| may be empty as it
may be only coded and known in \verb|microBurst()|).
\begin{matlab}
%}
abT = t1(end)-t1(1);
%{
\end{matlab}


Project along \verb|dx1| to form an intermediate
approximation of \verb|x|; run another application of the
micro-burst and form a second estimate of the slow vector
field.  
\begin{matlab}
%}
    xint = xm1(end,:) + (Dt/2-abT)*dx1;
    [t2,xm2] = microBurst(T+Dt/2, xint, bT);
    del = t2(end)-t2(end-1);
    dx2 = (xm2(end,:)-xm2(end-1,:))/del; 
    
    xint = xm1(end,:) + (Dt/2-abT)*dx2;
    [t3,xm3] = microBurst(T+Dt/2, xint, bT);
    del = t3(end)-t3(end-1);
    dx3 = (xm3(end,:)-xm3(end-1,:))/del; 
    
    xint = xm1(end,:) + (Dt-abT)*dx3;
    [t4,xm4] = microBurst(T+Dt, xint, bT);
    del = t4(end)-t4(end-1);
    dx4 = (xm4(end,:)-xm4(end-1,:))/del; 
%{
\end{matlab}
Check for round-off error.
\begin{matlab}
%}
    xt=[reshape(t2(end-1:end),[],1) xm2(end-1:end,:)];
    if norm(diff(xt))/norm(xt,'fro') < roundingTol
    warning(['significant round-off error in 2nd projection at T=' num2str(T)])
    end
%{
\end{matlab}

Use the weighted average of the estimates of the slow vector
field to take a macrostep.
\begin{matlab}
%}
    x(jT,:) = xm1(end,:) + Dt*(dx1 + 2*dx2 + 2*dx3 + dx4)/6; 
%{
\end{matlab}

Now end the if-statement that tests whether a projective
step saves simulation time.
\begin{matlab}
%}
    end
%{
\end{matlab}

If saving trusted microscale data, then populate the cell
arrays for the current loop iterate with the time-steps and
output of the first application of the micro-burst. Separate
bursts by~\verb|NaN|s.
\begin{matlab}
%}
    if saveMicro 
        tms{jT} = [reshape(t1,[],1); nan];
        xms{jT} = [xm1; nan(1,size(xm1,2))];
%{
\end{matlab}

If saving all microscale data, then repeat for the remaining
applications of the micro-burst.         
\begin{matlab}
%}
        if saveFullMicro 
            rm.t{jT} = [reshape(t2,[],1); nan;...
                        reshape(t3,[],1); nan;...
                        reshape(t4,[],1); nan];
            rm.x{jT} = [xm2; nan(1,size(xm2,2));...
                        xm3; nan(1,size(xm2,2));...
                        xm4; nan(1,size(xm2,2))];
%{
\end{matlab}

If saving Projective Integration estimates of the slow
vector field, then populate the corresponding cells with
times and estimates.
\begin{matlab}
%}
            if saveSvf 
                svf.t(4*jT-7:4*jT-4) = [t1(end); t2(end); t3(end); t4(end)];
                svf.dx(4*jT-7:4*jT-4,:) = [dx1; dx2; dx3; dx4];
            end
        end
    end
%{
\end{matlab}
Terminate the main loop:
\begin{matlab}
%}
end
%{
\end{matlab}


Overwrite \verb|x(1,:)| with the specified initial condition
\verb|tSpan(1)|.
\begin{matlab}
%}
x(1,:) = reshape(x0,1,[]); 
%{
\end{matlab}

For additional requested output, concatenate all the cells
of time and state data into two arrays. 
\begin{matlab}
%}
if saveMicro 
    tms = cell2mat(tms);
    xms = cell2mat(xms);
    if saveFullMicro
        rm.t = cell2mat(rm.t);
        rm.x = cell2mat(rm.x);
    end
end
%{
\end{matlab}


\subsection{If no output specified, then plot simulation}
\begin{matlab}
%}
if nArgs==0
    figure, plot(tSpan,x,'o:')
    title('Projective Simulation with PIRK4')
end
%{
\end{matlab}

This concludes \verb|PIRK4()|.
\begin{matlab}
%}
end
%{
\end{matlab}
\end{devMan}
%}

