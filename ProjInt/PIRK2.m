% PIRK2() implements second-order Projective Integration
% with a user-specified microsolver.  The macrosolver adapts
% the explicit second-order Runge--Kutta Improved Euler
% scheme. JM and AJR, Oct 2018 -- Oct 2020.  Execute with no
% arguments to see an example.  
%!TEX root = ../Doc/eqnFreeDevMan.tex

%{
\section{\texttt{PIRK2()}: projective integration of second-order accuracy}
\label{sec:PIRK2}
\localtableofcontents

\subsection{Introduction}

This Projective Integration scheme implements a macroscale
scheme that is analogous to the second-order Runge--Kutta
Improved Euler integration.

\begin{matlab}
%}
function [x, tms, xms, rm, svf] = PIRK2(microBurst, tSpan, x0, bT)
%{
\end{matlab}

\paragraph{Input}
If there are no input arguments, then this function applies
itself to the Michaelis--Menton example: see the code in
\cref{sec:pirk2eg} as a basic template of how to use.
\begin{itemize}
\item \verb|microBurst()|, a user-coded function that
computes a short-time burst of the microscale simulation. 
\begin{verbatim}
[tOut, xOut] = microBurst(tStart, xStart, bT)
\end{verbatim}
\begin{itemize}
\item Inputs:
  \verb|tStart|,~the start time of a burst of simulation;
  \(\verb|xStart|\),~the row \(n\)-vector of the starting
  state; \verb|bT|, \emph{optional}, the total time to
  simulate in the burst---if your \verb|microBurst()|
  determines the burst time, then replace~\verb|bT| in the
  argument list by~\verb|varargin|. 
\item Outputs:
  \verb|tOut|,~the column vector of solution times; and 
  \verb|xOut|,~an array in which each \emph{row} contains
  the system state at corresponding times.
\end{itemize}
Be wary that for very large scale separations (such as
\verb|MMepsilon<1e-5| in the Michaelis--Menten example),
microscale integration by error-controlled variable-step
routines (such as \verb|ode23/45|) often generate microscale
variations that ruin the projective extrapolation of
\verb|PIRK2()|.  In such cases, a fixed time-step microscale
integrator is much better (such as \verb|rk2Int()|). 

\item \verb|tSpan| is an \(\ell\)-vector of times at which
the user requests output, of which the first element is
always the initial time. \verb|PIRK2()| does not use
adaptive time-stepping; the macroscale time-steps are
(nearly) the steps between elements of \verb|tSpan|.

\item \verb|x0| is an \(n\)-vector of initial values at the
initial time \verb|tSpan(1)|. Elements of~\verb|x0| may be
\verb|NaN|: such \verb|Nan|s are carried in the simulation
through to the output, and often represent boundaries\slash
edges in spatial fields.

\item \verb|bT|, \emph{optional}, either missing, or
empty~(\verb|[]|), or a scalar: if a given scalar, then it
is the length of the micro-burst simulations---the minimum
amount of time needed for the microscale simulation to relax
to the slow manifold; else if missing or~\verb|[]|, then
\verb|microBurst()| must itself determine the length of a
burst.
\begin{matlab}
%}
if nargin<4, bT=[]; end
%{
\end{matlab}

\end{itemize}


\paragraph{Choose a long enough burst length}
Suppose: firstly, you have some desired relative
accuracy~\(\varepsilon\) that you wish to achieve (e.g.,
\(\varepsilon\approx0.01\) for two digit accuracy);
secondly, the slow dynamics of your system occurs at
rate\slash frequency of magnitude about~\(\alpha\); and
thirdly, the rate of \emph{decay} of your fast modes are
faster than the lower bound~\(\beta\) (e.g., if three fast
modes decay roughly like \(e^{-12t}, e^{-34t}, e^{-56t}\)
then \(\beta\approx 12\)).
\begin{figure}
\centering
\def\aD{\alpha\Delta}\def\bD{\beta\Delta}\def\dD{\delta/\Delta}
\caption{\label{fig:bTlength}Need macroscale step~\(\Delta\)
such that $|\aD|\lesssim\sqrt{6\varepsilon}$ for given
relative error~\(\varepsilon\) and slow rate~\(\alpha\), and
then $\dD\gtrsim\frac1{\bD}\log|\bD|$ determines the minimum
required burst length~\(\delta\) for every given fast
rate~\(\beta\).}
\tikzsetnextfilename{ProjInt/bTlength}
\begin{tikzpicture}
  \begin{loglogaxis}[xlabel={$\bD$}
  ,ylabel={$(\dD)_{\min}$}
  ,domain=2.7:1000 ,grid=both ]
  \addplot+[no marks]{ln(x)/x};
  \end{loglogaxis} 
\end{tikzpicture}
\end{figure}
Then set 
\begin{enumerate}
\item a macroscale time-step, \(\Delta=\verb|diff(tSpan)|\),
such that \(\alpha\Delta\approx\sqrt{6\varepsilon}\), and
\item a microscale burst length, \(\delta=\verb|bT| \gtrsim
\frac1\beta\log|\beta\Delta|\), see \cref{fig:bTlength}.
\end{enumerate}




\paragraph{Output}
If there are no output arguments specified, then a plot is 
drawn of the computed solution~\verb|x| versus \verb|tSpan|.
\begin{itemize}
\item  \verb|x|, an \(\ell \times n\) array of the
approximate solution vector. Each row is an estimated state
at the corresponding time in \verb|tSpan|.  The simplest
usage is then \verb|x = PIRK2(microBurst,tSpan,x0,bT)|.

However, microscale details of the underlying Projective
Integration computations may be helpful. \verb|PIRK2()|
provides up to four optional outputs of the microscale
bursts. 

\item \verb|tms|, optional, is an \(L\) dimensional column
vector containing the microscale times within the burst
simulations, each burst separated by~\verb|NaN|; 

\item \verb|xms|, optional, is an \(L\times n\) array of the
corresponding microscale states---each rows is an accurate
estimate of the state at the corresponding time~\verb|tms|
and helps visualise details of the solution. 

\item \verb|rm|, optional, a struct containing the
`remaining' applications of the microBurst required by the
Projective Integration method during the calculation of the
macrostep: \begin{itemize}
\item \verb|rm.t|~is a column vector of microscale times; and 
\item \verb|rm.x|~is the array of corresponding burst states.
\end{itemize}
The states \verb|rm.x| do not have the same physical
interpretation as those in \verb|xms|; the \verb|rm.x| are
required in order to estimate the slow vector field during
the calculation of the Runge--Kutta increments, and do
\emph{not} accurately approximate the macroscale dynamics.

\item  \verb|svf|, optional, a struct containing the
Projective Integration estimates of the slow vector field.
\begin{itemize}
\item \verb|svf.t| is a \(2\ell\) dimensional column vector
containing all times at which the Projective Integration
scheme is extrapolated along microBurst data to form a
macrostep. 
\item \verb|svf.dx| is a \(2\ell\times n\) array containing
the estimated slow vector field.
\end{itemize}

\end{itemize}







\subsection{If no arguments, then execute an example}
\label{sec:pirk2eg}
\begin{matlab}
%}
if nargin==0
%{
\end{matlab}
\paragraph{Example code for Michaelis--Menton dynamics} The
Michaelis--Menten enzyme kinetics is expressed as a
singularly perturbed system of differential equations for
\(x(t)\) and~\(y(t)\):
\begin{equation*}
\frac{dx}{dt}=-x+(x+\tfrac12)y \quad\text{and}\quad
\frac{dy}{dt}=\frac1\epsilon\big[x-(x+1)y\big]
\end{equation*}
(encoded in function \verb|MMburst()| in the next
paragraph). With initial conditions \(x(0)=1\) and
\(y(0)=0\), the following code computes and plots a solution
over time \(0\leq t\leq6\) for parameter
\(\epsilon=0.05\)\,.  Since the rate of decay is
\(\beta\approx 1/\epsilon\) we choose a burst length
\(\epsilon\log(\Delta/\epsilon)\) as here the macroscale
time-step \(\Delta=1\).
\begin{matlab}
%}
global MMepsilon
MMepsilon = 0.05
ts = 0:6
bT = MMepsilon*log( (ts(2)-ts(1))/MMepsilon )
[x,tms,xms] = PIRK2(@MMburst, ts, [1;0], bT);
figure, plot(ts,x,'o:',tms,xms)
title('Projective integration of Michaelis--Menten enzyme kinetics')
xlabel('time t'), legend('x(t)','y(t)')
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
\input{../ProjInt/odeOct.m}



\begin{devMan}



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
nArgOut=nargout(); 
saveMicro = (nArgOut>1); 
saveFullMicro = (nArgOut>3); 
saveSvf = (nArgOut>4); 
%{
\end{matlab}


Run a preliminary application of the microBurst on the given
initial state to help relax to the slow manifold. This is
done in addition to the microBurst in the main loop, because
the initial state is often far from the attracting slow
manifold. Require the user to input and output rows of the
system state.
\begin{matlab}
%}
x0 = reshape(x0,1,[]); 
[relax_t,relax_x0] = microBurst(tSpan(1),x0,bT);
%{
\end{matlab}

Use the end point of this preliminary microBurst as the
initial state for the loop of macro-steps.
\begin{matlab}
%}
tSpan(1) = relax_t(end); 
x(1,:)=relax_x0(end,:); 
%{
\end{matlab}

If saving information, then record the first application of
the microBurst. Allocate cell arrays for times and states
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
            svf.t = nan(2*nT-2,1); 
            svf.dx = nan(2*nT-2,length(x0)); 
        end
    end
end
%{
\end{matlab}


\paragraph{Loop over the macroscale time-steps}
Also set an initial rounding tolerance for checking.
\begin{matlab}
%}
roundingTol = 1e-8;
for jT = 2:nT
    T = tSpan(jT-1);
%{
\end{matlab}
If two applications of the microBurst would cover one entire
macroscale time-step, then do so (setting some internal
states to \verb|NaN|); else proceed to projective step.
\begin{matlab}
%}
    if ~isempty(bT) && 2*abs(bT)>=abs(tSpan(jT)-T) && bT*(tSpan(jT)-T)>0
        [t1,xm1] = microBurst(T, x(jT-1,:), tSpan(jT)-T);
        x(jT,:) = xm1(end,:);
        t2 = nan;   xm2 = nan(1,size(xm1,2));
        dx1 = xm2;  dx2 = xm2;
    else
%{
\end{matlab}

Run the first application of the microBurst; since this
application directly follows from the initial conditions, or
from the latest macrostep, this microscale information is
physically meaningful as a simulation of the system. Extract
the size of the final time-step. 
\begin{matlab}
%}
    [t1,xm1] = microBurst(T, x(jT-1,:), bT);
%{
\end{matlab}
To estimate the derivative by numerical differentiation, we
balance approximation error~\(\|\ddot x\|/dt\) with
round-off error~\(\|x\|\epsilon/dt\) by the optimal
time-step \(dt\approx\sqrt(\|x\|\epsilon/\|\ddot x\|)\).
Omit~\(\|\ddot x\|\) as we do not know it. Also,
limit~\(dt\) to at most the last tenth of the burst, and at
least one step.
\begin{matlab}
%}
    nt = length(t1);
    optdt = min(0.1*(t1(nt)-t1(1)),sqrt(max(rms(xm1))*1e-15));
    [~,kt] = min(abs(t1(nt)-optdt-t1(1:nt-1)));
    ktnt = [kt nt];
    del = t1(nt)-t1(kt);    
%{
\end{matlab}
Check for round-off error, and decrease tolerance so that
warnings are not repeated unless things get worse.
\begin{matlab}
%}
    xt = [reshape(t1(ktnt),[],1) xm1(ktnt,:)];
    if norm(diff(xt))/norm(xt,'fro') < roundingTol
    warning(['significant round-off error in 1st projection at T=' num2str(T)])
    roundingTol = roundingTol/10; 
    end
%{
\end{matlab}

Find the needed time-step to reach time \verb|tSpan(n+1)|
and form a first estimate \verb|dx1| of the slow vector
field.
\begin{matlab}
%}
    Dt = tSpan(jT)-t1(end);  
    dx1 = (xm1(nt,:)-xm1(kt,:))/del; 
%{
\end{matlab}

Project along \verb|dx1| to form an intermediate
approximation of~\verb|x|; run another application of the
microBurst and form a second estimate of the slow vector
field (assuming the burst length is the same, or nearly so).  
\begin{matlab}
%}
    xint = xm1(end,:) + (Dt-(t1(end)-t1(1)))*dx1;
    [t2,xm2] = microBurst(T+Dt, xint, bT);
%{
\end{matlab}
As before, choose~\(dt\) as best we can to estimate
derivative.
\begin{matlab}
%}
    nt = length(t2);
    optdt = min(0.1*(t2(nt)-t2(1)),sqrt(max(rms(xm2))*1e-15));
    [~,kt] = min(abs(t2(nt)-optdt-t2(1:nt-1)));
    ktnt = [kt nt];
    del = t2(nt)-t2(kt);    
    dx2 = (xm2(nt,:)-xm2(kt,:))/del; 
%{
\end{matlab}
Check for round-off error, and decrease tolerance so that
warnings are not repeated unless things get worse.
\begin{matlab}
%}
    xt = [reshape(t2(ktnt),[],1) xm2(ktnt,:)];
    if norm(diff(xt))/norm(xt,'fro') < roundingTol
    warning(['significant round-off error in 2nd projection at T=' num2str(T)])
    roundingTol = roundingTol/10; 
    end
%{
\end{matlab}

Use the weighted average of the estimates of the slow vector
field to take a macro-step.
\begin{matlab}
%}
    x(jT,:) = xm1(end,:) + Dt*(dx1+dx2)/2; 
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
output of the first application of the microBurst. Separate
bursts by~\verb|NaN|s.
\begin{matlab}
%}
    if saveMicro 
        tms{jT} = [reshape(t1,[],1); nan];
        xms{jT} = [xm1; nan(1,size(xm1,2))];
%{
\end{matlab}

If saving all microscale data, then repeat for the remaining
applications of the microBurst.         
\begin{matlab}
%}
        if saveFullMicro 
            rm.t{jT} = [reshape(t2,[],1); nan];
            rm.x{jT} = [xm2; nan(1,size(xm2,2))];
%{
\end{matlab}

If saving Projective Integration estimates of the slow
vector field, then populate the corresponding cells with
times and estimates.
\begin{matlab}
%}
            if saveSvf 
                svf.t(2*jT-3:2*jT-2) = [t1(end); t2(end)];
                svf.dx(2*jT-3:2*jT-2,:) = [dx1; dx2];
            end
        end
    end
%{
\end{matlab}
End the main loop over all the macro-steps.
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


\subsection{If no output specified, then plot the simulation}
\begin{matlab}
%}
if nArgOut==0
    figure, plot(tSpan,x,'o:')
    title('Projective Simulation with PIRK2')
end
%{
\end{matlab}

This concludes \verb|PIRK2()|.
\begin{matlab}
%}
end
%{
\end{matlab}
\end{devMan}
%}


