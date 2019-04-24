% PIG implements Projective Integration scheme with any
% inbuilt integrator or user-specified integrator for the
% slow-time macroscale, and with any inbuilt/user-specified
% microsolver. JM & AJR, Sept 2018 -- Apr 2019.
%!TEX root = ../Doc/eqnFreeDevMan.tex

%{
\section{\texttt{PIG()}: Projective Integration via a General macroscale integrator}
\label{sec:PIG}
\localtableofcontents

\subsection{Introduction}

This is a Projective Integration scheme when the macroscale
integrator is any specified coded method. The advantage is
that one may use \script's inbuilt integration functions,
with all their sophisticated error control and adaptive
time-stepping, to do the macroscale integration\slash
simulation.

By default, for the microscale simulations \verb|PIG()| uses `constraint-defined manifold
computing', \verb|cdmc()| (\cref{sec:cdmc}).  This algorithm,
initiated by \cite{Gear05}, uses a backward projection so
that the simulation time is unchanged after running the
microscale simulator. 

\begin{matlab}
%}
function [T,X,tms,xms,svf] = PIG(macroInt,microBurst,Tspan,x0 ...
                                ,restrict,lift,cdmcFlag)
%{
\end{matlab}

\paragraph{Inputs:}
\begin{itemize}

\item \verb|macroInt()|, the numerical method that the user
wants to apply on a slow-time macroscale. Either specify a
standard \script\ integration function (such as
\verb|'ode23'| or~\verb|'ode45|'), or code your own
integration function using standard arguments. That is, if
you code your own, then it must be 
\begin{equation*}
  \verb|[Ts,Xs] = macroInt(F,Tspan,X0)| 
\end{equation*}  
where \begin{itemize}
\item function \verb|F(T,X)| notionally evaluates the time
derivatives \(d\Xv/dt\) at any time; 
\item \verb|Tspan| is either the macro-time interval, or the
vector of macroscale times at which macroscale values are to
be returned; and
\item \verb|X0| are the initial values of~\(\Xv\) at time
\verb|Tspan(1)|.
\end{itemize}

Then the \(i\)th~\emph{row} of~\verb|Xs|, \verb|Xs(i,:)|, is
to be the vector~\(\Xv(t)\) at time \(t=\verb|Ts(i)|\).
Remember that in \verb|PIG()| the function \verb|F(T,X)| is
to be estimated by Projective Integration.

\item \verb|microBurst()| is a function that produces output
from the user-specified code for a burst of microscale
simulation. The function must internally specify\slash decide how long a
burst it is to use.  Usage
\begin{equation*}
\verb|[tbs,xbs] = microBurst(tb0,xb0)|
\end{equation*}
\emph{Inputs:} \verb|tb0| is the start time of a burst;
\verb|xb0|~is the \(n\)-vector microscale state at the start
of a burst. 
 
\emph{Outputs:} \verb|tbs|, the vector of solution times;
and \verb|xbs|, the corresponding microscale states.


\item \verb|Tspan|, a vector of macroscale times at which
the user requests output.  The first element is always the
initial time.  If \verb|macroInt| reports adaptively selected time
steps (e.g., \verb|ode45|), then \verb|Tspan| consists of an
initial and final time only.

\item \verb|x0|, the \(n\)-vector of initial microscale
values at the initial time~\verb|Tspan(1)|.
\end{itemize}

\paragraph{Optional Inputs:}
\verb|PIG()| allows for none, two or three additional inputs
after~\verb|x0|. If you distinguish distinct microscale and
macroscale states and your aim is to do Projective
Integration on the macroscale only, then lifting and
restriction functions must be provided to convert between
them. Usage \verb|PIG(...,restrict,lift)|:
\begin{itemize}
\item \verb|restrict(x)|, a function that takes an input
high-dimensional, \(n\)-D, microscale state~\xv\ and computes the
corresponding low-dimensional, \(N\)-D, macroscale state~\Xv; 
\item \verb|lift(X,xApprox)|, a function that converts an
input low-dimensional, \(N\)-D, macroscale state~\Xv\ to a
corresponding high-dimensional, \(n\)-D, microscale state~\xv, given
that \verb|xApprox| is a recently computed microscale state
on the slow manifold.
\end{itemize}
Either both \verb|restrict()| and \verb|lift()| are to be
defined, or neither. If neither are defined, then they are
assumed to be identity functions, so that \verb|N=n| in the
following. 

If desired, the default constraint-defined manifold
computing microsolver may be disabled, via
\verb|PIG(...,restrict,lift,cdmcFlag)|
\begin{itemize}
\item \verb|cdmcFlag|, \emph{any} seventh input to
\verb|PIG()|, will disable \verb|cdmc()|, e.g., the string
\verb|'cdmc off'|. 
\end{itemize}

If the \verb|cdmcFlag| is to be set without using a
\verb|restrict()| or \verb|lift()| function, then use empty
matrices~\verb|[]| for the restrict and lift functions.


\paragraph{Output}
Between zero and five outputs may be requested. If there are
no output arguments specified, then a plot is drawn of the
computed solution~\verb|X| versus~\verb|T|. Most often you
would store the first two output results of \verb|PIG()|,
via say \verb|[T,X] = PIG(...)|.
\begin{itemize}
\item  \verb|T|, an \(L\)-vector of times at which
\verb|macroInt| produced results.
\item \verb|X|, an \(L \times N\) array of the computed
solution: the \(i\)th~\emph{row} of~\verb|X|, \verb|X(i,:)|,
is to be the macro-state vector~\(\Xv(t)\) at time
\(t=\verb|T(i)|\).
\end{itemize}

However, microscale details of the underlying Projective
Integration computations may be helpful, and so \verb|PIG()|
provides some optional outputs of the microscale bursts, via
\verb|[T,X,tms,xms] = PIG(...)|

\begin{itemize}
\item \verb|tms|, optional, is an \(\ell\)-dimensional column
vector containing microscale times with bursts,
each burst separated by~\verb|NaN|; 

\item \verb|xms|, optional, is an \(\ell\times n\) array of
the corresponding microscale states. 
\end{itemize}
In some contexts it may be helpful to see directly how
Projective Integration approximates a reduced slow vector
field, via \verb|[T,X,tms,xms,svf] = PIG(...)| in which
\begin{itemize}
\item  \verb|svf|, optional, a struct containing the
Projective Integration estimates of the slow vector field.
\begin{itemize}
\item \verb|svf.T| is a \(\hat L\)-dimensional  column
vector containing all times at which the microscale
simulation data is extrapolated to form an estimate of
\(d\xv/dt\) in \verb|macroInt()|.
\item \verb|svf.dX| is a \(\hat L\times N\) array containing
the estimated slow vector field.
\end{itemize}
\end{itemize}
If \verb|macroInt()| is, for example, the forward Euler
method (or the Runge--Kutta method), then \(\hat L = L\) (or
\(\hat L = 4L \)). 





\subsection{If no arguments, then execute an example}
\label{sec:pigeg}

\begin{matlab}
%}
if nargin==0
%{
\end{matlab}
\begin{figure}
\centering
\caption{\label{fig:PIGsing}Projective Integration by
\texttt{PIG} of the example system~\eqref{eq:PIGsing} with \(\epsilon=10^{-3}\) (\cref{sec:pigeg}).  The macroscale solution~\(X(t)\) is
represented by just the blue circles.  The microscale bursts
are the microscale states \((x_1(t),x_2(t)) =
(\text{red},\text{yellow})\) dots.}
\includegraphics[scale=0.9]{PIGsing}
\end{figure}%
As a basic example, consider a microscale system of the
singularly perturbed system of differential equations
\begin{equation}
\frac{dx_1}{dt}=\cos(x_1)\sin(x_2)\cos(t) \quad\text{and}\quad
\frac{dx_2}{dt}=\frac1\epsilon\big[\cos(x_1)-x_2\big].
\label{eq:PIGsing}
\end{equation}
The macroscale variable is \(X(t)=x_1(t)\), and the
evolution \(dX/dt\) is unclear. With initial condition
\(X(0)=1\), the following code computes and plots a solution
of the system~\eqref{eq:PIGsing} over time \(0\leq t\leq6\)
for parameter \(\epsilon=10^{-3}\)(\cref{fig:PIGsing}).
Whenever needed by \verb|microBurst()|, the microscale
system~\eqref{eq:PIGsing} is initialised (`\verb|lift|ed')
using \(x_2(t) = x_2^{\text{approx}}\) (yellow dots in
\cref{fig:PIGsing}).

First we code the right-hand side function of the microscale
system~\eqref{eq:PIGsing} of \ode{}s.
\begin{matlab}
%}
epsilon = 1e-3;
dxdt=@(t,x) [ cos(x(1))*sin(x(2))*cos(t)
             ( cos(x(1))-x(2) )/epsilon ];
%{
\end{matlab}
Second, we code microscale bursts, here using the standard
\verb|ode45()|. We choose a burst length
\(2\epsilon\log(1/\epsilon)\) as the rate of decay is
\(\beta\approx 1/\epsilon\) and we do not know the
macroscale time-step invoked by \verb|macroInt()|, so
blithely assume \(\Delta\le1\) and then double the usual
formula for safety.
\begin{matlab}
%}
bT = 2*epsilon*log(1/epsilon)
if ~exist('OCTAVE_VERSION','builtin')
    micB='ode45'; else micB='rk2Int'; end
microBurst = @(tb0, xb0) feval(micB,dxdt,[tb0 tb0+bT],xb0);
%{
\end{matlab}
Third, code functions to convert between macroscale and
microscale states.
\begin{matlab}
%}
restrict = @(x) x(1);
lift = @(X,xApprox) [X; xApprox(2)];
%{
\end{matlab}

Fourth, invoke \verb|PIG| to use \script's \verb|ode23|\slash\verb|lsode|, say, on the
macroscale slow evolution. Integrate the micro-bursts over
\(0\leq t\leq6\) from initial condition \(\xv(0)=(1,0)\).
You could set \verb|Tspan=[0 -6]| to integrate backward in
macroscale time with forward microscale bursts \cite[]{Gear03b}.
\begin{matlab}
%}
Tspan = [0 6]; 
x0 = [1;0];
if ~exist('OCTAVE_VERSION','builtin')
    macInt='ode23'; else macInt='odeOct'; end
[Ts,Xs,tms,xms] = PIG(macInt,microBurst,Tspan,x0,restrict,lift);
%{
\end{matlab}
Plot output of this projective integration.
\begin{matlab}
%}
figure, plot(Ts,Xs,'o:',tms,xms,'.')
title('Projective integration of singularly perturbed ODE')
xlabel('time t')
legend('X(t) = x_1(t)','x_1(t) micro bursts','x_2(t) micro bursts')
%{
\end{matlab}
Upon finishing execution of the example, exit this function.
\begin{matlab}
%}
return
end%if no arguments
%{
\end{matlab}






\begin{devMan}
\subsection{The projective integration code}


If no lifting/restriction functions are provided, then
assign them to be the identity functions.
\begin{matlab}
%}
if nargin < 5 || isempty(restrict)
    lift=@(X,xApprox) X;
    restrict=@(x) x;
end
%{
\end{matlab}

Get the number of expected outputs and set logical indices
to flag what data should be saved. 
\begin{matlab}
%}
nArgOut = nargout();
saveMicro = (nArgOut>2);
saveSvf = (nArgOut>4);
%{
\end{matlab}

Find the number of time-steps at which output is expected,
and the number of variables.
\begin{matlab}
%}
nT = length(Tspan)-1;
nx = length(x0);
nX = length(restrict(x0));
%{
\end{matlab}

Reformulate the microsolver to use \verb|cdmc()|, unless
flagged otherwise. The result is that the solution from
microBurst will terminate at the given initial time. 
% Should be OK in Octave.
\begin{matlab}
%}
if nargin<7
    microBurst = @(t,x) cdmc(microBurst,t,x);
else
warning(['A ' class(cdmcFlag) ' seventh input to PIG().'...
 ' PIG will not use constraint-defined manifold computing.'])
end
%{
\end{matlab}



Execute a preliminary application of the microBurst on the initial
state. This is done in addition to the microBurst in
the main loop, because the initial state is often far
from the attracting slow manifold.
\begin{matlab}
%}
[relaxT,x0MicroRelax] = microBurst(Tspan(1),x0);
xMicroLast = x0MicroRelax(end,:).';
X0Relax = restrict(xMicroLast);
%{
\end{matlab}

Update the initial time.
\begin{matlab}
%}
Tspan(1) = relaxT(end); 
%{
\end{matlab}
Allocate cell arrays for times and states for any of the
outputs requested by the user. If saving information, then
record the first application of the microBurst. It
is unknown a priori how many applications of microBurst
will be required; this code may be run more efficiently if
the correct number is used in place of \verb|nT+1| as the
dimension of the cell arrays.
\begin{matlab}
%}
if saveMicro
    tms=cell(nT+1,1); xms=cell(nT+1,1);
    n=1;
    tms{n} = reshape(relaxT,[],1);
    xms{n} = x0MicroRelax;
    
    if saveSvf
        svf.T = cell(nT+1,1);
        svf.dX = cell(nT+1,1);
    else
        svf = [];
    end
else
    tms = []; xms = []; svf = [];
end
%{
\end{matlab}

\paragraph{Define a function of macro simulation}
The idea of \verb|PIG()| is to use the output from the
\verb|microBurst()| to approximate an unknown function
\verb|F(t,X)| that computes \(d\Xv/dt\). This approximation
is then used in the system\slash user-defined `coarse
solver' \verb|macroInt()|. The approximation is computed in 
the function
\begin{matlab}
%}
function [dXdt]=PIFun(t,X)
%{
\end{matlab}
Run a microBurst from the given macroscale initial values.
\begin{matlab}
%}
  x = lift(X,xMicroLast);
  [tTmp,xMicroTmp] = microBurst(t,reshape(x,[],1));
  xMicroLast = xMicroTmp(end,:).';
%{
\end{matlab}
Compute the standard Projective Integration approximation of
the slow vector field.
\begin{matlab}
%}
  X2 = restrict(xMicroTmp(end,:));
  X1 = restrict(xMicroTmp(end-1,:));
  dt = tTmp(end)-tTmp(end-1);
  dXdt = (X2 - X1).'/dt;
%{
\end{matlab}
Save the microscale data, and the Projective Integration
slow vector field, if requested.
\begin{matlab}
%}
  if saveMicro
      n=n+1;
      tms{n} = [reshape(tTmp,[],1); nan];
      xms{n} = [xMicroTmp; nan(1,nx)];
      if saveSvf
          svf.T{n-1} = t;
          svf.dX{n-1} = dXdt;
      end
  end
end% PIFun function
%{
\end{matlab}


\paragraph{Invoke the macroscale integration}
Integrate \verb|PIF()| with the user-specified simulator
\verb|macroInt()|. For some reason, in \script\ we need to
use a one-line function, \verb|PIF|, that invokes the above
macroscale function, \verb|PIFun|.  We also need to use
\verb|feval| because \verb|macroInt()| has multiple outputs.
\begin{matlab}
%}
PIF = @(t,x) PIFun(t,x);
[T,X] = feval(macroInt,PIF,Tspan,X0Relax.');
%{
\end{matlab}

Overwrite \verb|X(1,:)| and \verb|T(1)|, which a user
expects to be \verb|X0| and \verb|Tspan(1)| respectively,
with the given initial conditions.
\begin{matlab}
%}
X(1,:) = restrict(x0);
T(1) = Tspan(1);
%{
\end{matlab}

Concatenate all the additional requested outputs into arrays.
\begin{matlab}
%}
if saveMicro
    tms = cell2mat(tms);
    xms = cell2mat(xms);
    if saveSvf
        svf.T = cell2mat(svf.T);
        svf.dX = cell2mat(svf.dX);
    end
end
%{
\end{matlab}


\subsection{If no output specified, then plot the simulation}
\begin{matlab}
%}
if nArgOut==0
    figure, plot(T,X,'o:')
    title('Projective Simulation via PIG')
end
%{
\end{matlab}


This concludes \verb|PIG()|.
\begin{matlab}
%}
end
%{
\end{matlab}
\end{devMan}
%}
