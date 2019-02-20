% PIG implements Projective Integration scheme with any
% system or user-specified integrator for the slow-time
% macroscale, and a user-specified microsolver. JM,
% September 2018.
%!TEX root = ../Doc/eqnFreeDevMan.tex

%{
\section{\texttt{PIG()}: Projective Integration via a General macroscale integrator}
\label{sec:PIG}

This is an approximate Projective Integration scheme when
the macroscale integrator is any coded scheme. The advantage
is that one may use \script's inbuilt integration functions,
with all their sophisticated error control and adaptive
time-stepping, to do the macroscale simulation.

Unlike the \verb|PIRKn| functions, \verb|PIG()| does not
estimate the slow vector field at the times expected by any
user-specified scheme, but instead provides an estimate of
the slow vector field at a slightly different time, after an
application of the micro-burst simulator.  Consequently
\verb|PIG()| will incur an additional global error term
proportional to the burst length of the microscale
simulator. For that reason, \verb|PIG()| should be used with
\begin{itemize}
\item either very stiff problems, in which the burst length
of the micro-burst can be short, 
\item or the `constraint defined manifold' based micro-burst
provided by \verb|cdmc()|, that attempts to project the
variables onto the slow manifold without affecting the time.
\end{itemize}


\begin{matlab}
%}
function [t,x,tms,xms,svf] = PIG(macroInt,microBurst,tSpan,x0,lift,restrict)
%{
\end{matlab}

The inputs and outputs are a little different to the two
\verb|PIRKn| functions.
\paragraph{Inputs:}
\begin{itemize}
\item \verb|microBurst()| is a function that produces output
from the user-specified code for a burst of microscale
simulation. The function must know how long a burst it is to
use.  Usage
\begin{equation*}
\verb|[tbs,xbs] = microBurst(tb0,xb0)|
\end{equation*}
\emph{Inputs:} \verb|tb0| is the start time of a burst;
\verb|xb0|~is the vector state at the start of a burst. 
 
\emph{Outputs:} \verb|tbs|, the vector of solution times; and
\verb|xbs|, the corresponding states.

\item \verb|macroInt()|, the numerical method that the user
wants to apply on a slow-time macroscale. Either use a
standard \script\ integration function (such as \verb|ode23|
or~\verb|ode45|), or code this solver as a standard \script\
integration function. That is, if you code you own, then it
must be 
\begin{equation*}
  \verb|[ts,xs] = macroInt(f,tSpan,x0)| 
\end{equation*}  
where function \verb|f(t,x)| notionally evaluates the time
derivatives \(d\xv/dt\) at `any' time; \verb|tSpan|~is
either the macro-time interval, or the vector of times at
which a macroscale value is to be returned; and~\verb|x0|
are the initial values of~\(\xv\) at time \verb|tSpan(1)|.
Then the \(i\)th~\emph{row} of~\verb|xs|, \verb|xs(i,:)|, is
to be the vector~\(\xv(t)\) at time \(t=\verb|ts(i)|\).
Remember that in \verb|PIG()| the function \verb|f(t,x)| is
to be estimated by Projective Integration burst.

\item \verb|tSpan|, a vector of times at which the user
requests output, of which the first element is always the
initial time.  If \verb|macroInt| can adaptively select time
steps (e.g., \verb|ode45|), then \verb|tSpan| can consist of
an initial and final time only.

\item \verb|x0|, the vector of initial values at the initial
time~\verb|tSpan(1)|.
\end{itemize}



\paragraph{Output}
If there are no output arguments specified, then a plot is
drawn of the computed solution~\verb|x| versus~\verb|t|.
Most often you would only store the first two output results
of \verb|PIG()|, via say \verb|[t,x] = PIG(...)|.
\begin{itemize}
\item  \verb|t|, an \(\ell\)-vector of times at which
\verb|macroInt| produced results.
\item \verb|x|, an \(\ell \times n\) array of the computed
solution: the \(i\)th~\emph{row} of~\verb|x|, \verb|x(i,:)|,
is to be the vector~\(\xv(t)\) at time \(t=\verb|t(i)|\).

However, microscale details of the underlying Projective
Integration computations may be helpful, and so \verb|PIG()|
some optional outputs of the microscale bursts. 


\item \verb|tms|, optional, is an \(L\) dimensional column
vector containing microscale times of burst simulations,
each burst separated by~\verb|NaN|; 

\item \verb|xms|, optional, is an \(L\times n\) array of the
corresponding microscale states---this data is an accurate
simulation of the state and may help visualise more details
of the solution. 

\item  \verb|svf|, optional, a struct containing the
Projective Integration estimates of the slow vector field.
\begin{itemize}
\item \verb|svf.t| is a \(2\ell\) dimensional column vector
containing all times at which the Projective Integration
scheme is extrapolated along microsolver data to form a
macrostep. 
\item \verb|svf.dx| is a \(2\ell\times n\) array containing
the estimated slow vector field.
\end{itemize}
\end{itemize}



\subsection{If no arguments, then execute an example}
\label{sec:pigeg}
\begin{matlab}
%}
if nargin==0
%{
\end{matlab}
As a basic example, consider a singularly perturbed system
of differential equations for \(\xv(t)=(x_1(t),x_2(t))\):
\begin{equation*}
\frac{dx_1}{dt}=\cos(x_1)\sin(x_2)\cos(t) \quad\text{and}\quad
\frac{dx_2}{dt}=\frac1\epsilon\big[\cos(x_1)-x_2\big].
\end{equation*}
With initial conditions \(\xv(0)=(1,0)\), the following code
computes and plots a solution of the system over time
\(0\leq t\leq6\) for parameter \(\epsilon=10^{-3}\)\,. 

First we code the right-hand side function of the microscale
system of \ode{}s.
\begin{matlab}
%}
epsilon = 1e-3;
dxdt=@(t,x) [ cos(x(1))*sin(x(2))*cos(t)
             (cos(x(1))-x(2))/epsilon ];
%{
\end{matlab}
Second, we code microscale bursts, here using the standard
\verb|ode45()|. Since the rate of decay is \(\beta\approx
1/\epsilon\) we choose a burst length
\(2\epsilon\log(1/\epsilon)\) as here we do not know the
macroscale time step invoked by \verb|marcoInt()|, so
blithely use \(\Delta=1\), and then double the usual formula
for safety.
\begin{matlab}
%}
bT = 2*epsilon*log(1/epsilon)
microBurst = @(tb0,xb0) ode45(dxdt,[tb0 tb0+bT],xb0);
%{
\end{matlab}
Third, invoke \verb|PIG| to use \verb|ode23()|, say, on the
macroscale slow evolution. Integrate the micro-bursts over
\(0\leq t\leq6\) from initial condition \(\xv=(1,0)\). (You
could set \verb|tSpan=[0 -6]| to integrate backwards in time
with forward bursts.)
\begin{matlab}
%}
tSpan = [0 6]; 
lift = @(x) [x; 0.5];
restrict = @(x) x(1);
[ts,xs,tms,xms] = PIG('ode23',microBurst,tSpan,1, lift, restrict);
%{
\end{matlab}
Plot output of this projective integration.
\begin{matlab}
%}
figure, plot(ts,xs,'o:',tms,xms)
title('Projective integration of singular perturbed ODE')
xlabel('time t'), legend('x_1(t)','x_2(t)')
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

Find the number of time steps at which output is expected,
and the number of variables.
\begin{matlab}
%}
nT=length(tSpan)-1;
nx = length(lift(x0));
%{
\end{matlab}

Get the number of expected outputs and set logical indices
to flag what data should be saved. If no lifting/restriction operators were
set, assign them.
\begin{matlab}
%}
nArgs=nargout();
saveMicro = (nArgs>1);
saveSvf = (nArgs>2);
if nargin < 5 %no lift/restrict operators
    lift=@(x) x;
    restrict=@(x) x;
end
%{
\end{matlab}



Run a first application of the microBurst on the initial
conditions. This is done in addition to the microBurst in
the main loop, because the initial conditions are often far
from the attracting slow manifold.
\begin{matlab}
%}
x0 = reshape(x0,[],1);
[relax_t,x0_micro_relax] = microBurst(tSpan(1),lift(x0));
x0_relax = restrict(x0_micro_relax);
%{
\end{matlab}

Update the initial time.
\begin{matlab}
%}
tSpan(1) = relax_t(end); 
%{
\end{matlab}
Allocate cell arrays for times and states for any of the
outputs requested by the user. If saving information, then
record the first application of the microBurst. Note that it
is unknown a priori how many applications of the microBurst
will be required; this code may be run more efficiently if
the correct number is used in place of \verb|nT+1| as the
dimension of the cell arrays.
\begin{matlab}
%}
if saveMicro
    tms=cell(nT+1,1); xms=cell(nT+1,1);
    n=1;
    tms{n} = reshape(relax_t,[],1);
    xms{n} = x0_micro_relax;
    
    if saveSvf
        svf.t = cell(nT+1,1);
        svf.dx = cell(nT+1,1);
    end
end
%{
\end{matlab}

The idea of \verb|PIG()| is to use the output from the
microBurst to approximate an unknown function
\verb|ff(t,x)|, that describes the slow dynamics. This
approximation is then used in the system\slash user-defined
`coarse solver' \verb|macroInt()|. The approximation is
described in
\begin{matlab}
%}
function [dx]=genProjection(tt,xx)
%{
\end{matlab}
Run a microBurst from the given initial conditions.
\begin{matlab}
%}
    [t_tmp,x_micro_tmp] = microBurst(tt,reshape(lift(xx),[],1));
%{
\end{matlab}
Compute the standard Projective Integration approximation of
the slow vector field.
\begin{matlab}
%}
	del = t_tmp(end)-t_tmp(end-1);
	dx = ( restrict(x_micro_tmp(end,:))-restrict(x_micro_tmp(end-1,:)) )'/(del);
%{
\end{matlab}
Save the microscale data, and the Projective Integration
slow vector field, if requested.
\begin{matlab}
%}
	if saveMicro
        n=n+1;
        tms{n} = [reshape(t_tmp,[],1); nan];
        xms{n} = [x_micro_tmp; nan(1,nx)];
        if saveSvf
            svf.t{n-1} = tt;
            svf.dx{n-1} = dx;
        end
	end
end% function genProjection()
%{
\end{matlab}


Define the approximate slow vector field according to
Projective Integration.
\begin{matlab}
%}
ff=@(t,x) genProjection(t,x);
%{
\end{matlab}


Do Projective Integration of \verb|ff()| with the
user-specified microBurst.
\begin{matlab}
%}
[t,x]=feval(macroInt,ff,tSpan,x0_relax(end,:)');
%{
\end{matlab}

Overwrite \verb|x(1,:)| and \verb|t(1)|, which the user
expect to be \verb|x0| and \verb|tSpan(1)| respectively,
with the given initial conditions.
\begin{matlab}
%}
x(1,:) = x0';
t(1) = tSpan(1);
%{
\end{matlab}

For each additional requested output, concatenate all the
cells of time and state data into two arrays. Then, return
the two arrays in a cell.
\begin{matlab}
%}
if saveMicro
    tms = cell2mat(tms);
    xms = cell2mat(xms);
    if saveSvf
        svf.t = cell2mat(svf.t);
        svf.dx = cell2mat(svf.dx);
    end
end
%{
\end{matlab}


\subsection{If no output specified, then plot simulation}
\begin{matlab}
%}
if nArgs==0
    fifure, plot(t,x,'o:')
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
