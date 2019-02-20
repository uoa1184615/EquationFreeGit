%Linear example of PIRK. JM, Sept 18.
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{Example: PI using Runge--Kutta macrosolvers }
\label{sec:ExPIRK}
This script is a demonstration of the \verb|PIRK()| schemes,
that use a Runge--Kutta macrosolver, applied to a simple
linear system with some slow and fast directions. 

\begin{devMan}
Clear workspace and set a seed.
\begin{matlab}
%}
clear
rng(1)
%{
\end{matlab}

The majority of this example involves setting up details for
the microsolver. We use a simple function
\verb|gen_linear_system()| that outputs a function \({f(t,x)
= \mathbf{A}\vec x+\vec b}\), where \(\mathbf{A}\) has some
eigenvalues with large negative real part, corresponding to
fast variables and some eigenvalues with real part close to
zero, corresponding to slow variables. The function
\verb|gen_linear_system()| requires that we specify bounds
on the real part of the strongly stable eigenvalues,
\begin{matlab}
%}
fastband = [-5e2; -1e2]; 
%{
\end{matlab}
and bounds on the real part of the weakly stable/unstable
eigenvalues,
\begin{matlab}
%}
slowband = [-0.002; 0.002];
%{
\end{matlab}
We now generate a random linear system with seven fast and
three slow variables.
\begin{matlab}
%}
f = gen_linear_system(7,3,fastband,slowband); 
%{
\end{matlab}
Set the time step size and total integration time of the
microsolver.
\begin{matlab}
%}
dt = 0.001;
bT = 0.05;
%{
\end{matlab}
As a rule of thumb, the time steps \verb|dt| should satisfy
$\verb|dt|  \le 1/|\verb|fastband|(1)|$ and the time to
simulate with each application of the microsolver,
\verb|micro.bT|, should be larger than or equal to
$1/|\verb|fastband|(2)|$. We set the integration scheme to
be used in the microsolver. Since the time steps are so
small, we just use the forward Euler scheme
\begin{matlab}
%}
solver='fe'; 
%{
\end{matlab}
(Other options: \verb|'rk2'| for second order Runge--Kutta,
\verb|'rk4'| for fourth order, or any Matlab/Octave
integrator such as \verb|'ode45'|.)



A crucial part of the \textsc{pi} philosophy is that it does
not assume anything about the microsolver. For this reason,
the microsolver must be a `black box', which is run by
specifying an initial time and state, and a duration to
simulate for. All the details of the microsolver must be set
by the user. We generate and save a black box microsolver.
\begin{matlab}
%}
bbm = bbgen(solver,f,dt);
solver = bbm;
%{
\end{matlab}

Set the macroscale times at which we request output from the
\textsc{pi} scheme and the initial conditions.
\begin{matlab}
%}
tSpan=0: 1 : 30; 
IC = linspace(-10,10,10); 
%{
\end{matlab}



We implement the \textsc{pi} scheme, saving the coarse
states in \verb|x|, the `trusted' applications of the
microsolver in \verb|xmicro|, and the additional
applications of the microsolver in \verb|xrmicro|. Note that
the second and third outputs are optional and do not need to
be set.
\begin{matlab}
%}
[x, tms, xms, rm] = PIRK4(solver, bT, tSpan, IC); 
%{
\end{matlab}
For verification, we also compute the trajectories using a
standard solver.
\begin{matlab}
%}
[tt,ode45x] = ode45(f,tSpan([1,end]),IC);
%{
\end{matlab}

\begin{figure}
\caption{Demonstration of PIRK4(). From initial conditions,
the system rapidly trannsitions to an attracting invariant
manifold. The \textsc{pi} solution accurately tracks the
evolution of the variables over time while requiring only a
fraction of the computations of the standard
solver.}\label{fig:PIRK}
\includegraphics[scale=0.9]{../ProjInt/PIRK}
\end{figure}
\cref{fig:PIRK} plots the output.
\begin{matlab}
%}
tmsr = rm.t; xmsr = rm.x;
clf()
hold on
PI_sol=plot(tSpan,x,'bo');
std_sol=plot(tt,ode45x,'r');
plot(tms,xms,'k.');
plot(tmsr,xmsr,'g.');
legend([PI_sol(1),std_sol(1)],'PI Solution',...
    'Standard Solution','Location','NorthWest')
xlabel('Time');
ylabel('State');
%{
\end{matlab}
Save plot to a file.
\begin{matlab}
%}
set(gcf,'PaperPosition',[0 0 14 10])
print('-depsc2','PIRK')
%{
\end{matlab}
\end{devMan}
%}
