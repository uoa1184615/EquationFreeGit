%Exploration of cdmc. JM, Sept 18.
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{Explore: PI using constraint-defined manifold computing}
\label{sec:Excdmc}
In this example the PI-General scheme is applied to a singularly perturbed ordinary differential equation in which the time scale separation is not too strong. The resulting simulation is not accurate. In parallel, we run the same scheme but with \verb|cdmc()| used as a wrapper for the microsolver. This second implementation successfully replicates the true dynamics.
\begin{matlab}
%}
clear
%{
\end{matlab}

Set a weak time scale separation and model.
\begin{matlab}
%}
epsilon = 1e-2;
f=@(t,x) [cos(x(1))*sin(x(2))*cos(t); (cos(x(1))-x(2))/epsilon];
%{
\end{matlab}

Set the `naive' microsolver to be an integration using a standard solver, and set the standard time of simulation for the microsolver.
\begin{matlab}
%}
naiveSol = @(tIC, xIC,T) feval('ode45',f,[tIC tIC+T],xIC);
bT=5*epsilon;
%{
\end{matlab}
Create a second struct in which the solver is the output of \verb|cdmc()|.
\begin{matlab}
%}
cSol = @(t,x,T) cdmc(naiveSol,t,x,T);
%{
\end{matlab}

Set initial conditions, and the time to be covered by the macrosolver. Set the macrosolver to be used as a standard, non-stiff integration scheme.
\begin{matlab}
%}
IC = [1 3];
tspan=0:0.5:40;
macro.tspan = tspan;
macro.solver = 'ode45';
%{
\end{matlab}


Simulate using \verb|PIG()| with each of the above microsolvers. Generate a trusted solution using standard numerical methods.
\begin{matlab}
%}
[nt,nx] = PIG(naiveSol,bT,macro,IC);
[ct,cx] = PIG(cSol,bT,macro,IC);
[t45,xode45] = ode45(f,[tspan(1) tspan(end)],IC);
%{
\end{matlab}



Plot the output.
\begin{matlab}
%}
figure; set(gcf,'PaperPosition',[0 0 14 10])
hold on
nPI = plot(nt,nx,'bo');
PI=plot(ct,cx,'ko');
ODE45=plot(t45,xode45,'r-','LineWidth',2);

legend([nPI(1),PI(1),ODE45(1)],'Naive PIG Solution',...
    'PIG using cdmc','Accurate Solution')
xlabel('Time')
ylabel('State')
axis([0 40 0 3])
%{
\end{matlab}
The output is plotted in Figure~\ref{fig:PIGE}. The source of the error in the standard \verb|PIG()| scheme is the burst length \verb|bT|, that is significant on the slow time scale. Set \verb|bT| to \verb|20*epsilon| or \verb|50*epsilon|\footnote{this example is quite extreme: at bT=50*epsilon, it would be computationally much cheaper to simulate the entire length of tspan using the microsolver alone.} to worsen the error in both schemes. This example reflects a general principle, that most PI schemes will incur a global error term which is proportional to the simulation time of the microsolver and independent of the order of the microsolver. The \verb|PIRK()| schemes have been written to minimise, if not eliminate entirely, this error, but by design \verb|PIG()| works with any user-defined macrosolver and cannot reduce this error. The function \verb|cdmc()| reduces this error term by attempting to mimic the microsolver without advancing time. 
\begin{figure}
\includegraphics[width=0.8\textwidth]{ProjInt/PIGExplore}
\caption{Accurate simulation of a weakly stiff nonautonomous system by PIG() using cdmc(), and an inaccurate solution using a naive application of PIG().}\label{fig:PIGE}
\end{figure}
%}