%Basic example of PIG. JM, Sept 18.
%!TEX root = ../Doc/equationFreeDoc.tex
%{
\subsection{Example 2: PI using General macrosolvers }
\label{sec:ExPIG}
In this example the PI-General scheme is applied to a singularly perturbed ordinary differential equation. The aim is to allow a standard non-stiff numerical integrator, e.g. \verb|ode45()|, to be applied to a stiff problem on a slow, long time scale.

\begin{body}
\begin{matlab}
%}
clear
%{
\end{matlab}

Set time scale separation and model.
\begin{matlab}
%}
epsilon = 1e-4;
f=@(t,x) [cos(x(1))*sin(x(2))*cos(t); (cos(x(1))-x(2))/epsilon];
%{
\end{matlab}

Set the `black box' microsolver to be an integration using a standard solver, and set the standard time of simulation for the microsolver.
\begin{matlab}
%}
solver = @(tIC, xIC,T) feval('ode45',f,[tIC tIC+T],xIC);
bT=20*epsilon;
%{
\end{matlab}

Set initial conditions, and the time to be covered by the macrosolver. Set the macrosolver to be used as a standard, non-stiff integration scheme.
\begin{matlab}
%}
IC = [1 3];
tspan=[0 40];
macro.tspan = tspan;
macro.solver = 'ode45';
%{
\end{matlab}


Now time and integrate the above system over \verb|tspan| using \verb|PIG()| and, for comparison, a brute force implementation of \verb|ode45()|. Report the time taken by each method.
\begin{matlab}
%}
tic
[t,x,tms,xms] = PIG(solver,bT,macro,IC);
tPI=toc;
fprintf(['PI took %f seconds, using ode45 as the '...
    'macrosolver.\n'],tPI)
tic
[t45,xode45] = ode45(f,[tspan(1) tspan(end)],IC);
tODE45 = toc;

fprintf('Brute force ode45 took %f seconds.\n',tODE45)
%{
\end{matlab}



Plot the output on two figures, showing the truth and macrosteps on both, and all applications of the microsolver on the first figure.
\begin{matlab}
%}
figure; set(gcf,'PaperPosition',[0 0 14 10])
hold on
PIm=plot(tms,xms,'b.');
PI=plot(t,x,'bo');
ODE45=plot(t45,xode45,'r-','LineWidth',2);
legend([PI(1),ODE45(1),PIm(1)],'PI Solution',...
    'Standard Solution','PI microsolver')
xlabel('Time')
ylabel('State')
axis([0 40 0 3])

figure; set(gcf,'PaperPosition',[0 0 14 10])
hold on
PI2=plot(t,x,'bo');
ODE452=plot(t45,xode45,'r-','LineWidth',2);

legend([PI2(1),ODE452(1)],'PI Solution','Standard Solution')
xlabel('Time')
ylabel('State')
%{
\end{matlab}
The output is plotted in Figure~\ref{fig:PIG}.
\begin{figure}
%\includegraphics[width=0.5\textwidth]{../ProjInt/PIG}
\includegraphics[width=0.8\textwidth]{../ProjInt/PIGm}
\caption{Accurate simulation of a stiff nonautonomous system by PIG(). The microsolver is called on-the-fly by the macrosolver (here ode45).}\label{fig:PIG}
\end{figure}


Notes:
\begin{itemize}
\item the problem may be made more, or less, stiff by changing the time scale parameter \verb|epsilon|. \verb|PIG()| will handle a stiffer problem with ease; but if the problem is insufficiently stiff, then the algorithm will produce nonsense. This problem is handled by \verb|cdmc()|; see Section~\ref{sec:Excdmc}.
\item The mildly stiff problem in Example~\ref{sec:ExPIRK} may be efficiently solved by a standard solver, e.g. \verb|ode45()|. The stiff but low dimensional problem in this example can be solved efficiently by a standard stiff solver, e.g. \verb|ode15s()|. The real advantage of the PI schemes is in high dimensional stiff problems, that cannot be efficiently solved by most standard methods.
\end{itemize}

\end{body}
%}