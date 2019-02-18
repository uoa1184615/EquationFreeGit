%Generate a `black box' microsolver suitable for PI from a standard 
%numerical method, an ordinary differential equation, and a given upper
%bound on the time step.
%JM,Sept 2018.
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsubsection{\texttt{bbgen()}}
\label{sec:bbgen}
\verb|bbgen()| is a simple function that takes a standard numerical method and produces a black box solver of the type required by the PI schemes. 
\begin{matlab}
%}
function bb = bbgen(solver,f,dt)
%{
\end{matlab}
\paragraph{Input}
\begin{itemize}
\item \verb|solver|, a standard numerical solver for ordinary differential equations
\item \verb|f|, a function f(t,x) taking time and state inputs
\item \verb|dt|, a time step suitable for simulation with \verb|f|
\end{itemize}
\paragraph{Output}
\verb|bb = bb|\((t_{in},x_{in},T)\) a `black box' microsolver that initialises at \( (t_{in},x_{in}) \) and simulates forward a duration \(T\). 

\begin{funDescription}
\begin{matlab}
%}
bb = @(t_in,x_in,T) feval(solver,f,...
linspace(t_in,t_in+T,1+ceil(T/dt)),x_in);
end
%{
\end{matlab}
\end{funDescription}
%}