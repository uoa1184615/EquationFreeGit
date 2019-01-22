%Michaelis--Menton example of projective integrating
%fast-slow system.  This example simply introduces basic
%usage of the PIRK2() function. AJR, 29 Sep 2018
%!TEX root = ../Doc/equationFreeDoc.tex
%{
\begin{userExample}
\subsection{\texttt{egPIMM}: Example projective integration
of Michaelis--Menton kinetics}
\label{sec:egPIMM}
\localtableofcontents

The Michaelis--Menten enzyme kinetics is expressed as a
singularly perturbed system of differential equations for
\(x(t)\) and~\(y(t)\):
\begin{equation*}
\frac{dx}{dt}=-x+(x+\tfrac12)y \quad\text{and}\quad
\frac{dy}{dt}=\frac1\epsilon\big[x-(x+1)y\big].
\end{equation*}
As illustrated in \cref{fig:egPIMM2}, the slow
variable~\(x(t)\) evolves on a time scale of one, whereas
the fast variable~\(y(t)\) evolves on a time scale of the
small parameter~\(\epsilon\).

\subsubsection{Invoke projective integration}

Clear, and set the scale separation parameter~\(\epsilon\)
to something small like~\(0.01\). Here use \(\epsilon=0.1\)
for clearer graphs.
\begin{matlab}
%}
clear all, close all
global epsilon
epsilon = 0.1
%{
\end{matlab}

First, \cref{sec:egPIMMburst} encodes the computation of
bursts of the Michaelis--Menten system in a
function~\verb|MMburst()|. Second, here set macroscale times
of computation and interest into vector~\verb|ts|. Then,
invoke Projective Integration with \verb|PIRK2()| applied to
the burst function, say using bursts of simulations of
length~\(2\epsilon\), and starting from the initial
condition for the Michaelis--Menten system of
\((x(0),y(0))=(1,0)\) (off the slow manifold).
\begin{matlab}
%}
ts = 0:6
xs = PIRK2(@MMburst, 2*epsilon, ts, [1;0])
plot(ts,xs,'o:')
xlabel('time t'), legend('x(t)','y(t)')
pause(1)
%{
\end{matlab}
\cref{fig:egPIMM1} plots the macroscale results showing the
long time decay of the Michaelis--Menten system on the slow
manifold.
\begin{figure}
\centering
\caption{\label{fig:egPIMM1}Michaelis--Menten enzyme
kinetics simulated with the projective integration of
\texttt{PIRK2()}: macroscale samples.}
\includegraphics[scale=0.85]{../ProjInt/egPIMM1}
\end{figure}%
\cite{Sieber2018} [\S4] used this system as an example of 
their analysis of the convergence of Projective Integration.


\paragraph{Optional: request and plot the microscale bursts}
Because the initial conditions of the simulation are off the
slow manifold, the initial macroscale step appears to `jump'
(\cref{fig:egPIMM1}). To see the initial transient
attraction to the slow manifold we plot some microscale data
in \cref{fig:egPIMM2}. Two further output variables provide
this microscale burst information.
\begin{matlab}
%}
[xs,tMicro,xMicro] = PIRK2(@MMburst, 2*epsilon, ts, [1;0]);
figure, plot(ts,xs,'o:',tMicro,xMicro)
xlabel('time t'), legend('x(t)','y(t)')
pause(1)
%{
\end{matlab}
\cref{fig:egPIMM2} plots the macroscale and microscale
results---also showing that the initial burst is by default
twice as long. Observe the slow variable~\(x(t)\) is also
affected by the initial transient which indicates that other
schemes which `freeze' slow variables are less accurate.
\begin{figure}
\centering
\caption{\label{fig:egPIMM2}Michaelis--Menten enzyme
kinetics simulated with the projective integration of
\texttt{PIRK2()}: the microscale bursts show the initial
transients on a time scale of \(\epsilon=0.1\), and then the
alignment along the slow manifold.}
\includegraphics[scale=0.85]{../ProjInt/egPIMM2}
\end{figure}



\paragraph{Optional: simulate backwards in time}
\cref{fig:egPIMM3} shows that projective integration even
simulates backwards in time along the slow manifold using
short forward bursts. Such backwards macroscale simulations
succeed despite the fast variable~\(y(t)\), when backwards
in time, being viciously unstable. However, backwards
integration appears to need longer bursts,
here~\(3\epsilon\).
\begin{matlab}
%}
ts = 0:-1:-5
[xs,tMicro,xMicro] = PIRK2(@MMburst, 3*epsilon, ts, 0.2*[1;1]);
figure, plot(ts,xs,'o:',tMicro,xMicro)
xlabel('time t'), legend('x(t)','y(t)')
%{
\end{matlab}
\begin{figure}
\centering
\caption{\label{fig:egPIMM3}Michaelis--Menten enzyme
kinetics simulated backwards with the projective integration
of \texttt{PIRK2()}: the microscale bursts show the short
forward simulations used to project backwards in time at
\(\epsilon=0.1\).}
\includegraphics[scale=0.85]{../ProjInt/egPIMM3}
\end{figure}



\subsubsection{Code a burst of Michaelis--Menten enzyme kinetics}
\label{sec:egPIMMburst}
Say use \verb|ode23()| to integrate a burst of the
differential equations for the Michaelis--Menten enzyme
kinetics. Code differential equations in
function~\verb|dMMdt| with variables \(x=\verb|x(1)|\) and
\(y=\verb|x(2)|\). For the given start time~\verb|ti|, and
start state~\verb|xi|, \verb|ode23()| integrates the
differential equations for a burst time of~\verb|bT|, and
return the simulation data.
\begin{matlab}
%}
function [ts, xs] = MMburst(ti, xi, bT) 
    global epsilon
    dMMdt = @(t,x) [ -x(1)+(x(1)+0.5)*x(2)
          1/epsilon*( x(1)-(x(1)+1)*x(2) ) ];
    [ts, xs] = ode23(dMMdt, [ti ti+bT], xi);
end
%{
\end{matlab}
\end{userExample}
%}
