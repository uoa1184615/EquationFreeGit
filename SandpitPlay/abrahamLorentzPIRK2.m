% Abraham--Lorentz example of forward integration by bursts
% of backwards simulation.  AJR Oct 2023
%{
\section{\texttt{abrahamLorentz}:  projective integration of
Abraham--Lorentz system}
\label{sec:abrahamLorentz}
%\localtableofcontents

Investigate the Abraham--Lorentz system (4)--(6) of
\cite{Burby2020}. In terms of an electron's (\(\zeta=-1\))
non-dimensional position~\xv, velocity~\vv, and
acceleration~\av, the non-dimensional \ode{}s are 
\begin{subequations}\label{eqsAL}%
\begin{align}&
\dot\xv=\vv,  \label{eqx}
\\&
\epsilon_B\dot\vv=\av,  \label{eqv}
\\&
\epsilon_R\dot\av=\av+\vv\times\Bv(\xv),  \label{eqa}
\end{align}
\end{subequations}
for two possibly-small parameters \(\epsilon_R:=2r_0/(3cT)\)
(slight change for brevity so I do not have various
\(3/2\)~factors), and \(\epsilon_B:=1/(|\omega_c|T)\). When
convenient, I denote the components of the magnetic
field~\(\Bv(\xv)\) as~\((A,B,C)\).


For \(\epsilon_R\ll\epsilon_B\ll1\) the system has six slow
modes and three unstable modes (growth
rate~\(1/\epsilon_R\)).


\paragraph{Invoke projective integration}

Clear, and set the scale separation parameters.
Define some magnetic field function of position.
\begin{matlab}
%}
clear all, close all
global epsB epsR B
epsB=0.2
epsR=0.005
B=@(x) [-x(2)/2; x(1)/2; 1]
%{
\end{matlab}



\paragraph{Integrate forward in time}
Projective integration even simulates forward in time along
the slow manifold using short backward bursts
\cite[]{Gear03b, Frewen2009}. Such forward macroscale
simulations succeed despite the fast modes being viciously
unstable.  
\begin{matlab}
%}
ts = 0:epsB/2:1;
[xs,tMicro,xMicro] = PIRK2(@ALburst, ts, rand(9,1), 4*epsR);
figure
subplot(2,1,1)
plot(ts,xs(:,1:6),'o:',tMicro,xMicro(:,1:6))
if max(abs(xs(:)))>9, quasiLogAxes(Inf,1), end
xlabel('time t'), legend('x','y','z','u','v','w')
title('forward integration points, with backwards bursts')
subplot(2,1,2)
plot(ts,xs(:,7:9),'o:',tMicro,xMicro(:,7:9))
if max(abs(xs(:)))>9, quasiLogAxes(Inf,1), end
xlabel('time t'), legend('a1','a2','a3')
%ifOurCf2eps([mfilename '3'])
exportgraphics(gcf,mfilename+".pdf",'contenttype','vector')
%{
\end{matlab}




\paragraph{Code a burst of Abraham--Lorentz system}
Integrate a burst of length~\verb|bT| of the \ode{}s for the
Abraham--Lorentz system. Code \textsc{ode}s in
function~\verb|dALdt|   Starting at time~\verb|ti|, and
state~\verb|xi| (row), we here simply use \script's
\verb|ode23| to integrate a burst in time.
\begin{matlab}
%}
function [ts, xs] = ALburst(ti, xi, bT) 
global epsB epsR B
    dALdt = @(t,x) [ x(4:6)
                     x(7:9)/epsB
                ( x(7:9)+cross(x(4:6),B(x(1:3))) )/epsR ];
    [ts, xs] = ode23(dALdt, [ti ti-bT], xi);
end
%{
\end{matlab}

%}
