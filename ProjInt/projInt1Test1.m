%{
% Script to test the time integration function projInt1() 
% on discretisation of nonlinear diffusion PDE.  
% AJR, Oct 2017
%!TEX root = ../equationFreeDoc.tex
\subsection{A first test of basic projective integration}
\label{sec:ftbpi}

Seek to simulate the nonlinear diffusion \pde\
\begin{equation*}
\D tu=u\DD xu\quad \text{such that }u(\pm 1)=0,
\end{equation*}
with random positive initial condition.
\autoref{fig:pit1u} shows solutions are attracted to the parabolic \(u=a(t)(1-x^2)\) with slow algebraic decay \(\dot a=-2a^2\).
\begin{figure}
\centering
\caption{\label{fig:pit1u}field \(u(x,t)\) tests basic projective integration.}
\includegraphics[width=\linewidth]{ProjInt/projInt1Test1u9}
\end{figure}

Set the number of interior points in the domain~\([-1,1]\), and the macroscale time step. 
\begin{matlab}
%}
function projInt1Test1
n=9
ts=0:2:6
%{
\end{matlab}
Set the initial condition to parabola or some skewed random positive values.
\begin{matlab}
%}
x=linspace(-1,1,n+2)';
%u0=1-x(2:n+1).^2 +1e-9*randn(n,1);
u0=rand(n,1).*(1+x(2:n+1));
%{
\end{matlab}
Projectively integrate in time one step for the moment with: 
rank-two DMD projection; 
guessed microscale time-step; and 
guessed numbers of transient~(15) and slow~(7) steps.
\begin{matlab}
%}
[us,uss,tss]=projInt1(@dudt,u0,ts,2,2/n^2,[round(n^2/5.4) 7])
%{
\end{matlab}
Plot the macroscale predictions to draw \autoref{fig:pit1u}.
\begin{matlab}
%}
plot(x,[0*ts;us;0*ts],'o-')
xlabel('space x'),ylabel('u(x,t)')
%matlab2tikz('projInt1Test1u.ltx','noSize',true)
print('-depsc2',['projInt1Test1u' num2str(n)])
%{
\end{matlab}
Also plot a surface of the microscale bursts as shown in \autoref{fig:pit1micro}.
\begin{figure}
\centering
\caption{\label{fig:pit1micro}field \(u(x,t)\) during each of the microscale bursts used in the projective integration.}
\includegraphics[width=\linewidth]{ProjInt/projInt1Test1micro9}
\end{figure}
\begin{matlab}
%}
surf(tss,x,[0*tss;uss;0*tss])
ylabel('space x'),xlabel('time t'),zlabel('u(x,t)')
view([40 30])
print('-depsc2',['projInt1Test1micro' num2str(n)])
%{
\end{matlab}

End the main function (not needed for new enough Matlab).
\begin{matlab}
%}
end
%{
\end{matlab}

\paragraph{The nonlinear PDE discretisation}
Code the simple centred difference discretisation of the nonlinear diffusion \pde\ with zero boundary values.
\begin{matlab}
%}
function ut=dudt(u,t)
n=length(u);
h=2/(n+1);
j=2:n+1;
u=[0;u;0];
ut=u(j).*(u(j+1)-2*u(j)+u(j-1))/h^2;
end
%{
\end{matlab}
%}
