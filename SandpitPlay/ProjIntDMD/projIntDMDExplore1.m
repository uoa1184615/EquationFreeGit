% Explore the time integration function projIntDMD() 
% on discretisation of a nonlinear diffusion PDE.  
% AJR, Jan 2018
%!TEX root = doc.tex
%{
\section{\texttt{projIntDMDExplore1}: explore effect of varying parameters}
\label{sec:pi1eevp}
\localtableofcontents

Seek to simulate the nonlinear diffusion \pde\
\begin{equation*}
\D tu=u\DD xu\quad \text{such that }u(\pm 1)=0,
\end{equation*}
with random positive initial condition.

Set the number of interior points in the domain~\([-1,1]\), and the macroscale time-step. 
\begin{matlab}
%}
function projIntDMDExplore1
n=9
ts=0:2:6
dt=2/n^2
ICNoise=0.3
%{
\end{matlab}
Very strangely, the results from Matlab and Octave are different for the zero noise case!????  It should be deterministic.  Significantly different in that Matlab fails more often.

\begin{figure}
\caption{\label{fig:explor1icn0}errors in the projective integration of the nonlinear diffusion \pde\ from initial conditions that are on the slow manifold.  Plotted are stereo views of isosurfaces in parameter space: the first row is after the first projective step; the second row after the second step.  }
\centering
\includegraphics[width=\linewidth]{explore1icn0}
\end{figure}%
\autoref{fig:explor1icn0} shows the parameter variations when the system is already on the slow manifold.
The picture after two time-steps, bottom row, appears clearer than for one time-step.
The errors do not vary with rank provided that it is\({}\geq2\).
There is only a very weak dependence upon the length of the burst being analysed---and that could be due to reduction in the gap.
There is a weak dependence upon the transient-time, but only by a factor of two across the domain considered.


\begin{figure}
\caption{\label{fig:explor1icn3}errors in the projective integration of the nonlinear diffusion \pde\ from initial conditions with noise~\texttt{0.3*rand}.  Plotted are stereo views of isosurfaces in parameter space: the first row is after the first projective step; the second row after the second step.  }
\centering
\includegraphics[width=\linewidth]{explore1icn3}
\end{figure}%
With the addition of a noisy initial conditions, \autoref{fig:explor1icn3}, the rank has an effect, and the transient-time appears to be a slightly stronger influence.
I suspect this means that we need to allow the initial burst to have a longer transient time than subsequent bursts.
Initial conditions may typically need a longer `healing' time.
Thus code an extra \verb|timeStep| parameter.



Set the initial condition to parabola or some skewed random positive values.
Without noise this initial condition is already on the slow manifold so only little reason for transient time.
\begin{matlab}
%}
x=linspace(-1,1,n+2)';
u0=(0.5+ICNoise*rand(n+2,1)).*(1-x.^2);
%{
\end{matlab}

First find a reference solution of the microscale dynamics over all time.
\begin{matlab}
%}
[Us,Uss,Tss]=projIntDMD(@dudt,u0,ts,2,dt,[0 2]);
%{
\end{matlab}

Set up various combinations of parameters.
\begin{matlab}
%}
[rank,trant,slowt]=meshgrid(1:5,[1 2 4 6 8]*0.05,[2 4 8 12 16]*dt);
ps=[rank(:) trant(:) slowt(:)];
%{
\end{matlab}

Projectively integrate in time with various parameters.
\begin{matlab}
%}
errs=[]; relerrs=[];
for p=ps'
[us,uss,tss]=projIntDMD(@dudt,u0,ts,p(1),dt,p(2:3));
%{
\end{matlab}
Plot the macroscale predictions 
\begin{matlab}
%}
if 0
  clf,plot(x,Us,'o-',x,us,'x--')
  xlabel('space x'),ylabel('u(x,t)')
  pause(0.01)
end
%{
\end{matlab}
Accumulate errors as function of time.
\begin{matlab}
%}
err=sqrt(sum((us-Us).^2))
errs=[errs;err];
relerrs=[relerrs;err./sqrt(sum(Us.^2))];
%{
\end{matlab}

End the loop over parameters.
\begin{matlab}
%}
end
%{
\end{matlab}

Stereo view of isosurfaces of errors after both one and two time-steps.
The three surfaces are the quartiles of the errors, coloured accordingly, but with a little extra colour from position for clarity.
\begin{matlab}
%}
clf()
vals=nan(size(rank));
for k=1:2
  vals(:)=errs(:,k+1);
  q=prctile(vals(:),(0:4)*25)  
  for j=1:2, subplot(2,2,j+2*(k-1)),hold on
    for i=2:4 % draw three quartiles
      isosurface(rank,trant,slowt,vals,q(i) ...
      ,q(i)+0.03*(rank/10-trant+slowt))
    end,hold off
    xlabel('rank'),ylabel('transient'),zlabel('analysed')
    colorbar
    set(gca,'view',[57-j*5,30])
  end%j
end%k
%{
\end{matlab}
Save to file
\begin{matlab}
%}
print('-depsc2',['explore1icn' num2str(ICNoise*10)])
%{
\end{matlab}



End the main function (not needed for new enough Matlab).
\begin{matlab}
%}
end
%{
\end{matlab}





\paragraph{The nonlinear PDE discretisation}
Code the simple centred difference discretisation of the nonlinear diffusion \pde\ with constant (usually zero) boundary values.
\begin{matlab}
%}
function ut=dudt(t,u)
n=length(u);
dx=2/(n-1);
j=2:n-1;
ut=[0
    u(j).*(u(j+1)-2*u(j)+u(j-1))/dx^2
    0];
end
%{
\end{matlab}
%}
