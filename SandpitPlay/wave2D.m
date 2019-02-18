% Simulate the simple wave PDE in 2D on patches.
% First it checks the spectrum of the system.
% AJR, Nov 2018
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{wave2D}: example of a wave on patches in 2D}
\label{sec:wave2D}
\localtableofcontents

For \(u(x,y,t)\), test and simulate the simple wave PDE in 2D space:
\begin{equation*}
\DD tu=\delsq u\,.
\end{equation*}
This script shows one way to get started: a user's script
may have the following three steps (arrows indicate function
recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches2 
\item ode15s integrator \into patchSmooth2 \into wavePDE
\item process results
\end{enumerate}


\begin{body}
Establish global patch data struct to interface with a
function coding the wave \pde: to be solved
on \(2\pi\)-periodic domain, with \(9\times9\) patches,
spectral interpolation couples the patches, each patch of
half-size ratio~\(0.25\), and with \(5\times5\) points
within each patch.
\begin{matlab}
%}
clear all, close all
global patches
nSubP = 5;
nPatch = 9;
configPatches2(@wavePDE,[-pi pi], nan, nPatch, 0, 0.1, nSubP);
%{
\end{matlab}



\subsubsection{Check on the linear stability of the wave PDE}
Set a zero equilibrium as basis.
Then find the patch-interior points as the only ones to vary in order to construct the Jacobian.
\begin{matlab}
%}
disp('Check linear stability of the wave scheme')
uv0=zeros(nSubP,nSubP,nPatch,nPatch,2);
uv0([1 end],:,:,:,:,:)=nan;
uv0(:,[1 end],:,:,:,:)=nan;
i=find(~isnan(uv0));
%{
\end{matlab}
Now construct the Jacobian. Since linear wave \pde, use large perturbations.
\begin{matlab}
%}
small=1;
jac=nan(length(i));
sizejac=size(jac)
for j=1:length(i)
  uv=uv0(:);
  uv(i(j))=uv(i(j))+small;
  tmp=patchSmooth2(0,uv)/small;
  jac(:,j)=tmp(i);
end
%{
\end{matlab}
Now explore the eigenvalues a little: find the ten with the biggest real-part; if small enough, then the method may be good.
\begin{matlab}
%}
evals=eig(jac);
nEvals=length(evals)
[~,k]=sort(-abs(real(evals)));
evalsWithBiggestRealPart=evals(k(1:10))
if abs(real(evals(k(1))))>1e-4
    warning('eigenvalue failure: real-part > 1e-4')
    return, end
%{
\end{matlab}
Check eigenvalues close to true waves of the \pde\ 
(not yet the micro-discretised equations).
\begin{matlab}
%}
kwave=0:(nPatch-1)/2;
freq=sort(reshape(sqrt(kwave'.^2+kwave.^2),1,[]));
freq= freq(diff([-1 freq])>1e-9);
freqerr=[freq; min(abs(imag(evals)-freq))]
%{
\end{matlab}





\subsubsection{Execute a simulation}
Set a Gaussian initial condition using auto-replication of
the spatial grid: here \verb|u0| and~\verb|v0| are in the form required for computation: \(n_x\times n_y\times N_x\times N_y\).
\begin{matlab}
%}
x = reshape(patches.x,nSubP,1,[],1); 
y = reshape(patches.y,1,nSubP,1,[]);
u0 = exp(-x.^2-y.^2); 
v0 = zeros(size(u0));
%{
\end{matlab}
Initiate a plot of the simulation using only the microscale
values interior to the patches: set \(x\)~and \(y\)-edges to
\verb|nan| to leave the gaps. Start by showing the initial
conditions of \cref{fig:configPatches2ic} while the
simulation computes.
To mesh/surf plot we need to `transpose' to size \(n_x\times N_x\times n_y\times N_y\), then reshape to size  \(n_x\cdot N_x\times n_y\cdot N_y\).\begin{matlab}
%}
x=patches.x; y=patches.y;
x([1 end],:)=nan; y([1 end],:)=nan;
u = reshape(permute(u0,[1 3 2 4]), [numel(x) numel(y)]);
usurf = surf(x(:),y(:),u');
axis([-3 3 -3 3 -0.5 1]), view(60,40)
xlabel('space x'), ylabel('space y'), zlabel('u(x,y)')
drawnow
set(gcf,'paperposition',[0 0 14 10])
print('-depsc','wave2Dic.eps')
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:wave2Dic}initial
field~\(u(x,y,t)\) at time \(t=0\) of the patch scheme
applied to the simple wave~\pde:
\cref{fig:wave2Dt6} plots the computed field at time
\(t=6\).}
\includegraphics[scale=0.8]{Patch/wave2Dic}
\end{figure}
Integrate in time using standard functions.
\begin{matlab}
%}
disp('Wait while we simulate u_t=v, v_t=u_xx+u_yy')
[ts,uvs] = ode15s(@patchSmooth2,[0 1],[u0(:);v0(:)]);
%{
\end{matlab}
Animate the computed simulation to end with
\cref{fig:wave2Dt6}.  Subsample to plot at most 200 times.
\begin{matlab}
%}
di = ceil(length(ts)/200);
for i = [1:di:length(ts)-1 length(ts)]
  uv = patchEdgeInt2(uvs(i,:)); 
  uv = reshape(permute(uv,[1 3 2 4 5]), [numel(x) numel(y) 2]);
  usurf.ZData = uv(:,:,1)';
  title(['wave PDE on patches: time = ' num2str(ts(i))])
  pause(0.1)
end
title('')
set(gcf,'paperposition',[0 0 14 10])
print('-depsc',['wave2Dt' num2str(ts(end)) '.eps'])
%{
\end{matlab}
\begin{figure}
\centering
\caption{\label{fig:wave2Dt6}field~\(u(x,y,t)\) at
time \(t=6\) of the patch scheme applied to the simple wave~\pde\ with initial condition in
\cref{fig:wave2Dic}.}
\includegraphics[scale=0.8]{Patch/wave2Dt6}
\end{figure}



\subsubsection{Example of simple wave PDE inside patches}
As a microscale discretisation of \(u_{tt}=\delsq(u)\),
so code \(\dot u_{ijkl}=v_{ijkl}\) and \(\dot v_{ijkl} =\frac1{\delta x^2}
(u_{i+1,j,k,l} -2u_{i,j,k,l} +u_{i-1,j,k,l}) + \frac1{\delta y^2}
(u_{i,j+1,k,l} -2u_{i,j,k,l} +u_{i,j-1,k,l})\).
\begin{matlab}
%}
function uvt = wavePDE(t,uv,x,y)
  if ceil(t+1e-7)-t<2e-2, simTime=t, end %track progress
  dx=diff(x(1:2));  dy=diff(y(1:2));   % micro-scale spacing
  i=2:size(uv,1)-1;  j=2:size(uv,2)-1; % interior patch-points
  uvt = nan(size(uv));  % preallocate storage
  uvt(i,j,:,:,1) = uv(i,j,:,:,2);
  uvt(i,j,:,:,2) = diff(uv(:,j,:,:,1),2,1)/dx^2 ...
                  +diff(uv(i,:,:,:,1),2,2)/dy^2;
end
%{
\end{matlab}
\end{body}
%}