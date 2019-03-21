% Example script to simulate an ensemble of solutions for
% heterogeneous diffusion in 1D on patches as an example
% application of patches in space.   JB & AJR, Nov 2017 --
% Mar 2019
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section[\texttt{EnsembleAverageExample}: simulate an ensemble of solutions for heterogeneous diffusion in 1D \ldots]
{\texttt{EnsembleAverageExample}: simulate an ensemble of solutions for heterogeneous diffusion in 1D on patches}
\label{sec:EnsembleAverageExample}
\localtableofcontents

This example is an extension of the homogenisation example
of \cref{sec:HomogenisationExample} for heterogeneous
diffusion. In cases where the periodicity of the 
heterogeneous diffusion is known, then
\cref{sec:HomogenisationExample} provides a efficient patch
dynamics simulation. However, if the diffusion is not
completely known or is stochastic, then we cannot choose
ideal patch and core sizes as described by
\cite{Bunder2013b} and applied in
\cref{sec:HomogenisationExample}. In this case,
\cite{Bunder2013b} recommend constructing an ensemble of
diffusivity configurations and then computing an ensemble
of field solutions, finally averaging over the ensemble of
fields to obtain the ensemble averaged field solution.

For a first comparison, we present a very similar example to
that presented in \cref{sec:HomogenisationExample}, but
whereas \cref{sec:HomogenisationExample} simulates using
only one diffusivity configuration, here we simulate over an
ensemble. For example, \cref{fig:HomogenisationCtsUEnsAve}
is similar to \cref{fig:HomogenisationCtsU}, but the former
is an ensemble average of an ensemble of eight different
simulations with different diffusivity configurations and
the latter is simulated from just one diffusivity
configuration. The main difference between these two is that
the  average over the ensemble removes any heterogeneity in
the solution.

Much of this script is similar to that of
\cref{sec:HomogenisationExample}, but with some additions to
manage the ensemble. The first part of the script implements
the following gap-tooth scheme (left-right arrows denote
function recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches1 
\item ode15s \into patchSmooth1 \into heteroDiff
\item process results
\end{enumerate}
\begin{figure}
\centering \caption{\label{fig:HomogenisationCtsUEnsAve}the
diffusing field~\(u(x,t)\) in the patch (gap-tooth) scheme
applied to microscale heterogeneous diffusion with an
ensemble average. The ensemble average smooths out the
heterogeneous diffusion.}
\includegraphics[scale=0.9]{HomogenisationCtsUEnsAve}
\end{figure}%

Consider a lattice of values~\(u_i(t)\), with lattice
spacing~\(dx\), and governed by the heterogeneous diffusion 
\begin{equation}
\dot u_i=[c_{i-1/2}(u_{i-1}-u_i)+c_{i+1/2}(u_{i+1}-u_i)]/dx^2.
\label{eq:HomogenisationExample}
\end{equation}
In this 1D space, the macroscale, homogenised, effective
diffusion should be the harmonic mean of these coefficients.
But we do not have full knowledge of these coefficients.

\subsection{Script to simulate via stiff or projective integration}
Say we only know four diffusivities in our diffusion problem, 
as defined here (which are the same as those given in 
\cref{sec:HomogenisationExample}).
\begin{matlab}
%}
clear all
mPeriod = 4
rng('default'); rng(1);
cDiff = exp(4*rand(mPeriod,1))
cHomo = 1/mean(1./cDiff)
%{
\end{matlab}

The chosen parameters are the same as
\cref{sec:HomogenisationExample}, but here we also introduce
the Boolean \verb|patches.EnsAve| which determines whether
or not we construct an ensemble average of diffusivity
configurations. Setting \verb|patches.EnsAve=0| simulates
the same problem as in \cref{sec:HomogenisationExample}.
\begin{matlab}
%}
global patches
nPatch = 9
ratio = 0.2
nSubP = 11
Len = 2*pi;
ordCC = 4;
patches.nCore = 3; 
patches.ratio = ratio*(nSubP - patches.nCore)/(nSubP - 1);
configPatches1(@heteroDiff,[0 Len],nan,nPatch ...
    ,ordCC,patches.ratio,nSubP);
patches.EnsAve = 1;
%{
\end{matlab}

In the case of ensemble averaging, \verb|nVars| is the size
of the ensemble (for the case of no ensemble averaging
\verb|nVars| is the number of different field variables,
which in this example is \(\verb|nVars|=1\)) and we use the
ensemble described by \cite{Bunder2013b} which includes all
reflected and translated configurations of
\verb|patches.cDiff|. We must increase the size of the
diffusivity matrix to \((\verb|nSubP-1)|\times
\verb|nPatch|\times \verb|nVars|\).
\begin{matlab}
%}
patches.cDiff = cDiff((mod(round(patches.x(1:(end-1),:) ...
  /(patches.x(2)-patches.x(1))-0.5),mPeriod)+1));
if patches.EnsAve    
  nVars = mPeriod+(mPeriod>2)*mPeriod;
  patches.cDiff = repmat(patches.cDiff,[1,1,nVars]);    
  for sx = 2:mPeriod
    patches.cDiff(:,:,sx) = circshift( ...
      patches.cDiff(:,:,sx-1),[sx-1,0]);
   end;
   if nVars>2
     patches.cDiff(:,:,(mPeriod+1):end) = flipud( ...
       patches.cDiff(:,:,1:mPeriod)); 
   end;
end
%{
\end{matlab}

\paragraph{Conventional integration in time}
Set an initial condition, and here integrate forward in time
using a standard method for stiff systems. Integrate the
interface \verb|patchSmooth1| (\cref{sec:patchSmooth1}) to
the microscale differential equations.
\begin{matlab}
%}
u0 = sin(patches.x)+0.2*randn(nSubP,nPatch);
if patches.EnsAve
  u0 = repmat(u0,[1,1,nVars]);
end
[ts,ucts] = ode15s(@patchSmooth1, [0 2/cHomo], u0(:));
ucts = reshape(ucts,length(ts),length(patches.x(:)),[]);
%{
\end{matlab}

Plot the ensemble averaged simulation in
\cref{fig:HomogenisationCtsUEnsAve}.
\begin{matlab}
%}
if patches.EnsAve % calculate the ensemble average
  uctsAve = mean(ucts,3);
else
  uctsAve = ucts;
end
figure(1),clf
xs = patches.x;  xs([1 end],:) = nan;
mesh(ts,xs(:),uctsAve'),  view(60,40)
xlabel('time t'), ylabel('space x'), zlabel('u(x,t)')
set(gcf,'PaperPosition',[0 0 14 10]);% cm
print('-depsc2','Figs/HomogenisationCtsUEnsAve')
%{
\end{matlab}

\paragraph{Use projective integration in time}
\begin{figure}
\centering \caption{\label{fig:HomogenisationUEnsAve}field
\(u(x,t)\) shows basic projective integration of patches of
heterogeneous diffusion with an ensemble average: different
colours correspond to the times in the legend. Once
transients have decayed, this field solution is smooth due
to the ensemble average.}
\includegraphics[scale=0.9]{HomogenisationUEnsAve}
\end{figure}%
Now take \verb|patchSmooth1|, the interface to the time
derivatives, and wrap around it the projective integration
\verb|PIRK2| (\cref{sec:PIRK2}), of bursts of simulation
from \verb|heteroBurst| (\cref{sec:heteroBurst}), as
illustrated by \cref{fig:HomogenisationUEnsAve}. The rest of
this code follows that of \cref{sec:HomogenisationExample},
but as we now evaluate an ensemble of field solutions, our
final step is always an ensemble average. 
 
Mark that edge of patches are not to be used in the
projective extrapolation by setting initial values to \nan.
\begin{matlab}
%}
u0([1 end],:) = nan;
%{
\end{matlab}
Set the desired macro- and microscale time-steps over the
time domain.
\begin{matlab}
%}
ts = linspace(0,2/cHomo,7)
bT = 3*( ratio*Len/nPatch )^2/cHomo
addpath('../ProjInt','../SandpitPlay/RKint')
[us,tss,uss] = PIRK2(@heteroBurst, ts, u0(:), bT);
%{
\end{matlab}
Plot an average of the ensemble of macroscale predictions 
to draw \cref{fig:HomogenisationUEnsAve}.
\begin{matlab}
%}
usAve = mean(reshape(us,size(us,1),length(xs(:)),nVars),3); 
ussAve = mean(reshape(uss,length(tss),length(xs(:)),nVars),3);
figure(2),clf
plot(xs(:),usAve','.')
ylabel('u(x,t)'), xlabel('space x')
legend(num2str(ts',3))
set(gcf,'PaperPosition',[0 0 14 10]);% cm
print('-depsc2','Figs/HomogenisationUEnsAve')
%{
\end{matlab}
Also plot a surface detailing the ensemble average
microscale bursts as shown
\cref{fig:HomogenisationMicroEnsAve}.
\begin{figure}
\centering
\caption{\label{fig:HomogenisationMicroEnsAve}stereo pair of
ensemble averaged fields~\(u(x,t)\) during each of the
microscale bursts used in the projective integration.}
\includegraphics[scale=0.9]{HomogenisationMicroEnsAve}
\end{figure}
\begin{matlab}
%}
figure(3),clf
for k = 1:2, subplot(1,2,k)
  surf(tss,xs(:),ussAve',  'EdgeColor','none')
  ylabel('x'), xlabel('t'), zlabel('u(x,t)')
  axis tight, view(126-4*k,45)
end
set(gcf,'PaperPosition',[0 0 14 6]);% cm
print('-depsc2','Figs/HomogenisationMicroEnsAve')
%{
\end{matlab}
End of the script.


\subsection{\texttt{heteroDiff()}: heterogeneous diffusion}
\label{sec:heteroDiff}
This function codes the \pde\ of the heterogeneous diffusion
inside the patches, coding it as \ode{}s on a microscale
lattice.  For 2D input arrays~\verb|u| and~\verb|x| (via
edge-value interpolation of \verb|patchSmooth1|,
\cref{sec:patchSmooth1}), computes the time
derivative~\cref{eq:HomogenisationExample} at each point in
the interior of a patch, output in~\verb|ut|.  The column
vector (or possibly array) of diffusion coefficients~\(c_i\)
have previously been stored in struct~\verb|patches|.
\begin{matlab}
%}
function ut = heteroDiff(t,u,x)
  global patches
  dx = diff(x(2:3)); % space step
  i = 2:size(u,1)-1; % interior points in a patch
  ut = nan(size(u)); % preallocate output array
  ut(i,:,:) = diff(patches.cDiff.*diff(u))/dx^2; 
end% function
%{
\end{matlab}



\subsection{\texttt{heteroBurst()}: a burst of heterogeneous diffusion}
\label{sec:heteroBurst}
This code integrates in time the derivatives computed by
\verb|heteroDiff| from within the patch coupling of
\verb|patchSmooth1|.  Try four possibilities:
\begin{itemize}
\item \verb|ode23| generates `noise' that is unsightly at
best and may be ruinous;
\item \verb|ode45| is similar to \verb|ode23|, but with
reduced noise;
\item \verb|ode15s| does not cater for the \nan{}s in some
components of~\verb|u|;
\item \verb|rk2int| simple specified step integrator, but
may require inefficiently small time-steps.
\end{itemize}
\begin{matlab}
%}
function [ts, ucts] = heteroBurst(ti, ui, bT) 
  switch '45'
  case '23',  [ts,ucts] = ode23(@patchSmooth1,[ti ti+bT],ui(:));
  case '45',  [ts,ucts] = ode45(@patchSmooth1,[ti ti+bT],ui(:));
  case '15s', [ts,ucts] = ode15s(@patchSmooth1,[ti ti+bT],ui(:));
  case 'rk2', ts = linspace(ti,ti+bT,200)';
              ucts = rk2int(@patchSmooth1,ts,ui(:));
  end
end
%{
\end{matlab}
Fin.
%}
