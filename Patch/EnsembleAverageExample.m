% Example script to simulate an ensemble of solutions for
% heterogeneous diffusion in 1D on patches as an example
% application of patches in space.   JB & AJR, Nov 2017 --
% Mar 2019
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section[\texttt{ensembleAverageExample}: simulate an ensemble of solutions for heterogeneous diffusion in 1D \ldots]
{\texttt{ensembleAverageExample}: simulate an ensemble of solutions for heterogeneous diffusion in 1D on patches}
\label{sec:EnsembleAverageExample}
\localtableofcontents


\subsection{Introduction}
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
diffusivity configurations and then computing an ensemble of
field solutions, finally averaging over the ensemble of
fields to obtain the ensemble averaged field solution.

For a first comparison, we present a very similar example to
that of \cref{sec:HomogenisationExample}, but whereas
\cref{sec:HomogenisationExample} simulates using only one
diffusivity configuration, here we simulate over an
ensemble. For example, \cref{fig:HomogenisationCtsU} is
similar to \cref{fig:HomogenisationCtsUEnsAve}, but the
latter is an average of an ensemble of eight different
simulations with different diffusivity configurations,
whereas the former is simulated from just one diffusivity
configuration. The main difference between these two is that
the average over the ensemble caters for the heterogeneity
in the problem.

Much of this script is similar to that of
\cref{sec:HomogenisationExample}, but with additions to
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
ensemble average. The ensemble average caters for the
heterogeneity.}
\includegraphics[scale=0.9]{HomogenisationCtsUEnsAve}
\end{figure}%

Consider a lattice of values~\(u_i(t)\), with lattice
spacing~\(dx\), and governed by the heterogeneous diffusion 
\begin{equation}
\dot u_i=[c_{i-1/2}(u_{i-1}-u_i)+c_{i+1/2}(u_{i+1}-u_i)]/dx^2.
\label{eq:HomogenisationExample2}
\end{equation}
In this 1D space, the macroscale, homogenised, effective
diffusion should be the harmonic mean of these coefficients.
But suppose we do not know this.

\subsection{Script to simulate via stiff or projective integration}
Say there are four different diffusivities in our diffusive
medium, as defined here.
\begin{matlab}
%}
clear all
mPeriod = 4
rand('seed',1);
c = exp(4*rand(mPeriod,1))
cHomo = 1/mean(1./c)
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
reflected and translated configurations of \verb|patches.c|.
 Hence we increase the size of the diffusivity matrix to
\((\verb|nSubP-1)|\times \verb|nPatch|\times \verb|nVars|\).
\begin{matlab}
%}
patches.c = c((mod(round(patches.x(1:(end-1),:) ...
  /(patches.x(2)-patches.x(1))-0.5),mPeriod)+1));
if patches.EnsAve    
  nVars = mPeriod+(mPeriod>2)*mPeriod;
  patches.c = repmat(patches.c,[1,1,nVars]);    
  for sx = 2:mPeriod
    patches.c(:,:,sx) = circshift( ...
      patches.c(:,:,sx-1),[sx-1,0]);
   end;
   if nVars>2
     patches.c(:,:,(mPeriod+1):end) = flipud( ...
       patches.c(:,:,1:mPeriod)); 
   end;
end
%{
\end{matlab}

\paragraph{Conventional integration in time}
Set an initial condition, and here integrate forward in time
using a standard method for stiff systems. Integrate the
interface \verb|patchSmooth1()| (\cref{sec:patchSmooth1}) to
the microscale differential equations.
\begin{matlab}
%}
u0 = sin(patches.x)+0.2*randn(nSubP,nPatch);
if patches.EnsAve
  u0 = repmat(u0,[1,1,nVars]);
end
if ~exist('OCTAVE_VERSION','builtin')
    [ts,ucts] = ode15s( @patchSmooth1, [0 2/cHomo], u0(:));
else % octave version is slower
    [ts,ucts] = odeOcts(@patchSmooth1, [0 2/cHomo], u0(:));
end
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
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 10])
%print('-depsc2','ensAveExCtsU')
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
Now consider the interface, \verb|patchSmooth1()|, to the
time derivatives, and wrap around it the projective
integration \verb|PIRK2| (\cref{sec:PIRK2}), of bursts of
simulation from \verb|heteroBurst| (\cref{sec:heteroBurst}),
as illustrated by \cref{fig:HomogenisationUEnsAve}. The rest
of this code follows that of
\cref{sec:HomogenisationExample}, but as we now evaluate an
ensemble of field solutions, our final step is always an
ensemble average. 
 
Mark that edge of patches are not to be used in the
projective extrapolation by setting initial values to \nan.
\begin{matlab}
%}
disp('Now start Projective Integration')
u0([1 end],:) = nan;
%{
\end{matlab}
Set the desired macro- and microscale time-steps over the
time domain.
\begin{matlab}
%}
ts = linspace(0,2/cHomo,7)
bT = 3*( ratio*Len/nPatch )^2/cHomo
addpath('../ProjInt')
[us,tss,uss] = PIRK2(@heteroBurst, ts, u0(:), bT);
%{
\end{matlab}
\cref{fig:HomogenisationUEnsAve} shows an average of the
ensemble of macroscale predictions.
\begin{matlab}
%}
usAve = mean(reshape(us,size(us,1),length(xs(:)),nVars),3); 
ussAve = mean(reshape(uss,length(tss),length(xs(:)),nVars),3);
figure(2),clf
plot(xs(:),usAve','.')
ylabel('u(x,t)'), xlabel('space x')
legend(num2str(ts',3))
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 10])
%print('-depsc2','ensAveExU')
%{
\end{matlab}
Also plot a surface detailing the ensemble average
microscale bursts,
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
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 10])
%print('-depsc2','ensAveExMicro')
%{
\end{matlab}
End of the script.


\cref{sec:heteroDiff,sec:heteroBurst} list the addtional
functions used by this script.  Fin.
%}
