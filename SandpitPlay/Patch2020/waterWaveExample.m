% Simulate water waves on patches as an example of wave
% systems in the form h_t=u_x+... and u_u=h_x+...
% AJR, Nov 2017 -- Jul 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{waterWaveExample}: simulate a water wave PDE on patches}
\label{sec:waterWaveExample}
\localtableofcontents

\cref{fig:ps1WaveCtsUH} shows an example simulation in time
generated by the patch scheme applied to an ideal wave \pde\
\cite[]{Cao2013}. The inter-patch coupling is realised by
spectral interpolation of the mid-patch values to the patch
edges.
\begin{figure}
\centering \caption{\label{fig:ps1WaveCtsUH}water
depth~\(h(x,t)\) (above) and velocity field~\(u(x,t)\)
(below) of the gap-tooth scheme applied to the ideal linear
wave \pde~\cref{eq:genwaveqn} with \(f_1=f_2=0\). The
microscale random component to the initial condition
persists in the simulation---but the macroscale wave still
propagates.}
\includegraphics[scale=0.85]{waterWaveExample1CtsUH}
\end{figure}

This approach\ifcsname r@sec:idealWavePDE\endcsname, based
upon the differential equations coded in
\cref{sec:idealWavePDE},\fi\ may be adapted by a user to a
wide variety of 1D wave and wave-like systems.  For example,
the differential equations \ifcsname
r@sec:waterWavePDE\endcsname of \cref{sec:waterWavePDE} \fi
that describe the nonlinear microscale simulator of the
nonlinear shallow water wave \pde\ derived from the
Smagorinski model of turbulent flow \cite[]{Cao2012,
Cao2014b}.

\begin{devMan}
Often, wave-like systems are written in terms of two
conjugate variables, for example, position and momentum
density, electric and magnetic fields, and water
depth~\(h(x,t)\) and mean longitudinal velocity~\(u(x,t)\)
as herein. The approach developed in this section applies to
any wave-like system in the form
\begin{equation}
\D th=-c_1\D xu+f_1[h,u]
\quad\text{and}\quad
\D tu=-c_2\D xh+f_2[h,u],
\label{eq:genwaveqn}
\end{equation}
where the brackets indicate that the two nonlinear
functions~\(f_1\) and~\(f_2\) may involve various spatial
derivatives of the fields~\(h(x,t)\) and~\(u(x,t)\). For
example, \cref{sec:waterWavePDE} encodes a nonlinear
Smagorinski model of turbulent shallow water
\cite[e.g.]{Cao2012, Cao2014b} along an inclined flat bed:
let $x$~measure position along the bed and in terms of fluid
depth~$h(x,t)$ and depth-averaged longitudinal
velocity~$u(x,t)$ the model \pde{}s are
\begin{subequations}\label{eqs:patch:N}%
\begin{align}
\D th&=-\D x{(hu)}\,,\label{patch:Nh}
\\
\D tu&=0.985\left(\tan\theta-\D xh\right)
-0.003\frac{u|u|}{h} -1.045u\D xu +0.26h|u|\DD xu\,,
\label{patch:Nu}
\end{align}
\end{subequations}
where~$\tan\theta$ is the slope of the bed. The
\pde~\cref{patch:Nh} represents conservation of the fluid.
The momentum \pde~\cref{patch:Nu} represents  the effects of
turbulent bed drag~$u|u|/h$, self-advection~$u\D xu$,
nonlinear turbulent dispersion~$h|u|\DD xu$, and
gravitational hydrostatic forcing~$(\tan\theta-\D xh)$.
\cref{fig:ps1WaterWaveCtsUH} shows one simulation of this
system---for the same initial condition as
\cref{fig:ps1WaveCtsUH}.

\begin{figure}
\centering \caption{\label{fig:ps1WaterWaveCtsUH}water
depth~\(h(x,t)\) (above) and velocity field~\(u(x,t)\)
(below) of the gap-tooth scheme applied to the Smagorinski
shallow water wave \pde{}s~\cref{eqs:patch:N}. The
microscale random initial component decays where the water
speed is non-zero due to `turbulent' dissipation.}
\includegraphics[scale=0.85]{waterWaveExample2CtsUH}
\end{figure}%


For such wave-like systems, let's implement both a staggered
microscale grid and also staggered macroscale patches, as
introduced by \cite{Cao2014a} in their Figures~3 and~4,
respectively.




\subsection{Script code to simulate wave systems}
\label{sec:sc2waves}

This example script implements the following patch\slash
gap-tooth scheme (left-right arrows denote function
recursion).
\begin{enumerate}\def\itemsep{-1.5ex}
\item configPatches1, and add micro-information 
\item ode15s \into patchSmooth1 \into idealWavePDE
\item process results
\item ode15s \into patchSmooth1 \into waterWavePDE
\item process results
\end{enumerate}
Establish the global data struct~\verb|patches| for the
\pde{}s~\cref{eq:genwaveqn} (linearised) solved on
\(2\pi\)-periodic domain, with eight patches, each patch of
half-size ratio~\(0.2\), with eleven micro-grid points
within each patch, and spectral interpolation~(\(-1\)) of
`staggered' macroscale patches to provide the edge-values of
the inter-patch coupling conditions.
\begin{matlab}
%}
clear all
global patches
nPatch = 8
ratio = 0.2
nSubP = 11 %of the form 4*n-1
Len = 2*pi;
configPatches1(@idealWavePDE,[0 Len],nan,nPatch,-1,ratio,nSubP);
%{
\end{matlab}

Identify which micro-grid points are \(h\)~or~\(u\) values
on the staggered micro-grid. Also store the information in
the struct~\verb|patches| for use by the time derivative
function.
\begin{matlab}
%}
uPts = mod( bsxfun(@plus,(1:nSubP)',(1:nPatch)) ,2);
hPts = find(uPts==0);
uPts = find(uPts==1);
patches.hPts = hPts; patches.uPts = uPts;
%{
\end{matlab}

Set an initial condition of a progressive wave, and check
evaluation of the time derivative. The capital
letter~\verb|U| denotes an array of values merged from
both~\(u\) and~\(h\) fields on the staggered grids (here
with some optional microscale wave noise).
\begin{matlab}
%}
U0 = nan(nSubP,nPatch);
U0(hPts) = 1+0.5*sin(patches.x(hPts));
U0(uPts) = 0+0.5*sin(patches.x(uPts));
U0 = U0+0.02*randn(nSubP,nPatch);
%{
\end{matlab}

\paragraph{Conventional integration in time} Integrate in
time using standard \script\ stiff integrators. Here do the
two cases of the ideal wave and the water wave equations in
the one loop.
\begin{matlab}
%}
for k = 1:2
%{
\end{matlab}
When using \verb|ode15s|\slash\verb|lsode| we subsample the
results because micro-grid scale waves do not dissipate and
so the integrator takes very small time-steps for all time.
\begin{matlab}
%}
if ~exist('OCTAVE_VERSION','builtin')
    [ts,Ucts] = ode15s( @patchSmooth1,[0 4],U0(:));
    ts = ts(1:5:end);
    Ucts = Ucts(1:5:end,:);
else % octave version is slower
    [ts,Ucts] = odeOcts(@patchSmooth1,[0 4],U0(:));
end
%{
\end{matlab}
Plot the simulation.
\begin{matlab}
%}
  figure(k),clf
  xs = squeeze(patches.x);  xs([1 end],:) = nan;
  mesh(ts,xs(hPts),Ucts(:,hPts)'),hold on
  mesh(ts,xs(uPts),Ucts(:,uPts)'),hold off
  xlabel('time t'), ylabel('space x'), zlabel('u(x,t) and h(x,t)')
  axis tight, view(70,45)
%{
\end{matlab}
Optionally save the plot to file.
\begin{matlab}
%}
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 10])
%print('-depsc2',[mfilename num2str(k) 'CtsUH'])
%{
\end{matlab}

For the second time through the loop, change to the
Smagorinski turbulence model~\cref{eqs:patch:N} of shallow
water flow, keeping other parameters and the initial
condition the same. 
\begin{matlab}
%}
  patches.fun = @waterWavePDE;
end
%{
\end{matlab}

\paragraph{Could use projective integration} As yet a simple
implementation appears to fail, so it needs more exploration
and thought.
End of the main script.



\input{../Patch/idealWavePDE.m}

\input{../Patch/waterWavePDE.m}

Fin.
\end{devMan}
%}