% chanDispTimescale explores the timescales in the patch scheme applied to 2D shear dispersion in a long
% channel with 1D patches.    AJR, Dec 2020
%%%%%%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{chanDispTimescale}: example of estimating subpatch micro-time-scales by one patch of sdvection-diffusion in a channel}
\label{sec:chanDispTimescale}
\localtableofcontents

\cref{sec:chanDispSpmd} describes the physical problem.
\cref{sec:chanDispMicro} codes the micro-scale discretisation for the patch scheme computation.

Create and analyse the Jacobian to find slowest decay rates of sub-patch dynamics---confirmed with a simulation.
But not yet clear how to use the information??



First, clear all to remove any existing globals, old
composites, etc. 
\begin{matlab}
%}
clear all
%{
\end{matlab}
The micro-scale \pde\ is evaluated at positions~\(y_j\)
across the channel, \(|y|<1\)\,.  The even indexed points
are the collocation points for the \pde, whereas the even
indexed points are the half-grid points for specification of
\(y\)-diffusivities.
\begin{matlab}
%}
ny = 7
y = linspace(-1,1,2*ny+1);
yj = y(2:2:end);
%{
\end{matlab}
Set micro-scale advection (array~1) and diffusivity
(array~2) with (roughly) parabolic shape
\cite[e.g.]{Watt94b, MacKenzie03}. Here modify the parabola
by a heterogeneous log-normal factor with specified period
along the channel: modify the strength of the heterogeneity
by the coefficient of~\verb|randn| from zero to perhaps one:
coefficient~\(0.3\) appears a good moderate value. Remember
that \verb|configPatches1| reshapes \verb|cHetr| to~2D.
\begin{matlab}
%}
mPeriod = 4
cHetr = shiftdim([3/2 1],-1).*(1-y.^2) ...
        .*exp(0.99*randn([mPeriod 2*ny+1 2]));
%{
\end{matlab}


Configure the patch scheme with some arbitrary choices of
domain, size ratios, etc---as for full problem.  But then use just one patch in order to focus on the sub-patch dynamics, and scale domain length down to match.   The order of interpolation is irrelevant when there is one patch so just use second-order, say. 
Set \verb|patches| information to be global so the info can
be used without being explicitly passed as
arguments.  
\begin{matlab}
%}
global patches
nPatch=15
nSubP=2+mPeriod
ratio=0.4
Len=nPatch/ratio 
patches = configPatches1(@chanDispMicro, [0 Len/nPatch], nan ...
    , 1, 2, ratio, nSubP, 'EdgyInt',true  ...
    ,'hetCoeffs',cHetr );
%{
\end{matlab}
Additional parameters to \verb|patches|.
\begin{matlab}
%}
Peclet = 10, patches.Pe = Peclet;
%{
\end{matlab}


\subsection{Eigenvalues of the Jacobian}
Set zero to be the reference equilibrium in this linear problem.
Put NaNs on the patch-edges.
\begin{matlab}
%}
c0 = zeros(size(patches.x+0*yj));
c0([1 end],:,:)=nan;
%{
\end{matlab}
The dynamic variables are the locations of non-NaNs. 
\begin{matlab}
%}
i=find(~isnan(c0));
nJac=length(i);
%{
\end{matlab}
Construct the Jacobian column-wise from the transform of a complete set of unit basis vectors (as this is linear problem).
\begin{matlab}
%}
Jac=nan(nJac);
for j=1:nJac
  cj=c0; cj(i(j))=1;
  dcjdt=patchSmooth1(0,cj);
  Jac(:,j)=dcjdt(i);
end
%{
\end{matlab}

\begin{matlab}
%}
evals=eig(Jac);
[~,i]=sort(-real(evals));
decayRates = evals(i(1:5))
betaRateEig = - decayRates(2)
%{
\end{matlab}
Possibly plot pseudo-spectra, just for fun!
\begin{matlab}
%}
if exist('psa')==2
disp('plotting spectra and pseudo-spectra (Trefethen, 1999)')
figure(2),psa(Jac)
else 
disp('download psa() (Trefethen, 1999) to plot pseudo-spectra')
end
%{
\end{matlab}



\subsection{Simulate heterogeneous advection-diffusion}
Set random initial conditions of a simulation.
\begin{matlab}
%}
    c0 = rand(size(patches.x+0*yj));
    dc0dt = patchSmooth1(0,c0);
%{
\end{matlab}

Integrate in time, either via the automatic \verb|ode23| or
via \verb|RK2mesoPatch| which reduces communication between
patches. By default, \verb|RK2mesoPatch| does ten
micro-steps for each specified meso-step in~\verb|ts|.  For
stability: with noise up to~\(0.3\), need micro-steps less
than~\(0.005\); with noise~\(1\), need micro-steps less
than~\(0.0015\).
\begin{matlab}
%}
ts=linspace(0,2); 
%    [cs,uerrs] = RK2mesoPatch(ts,c0);
[ts,cs] = ode15s(@patchSmooth1,ts,c0(:));
cs=reshape(cs,[length(ts) size(c0)]);
cs=cs(:,2:end-1,:);
%finalCs=squeeze(cs(end,:,:))
stdy=squeeze(std(cs,[],2));
figure(1),semilogy(ts,stdy,'.')
j=[round(length(ts)*0.66) length(ts)];
betaRates = - diff(log(stdy(j,:)))./diff(ts(j));
betaRateSimulation = [mean(betaRates) std(betaRates)]
%{
\end{matlab}


Fin.
%}
