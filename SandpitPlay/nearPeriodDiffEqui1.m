% Explore heterogeneous diffusion in 1D on patches to
% compare with approach in arxiv:2308.07563 by authors CDE. 
% The microscale period epsilon is to be a little different
% from the patch size eta. Then explore accuracy via 
% forced equilibria. AJR, Aug 2023
%{
\documentclass[11pt,a5paper]{article}
\IfFileExists{ajr.sty}{\usepackage{ajr}}{}
\usepackage{amsmath,mybiblatex,defns,matWeave,mycleveref}



\title{\texttt{nearPeriodDiffEqui1}: errors in patch scheme
for equilibrium of a 1D heterogeneous diffusion with nearly
correct period}
%\label{sec:nearPeriodDiffEqui1}
\author{A.J. Roberts}
\date{\today}

\begin{document}

\maketitle


Explore heterogeneous diffusion in 1D on patches to compare
with approach in arxiv:2308.07563 by authors CDE. The
microscale period~\(\epsilon\) is to be a little different
from the cell (patch) size~\(\eta\). Then explore accuracy
via forced equilibria. Here we use cells that are patches in
the \emph{equation-free patch scheme}
\cite[e.g.,][]{Roberts2023a, Bunder2020a, Samaey08}. We
invoke functions from the \emph{Equation-Free Toolbox}
\cite[]{Maclean2020a}.

Suppose the spatial microscale lattice is at points~\(x_i\),
with constant spacing~\(dx\). With dependent
variables~\(u_i(t)\), let's seek the equilibriun of the
forced the microscale lattice diffusion system
\begin{equation}
\D t{u_{i}}= \frac1{dx^2}\delta[c_{i-1/2}\delta u_{i}] -u_i+f_i,
\label{eq:hetroDiff}
\end{equation}
in terms of the
centred difference operator~\(\delta\), and for some
time-constant forcing~\(f_i\). The system has a microscale
heterogeneity via the coefficients~\(c_{i+1/2}\) which has
periodicity~\(\epsilon\) in~\(x\). Instead of varying
cell-size for fixed heterogeneity period, here we explore
results for various periods~\(\epsilon\) at fixed
cell-size~\(\eta\).


\begin{figure}
\centering
\caption{\label{fig:ueq}example solutions of heterogeneous
diffusion equilibrium with forcing \(f(x)=\sin x\).   By
symmetry the plot shows half the domain.  As shown, solve
with \(10,30,90\)~patches\slash cells.  The bottom
\(90\)~patch equilibrium is the exact reference solution. 
The top two equilibria are for patch ratios \(r=1/9,1/3\)
respectively: computation is done only on the fraction~\(r\)
of the domain.}
\includegraphics[width=\textwidth]{nearPeriodDiffEqui1ueq}
\end{figure}


\subsection{Code various numbers of patches over domain}
Establish system of length~\(2\pi\). Explore various number
of patches.
\begin{matlab}
%}
clear all
Xlim = [-pi pi]
lMax = 3
nPatches = 10 *3.^(0:lMax-1)
maxDetune = 100
%{
\end{matlab}
Set up microgrid parameters, and set strength of
heterogeneity (abs-value less than one).  CDE used cell size
\(\eta/\epsilon\in[1,50]\), and at least \(4\,096\) points,
so here with, say, \(90\)~cells that is over \(45\)~points
per cell (per patch).   For computational speed, use less.
\begin{matlab}
%}
mPerPatch = 12 
eta = diff(Xlim)/nPatches(lMax)
dx = eta/mPerPatch
heteroAmp = 0.9 % 0.9 is close to CDE's (4.1a)
%{
\end{matlab}

\paragraph{Loop over cell to heterogeneity ratios}
The micro-scale heterogeneity must be \(2\pi\)-periodic over
the macroscale domain, so in general has to be of the
following form for integer~\verb|nDetune| (positive means
smaller period, as shown by CDE, negative means larger).
\begin{matlab}
%}
RmsErrs = [];  etaEps = [];
for nDetune = -9:maxDetune
epsilon = 2*pi/(2*pi/eta+nDetune)
%{
\end{matlab}
\begin{figure}
\centering
\caption{\label{fig:Err}errors in equilibrium as function of
cell-size~\(\eta\) to heterogeneity
periodicity~\(\epsilon\): blue dots are \(r=1/9\); red dots
are \(r=1/3\).  The errors are \textsc{rms} of the
difference between solutions, such as in \cref{fig:ueq}, in
the five common patches/cells.  The errors are essentially
zero in the vicinity of the ideal integer ratio
for~\(\eta/\epsilon\).  The low error valleys appear broader
for larger ratio~\(r\).}
\includegraphics[width=\textwidth]{nearPeriodDiffEqui1Err}
\end{figure}
Consequently, the heterogeneity 
\begin{align*}
\cos(2\pi x/\epsilon) &
= \cos\left[ 2\pi x/\eta +2\pi(1/\epsilon-1/\eta)x\right]
\\&
= \cos(2\pi x/\eta)\cos(kx) -\sin(2\pi x/\eta)\sin(kx)
\\&\text{for wavenumber }k:=2\pi(1/\epsilon-1/\eta).
\end{align*}
So the discrepancy between cell-size and
heterogeneous-period can be viewed as a modulation of
precise cell-periodicity by `macroscale' modulation of
wavenumber~\(k\).   That is, the discrepancy may be viewed
as an example of a ``functionally graded material''. In the
patch scheme such modulations are resolved on patches of
spacing~\(H\) provided their wavenumber \(k<\pi/H\). That
is, patch-resolution requires
\begin{align*}&
2\pi(1/\epsilon-1/\eta)<\pi/H
\\&\iff \eta/\epsilon-1<(\eta/H)/2
\\&\iff \eta/\epsilon <1+r/2
\end{align*}
where \(r:=\eta/H\) is the patch ratio. For example, for
\(r=\tfrac13, \tfrac19\) we need \(\tfrac\eta\epsilon <
\tfrac76, \tfrac{19}{18} \approx 1.17, 1.06\) in order to
realise accuracy---see \cref{fig:Err}.

The sides of the valleys for \(r=1/3\) in \cref{fig:Err}
appear to be affected by higher-harmonic structures
developing in the equilibrium solution, and these structures
are then relatively poorly resolved. 







\subsection{Code to create the patch schemes}
\label{sec:sc2heterodiff}


Establish the global data struct~\verb|patches| for the
microscale heterogeneous lattice diffusion
system~\cref{eq:hetroDiff} solved on \(2\pi\)-periodic
domain.  Use spectral interpolation for best accuracy.
\begin{matlab}
%}
global patches
us=[];
for nPatch = nPatches
    nSubP = mPerPatch+2;
    configPatches1(@heteroDiff,Xlim,'periodic',nPatch ...
        ,0,dx,nSubP,'EdgyInt',true);
%{
\end{matlab}
Set the microscale heterogeneity with harmonic mean one, and
set the forcing so solution should have amplitude near one.
Choose forcing anti-symmetric to ensure solution effectively
satisfies BCs of \(u=0\) at \(x=0,\pi\).
\begin{matlab}
%}
  xMid = ( patches.x(1:end-1,:,:,:)+patches.x(2:end,:,:,:) )/2;
  patches.cs = 1./(1+ heteroAmp*cos(2*pi*xMid/epsilon) );
%  xMid=squeeze(xMid)
%  xMicro=squeeze(patches.x)
%  cMicro = squeeze(patches.cs)
  patches.f = 2*sin(patches.x); % sign or sin
  xs=squeeze(patches.x);
%{
\end{matlab}



\paragraph{Solve for equilibrium}
The linear system is of the form \(\vec f(\vec u) =J\vec u
+\vec f_0\)\,. Get the constant term in the system by
evaluating at zero.
\begin{matlab}
%}
  u0 = 0*patches.x;
  f0 = patchSys1(1,u0);
%{
\end{matlab}
Put NaNs on the patch-edges, to then find all micro-grid
points interior to the patches, and hence are the variables
in~\(\vec u\).
\begin{matlab}
%}
  u0([1 end],:,:,:)=nan; 
  i=find(~isnan(u0));
  nJac=length(i);
%{
\end{matlab}
Create Jacobian~\(J\) column by column: since linear we
numerically differentiate with unit vectors.   
\begin{matlab}
%}
  Jac=nan(nJac);
  for j=1:nJac
    u0(i)=((1:nJac)==j);
    dudt= (patchSys1(0,u0)-f0);
    Jac(:,j)=dudt(i);
  end
  assert(rcond(Jac)>1e-9,'Jacobian seems too ill-conditioned')
%{
\end{matlab}
Solve linear system.
\begin{matlab}
%}
  u0(i) = -sparse(Jac)\f0(i);
  ueq = squeeze(u0);
%{
\end{matlab}
Check the residual.
\begin{matlab}
%}
  res = patchSys1(1,u0);
  normRes = norm(res(i));
  assert(normRes<1e-8,"norm of the residual is too big")
%{
\end{matlab}
Plot one example equilibrium for comparison
\begin{matlab}
%}
    if nDetune==1
    figure(1)
    subplot(lMax,1,find(nPatch==nPatches))
    j=find(xs(2,:)>0);
    plot(xs(:,j),ueq(:,j))
    ylim([min(0,min(min(ueq(:,j)))) max(1.2,max(max(ueq(:,j))))]), 
    xlim([0 pi])
    ylabel("$u(x)$")
    if nPatch==nPatches(end), 
        xlabel("$x$"), drawnow
        %set(gca,'position',[.2 .2 r r])
        exportgraphics(gcf,mfilename+"ueq.pdf" ...
        ,'ContentType','vector')
        end%if
    end%if nDetune
%{
\end{matlab}
Use common patches for quantitative comparison
\begin{matlab}
%}
    if nPatch==nPatches(1), J=find(xs(2,:)>0);
    else J=3*J-1;
    end
    us = cat(3,us,ueq(:,J));
%{
\end{matlab}

End of the for-loop over the number of patches. Here compute
the errors in the patch solutions compared to the
full-domain exact solution.
\begin{matlab}
%}
end%for nPatch
uErrs = us-us(:,:,end);
rmsErrs = reshape( rms(uErrs,[1 2],'omitnan') ,1,[]);
log10rmsErrs = log10(rmsErrs)
RmsErrs = [RmsErrs; rmsErrs];
etaEps = [etaEps; eta/epsilon];
%{
\end{matlab}

At end of loop over detuning parameters, plot errors as
function of the cell to periodicity ratio.
\begin{matlab}
%}
end%for nDetune
figure(2)
semilogy(etaEps,RmsErrs,'.:')
xlabel("$\eta/\epsilon$")
ylabel("RMS error")
set(gca,'position',[.2 .2 .64 .64])
exportgraphics(gcf,mfilename+"Err.pdf" ...
    ,'ContentType','vector')
%{
\end{matlab}

End of the main script.



\subsection{\texttt{heteroDiff()}: heterogeneous diffusion}
\label{sec:heteroDiff}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 2D input arrays~\verb|u|
and~\verb|x| (via edge-value interpolation of
\verb|patchSys1|, \cref{sec:patchSys1}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
array of diffusivities~\(c_i\) have previously
been stored in struct~\verb|patches.cs|.
\begin{matlab}
%}
function ut = heteroDiff(t,u,patches)
  dx = diff(patches.x(2:3));   % space step
  i = 2:size(u,1)-1;   % interior points in a patch
  ut = nan+u;          % preallocate output array
  ut(i,:,:,:) = diff(patches.cs(:,1,:,:).*diff(u))/dx^2 ...
                -u(i,:,:,:) +patches.f(i,:,:,:); 
end% function
%{
\end{matlab}
\end{document}
%}
