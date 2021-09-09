% Simulate heterogeneous Landau--Lifshitz in 1D on patches
% with projective integration as an example application of
% patches in space.  From Leitenmaier & Runborg,
% http://arxiv.org/abs/2108.09463
% AJR, Sep 2021
%{
\section{\texttt{homoLanLif1D}: computational homogenisation
of a 1D heterogeneous Landau--Lifshitz by simulation on
small patches}
\label{sec:homoLanLif1D}
%\localtableofcontents

See the \verb|Patch| folder for an introduction.

Here try Projective Integration---which unfortunately is
`caught between a rock and a hard place'. The `rock' is that
the macroscale waves have non-small frequency, proportional
to the wavenumber squared: the macroscale step has to
resolve such frequencies for stablity. Whereas the `hard
place' is that the sub-patch microscale dissipation is also
proportional to wavenumber squared: the microscale burst has
to be long relative enough to dissipate such micro-waves.
For effective Projective integration there needs to be a big
enough gap between the burst length and the macro-step---so
parameter~\(\epsilon\) has to be genuinely small---the
\(\epsilon=1/200\) here is not small enough for Projective
Integration to be effective for small
dissipation~\(\alpha\).  



\paragraph{Parameters}
There are two closely related examples
\cite[pp.6,27]{Leitenmaier2021}, that we distinguish here
with parameter \verb|ex5p1|: set to either zero or one. The
Landau--Lifshitz dissipation parameter~\(\alpha\) should be
small.  However, need a bigger gap between the oscillations
of the macroscale and the decay of the microscale, so
set relatively large damping parameter~\(\alpha\).
\begin{matlab}
%}
clear all
global alpha ex5p1
ex5p1 = 0; % set to 1 for L&O example of p.27
alpha = 0.5 % phenomenological damping parameter 
%{
\end{matlab}
The physical microscale periodicity of the heterogeneity
is~\(\epsilon\) (\(\epsilon\)~is \emph{not} the patch scale
ratio):
\begin{matlab}
%}
epsilon = 1/200/(1+ex5p1) %pp.6,27
%{
\end{matlab}



\subsection{Script code to simulate heterogeneous diffusion systems}
\label{sec:sc2heterodiff}

First establish the microscale heterogeneity has
micro-period~\verb|mPeriod| on the lattice with values of
the column vector from \cite{Leitenmaier2021} [pp.6,27]. 
Later, the heterogeneity is repeated to fill each patch.
\begin{matlab}
%}
dx = 1/2000 %1/6000 %p.27
mPeriod = round(epsilon/dx)
a = 1 + 0.5*sin(2*pi*(0.5:mPeriod)'/mPeriod); %p.6
%{
\end{matlab}

Establish the global data struct~\verb|patches| for the
microscale heterogeneous lattice diffusion
system~\cref{eq:hetroDiff} solved on \(1\)-periodic domain,
with maybe 24~patches, but 11~is enough, here each patch of
size ratio to fit one period of the heterogeneity in each
patch, and spectral inter-patch interpolation to provide the
patch edge-values. Invoking \verb|EdgyInt| means the
edge-values come from interpolating the opposite
next-to-edge values of the patches (not the mid-patch
values).  
\begin{matlab}
%}
global patches
nPatch = 11 %24 %p.6, odd is slightly cleaner
nSubP = mPeriod+2 
ratio = nPatch*epsilon 
configPatches1(@heteroLanLif1D,[0 1],nan,nPatch ...
    ,0,ratio,nSubP,'EdgyInt',true ...
    ,'hetCoeffs',a);
assert(abs(dx-diff(patches.x(2:3)))<1e-10 ...
      ,'microscale grid spacing error')
%{
\end{matlab}


\paragraph{Simulate}
Set the initial conditions of a simulation to be that of
\cite{Leitenmaier2021} [pp.6], except possibly perturbed by
random microscale noise.  Scale the initial conditions so
that \(|\Mv(x,0)\|=1\)\,.
\begin{matlab}
%}
u0 = 0.5+exp(-0.1*cos(2*pi*(patches.x-0.32)));
v0 = 0.5+exp(-0.2*cos(2*pi*patches.x)) +0.03*randn(size(patches.x));
w0 = 0.5+exp(-0.1*cos(2*pi*(patches.x-0.75)));
M0 = [ u0 v0 w0 ]./sqrt(u0.^2+v0.^2+w0.^2); 
dM0dt = patchSys1(0,M0(:)); 
%{
\end{matlab}




\paragraph{Use projective integration in time}
Based upon \verb|homogenisationExample.m| of
\cref{sec:HomogenisationExample}. Now take
\verb|patchSys1|, the interface to the time derivatives,
and wrap around it the projective integration \verb|PIRK2|
(\cref{sec:PIRK2}), of bursts of simulation, as illustrated
by \cref{fig:HomogenisationU}.

This second part of the script implements the following
design, where the micro-integrator could be, for example,
\verb|rk2int| (as here \verb|ode23| is too microscale
rough).
\begin{enumerate} \def\itemsep{-1.5ex}
\item configPatches1 (done in first part)
\item PIRK2 \into heteroBurst \into rk2int \into
patchSys1 \into heteroLanLif1D
\item process results
\end{enumerate}
Mark that edge of patches are not to be used in the
projective extrapolation by setting initial values to \nan.
\begin{matlab}
%}
M0([1 end],:,:,:) = nan;
%{
\end{matlab}
Set the desired macro- and microscale time-steps over the
time domain: macroscale steps need to resolve the macroscale
oscillations of \cref{fig:homoLanLif1DSpec} with frequencies
up to \(7.4(N-1)^2\) so set to the reciprocal; the burst
time needs to complete the decay of the fast sub-patch
modes, rate~\(\alpha\cdot10^6\), so set the burst time to
thrice the reciprocal. 
\begin{matlab}
%}
dT = 1/(7.4*(nPatch-1)^2)
ts = 0:dT:0.1;
bT = 9e-6/alpha % 3e-6 is too small
addpath('../ProjInt')
microBurst = @(ti, ui, bT) ...
    rk2Int(@patchSys1,[ti ti+bT],ui(:));
[Ms,tss,Mss] = PIRK2(microBurst, ts, M0(:), bT);
%{
\end{matlab}


Reshape results for processing.  For simplicity, set edge
values to \verb|nan|s. For the field values (which are rows
in~\verb|Ms|) we need to reshape, permute,  and reshape
again. 
\begin{matlab}
%}
xs = squeeze(patches.x);  
Ms = reshape(Ms,length(ts),nSubP,3,nPatch);
Ms(:,[1 end],:,:) = nan; % nan patch edges
Ms = reshape( permute(Ms,[2 4 1 3]) ,[],length(ts),3);
%{
\end{matlab}
Check on constancy of~\(|\Mv(x,t)|\) in time.  The mean and
standard deviation appears to show that, with~\verb|ode15s|,
they are constant to errors typically~\(10^{-5}\).
\begin{matlab}
%}
Mabs = sqrt( sum(Ms.^2,3) );
meanMabs = mean(Mabs(:),'omitnan')
stdevMabs = std(Mabs(:),'omitnan')
%{
\end{matlab}



\paragraph{Plot space-time surface of the simulation}
Choose whether to save some plots, or not.
\begin{matlab}
%}
%%
global OurCf2eps
OurCf2eps = false;
%{
\end{matlab}

Subsampled surface over the macroscale duration of the
simulation to show the propagation of the macroscale
modes over the heterogeneous lattice.
\begin{matlab}
%}
figure(1),clf
if length(ts)>50
  [~,j]=min(abs(ts(:)-linspace(ts(1),ts(end),50)));
  else j=1:length(ts); end
uvw='uvw';
for p=1:3
  subplot(2,2,p)
  mesh(ts(j),xs(:),Ms(:,j,p))
  view(60,40), colormap(0.8*hsv)
  xlabel('time t'), ylabel('space x')
  zlabel([uvw(p) '(x,t)']) 
end
%{
\end{matlab}
Final time plot to compare with Fig.~2.1 of \cite{Leitenmaier2021}.
\begin{matlab}
%}
subplot(2,2,4)
plot(xs(:),squeeze(Ms(:,end,:)),'.')
xlabel('space x'), legend(uvw(1),uvw(2),uvw(3))
title(['time = ' num2str(ts(end),4)])
ifOurCf2eps([mfilename 'uvw'])
%{
\end{matlab}


Fin.
%}

