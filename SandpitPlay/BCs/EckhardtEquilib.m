% Find an equilibrium in forced heterogeneous diffusion in
% 1D on patches as an example application of patches in
% space.  Use fsolve as quick to code, and also works for
% nonlinear, although slower to compute for this linear
% example.  Adapted from the second example of Eckhardt
% (2210.04536, sec 6.2.1).  Implement Dirichlet BCs using
% new facilities in the patch toolbox.   AJR, 1 Feb 2023
%!TEX root = doc.tex
%{
\section{\texttt{EckhardtEquilib}: find an equilibrium of a 1D
heterogeneous diffusion via small patches}
\label{sec:EckhardtEquilib}

\cref{sec:Eckhardt2210eg2,sec:heteroDiffF} describe details
of the problem and more details of the following
configuration.  The aim is to find the equilibrium,
\cref{fig:EckhardtEquilib}, of the forced heterogeneous
system with a forcing corresponding to that applied at time
\(t=1\).  Computational efficiency comes from only computing
the microscale heterogeneity on small spatially sparse
patches, potentially much smaller than those shown in
\cref{fig:EckhardtEquilib}.
\begin{figure}
\centering
\begin{tabular}{@{}cc@{}}
\parbox[t]{10em}{\caption{\label{fig:EckhardtEquilib}%
Equilibrium of the heterogeneous diffusion problem with
forcing the same as that applied at time \(t=1\), and for
relatively large \(\epsilon=0.04\) so we can see the
patches.  By default this code sets \(\epsilon=0.004\)
whence the microscale heterogeneity and patches are tiny.
}} &
\def\extraAxisOptions{mark size=1pt}
\raisebox{-\height}{\input{Figs/EckhardtEquilib}}
\end{tabular}
\end{figure}



\paragraph{First configure the patch system}
Establish the microscale heterogeneity has micro-period
\verb|mPeriod| on the lattice, and coefficients to match
\cite{Eckhardt2022} [\S6.2.1]. 
\begin{matlab}
%}
clear all
global patches
%global OurCf2eps, OurCf2eps=true %option to save plots
mPeriod = 6
y = linspace(0,1,mPeriod+1)';
a = 1./(2-cos(2*pi*y(1:mPeriod)))
global microTimePeriod; microTimePeriod=0;
%{
\end{matlab}

Set the number of patches, the number of periods per patch,
and the spatial period~\(\epsilon\), via
integer~\(1/\epsilon\).
\begin{matlab}
%}
nPatch = 7
nPeriodsPatch = 1 % any integer
rEpsilon = 25 % 25 for graphic, up to 2000 say
dx = 1/(mPeriod*rEpsilon+1)
nSubP = nPeriodsPatch*mPeriod+2
%{
\end{matlab}

Establish the global data struct~\verb|patches| for the
microscale heterogeneous lattice diffusion system
\cref{eq:hetroDiffF} solved on domain~\([0,1]\), with
Chebyshev-like distribution of patches, and say fourth order
interpolation to provide the edge-values.  Use `edgy'
interpolation. 
\begin{matlab}
%}
ordCC = 4
configPatches1(@heteroDiffF,[0 1],'chebyshev',nPatch ...
    ,ordCC,dx,nSubP,'EdgyInt',true,'hetCoeffs',a);
%{
\end{matlab}

Set the forcing coefficients, either the original parabolic,
or exp-sinusoidal.  At time \(t=1\) the resultant forcing we
actually apply here is simply the sum of the two components.
\begin{matlab}
%}
if 0 % given forcing
  patches.f1 = 2*( patches.x-patches.x.^2 );
  patches.f2 = 2*0.5+0*patches.x;
else% simple exp-sine forcing 
  patches.f1 = sin(pi*patches.x).*exp(patches.x);
  patches.f2 = pi/2*sin(pi*patches.x).*exp(patches.x);
end%if
%{
\end{matlab}


\paragraph{Find equilibrium with fsolve} We seek the
equilibrium for the forcing that applies at time \(t=1\) (as
if that specific forcing were applying for all time). For
this linear problem, it is computationally quicker using a
linear solver, but \verb|fsolve| is quicker in human time,
Start the search from a zero field.
\begin{matlab}
%}
u = 0*patches.x;
%{
\end{matlab}
But set patch-edge values to \verb|Nan| in order to use
\verb|patches.i| to index the interior sub-patch points as
they are the variables.
\begin{matlab}
%}
u([1 end],:,:,:) = nan;
patches.i = find(~isnan(u));
%{
\end{matlab}
Seek the equilibrium, and report the norm of the residual,
via the generic patch system wrapper \verb|theRes|
(\cref{sec:theRes}).
\begin{matlab}
%}
[u(patches.i),res] = fsolve(@theRes,u(patches.i));
normRes = norm(res)
%{
\end{matlab}

\paragraph{Plot the equilibrium} see \cref{fig:EckhardtEquilib}.
\begin{matlab}
%}
clf, plot(squeeze(patches.x),squeeze(u),'.')
xlabel('space $x$'),ylabel('equilibrium $u(x)$')
ifOurCf2tex(mfilename)%optionally save
%{
\end{matlab}
%}