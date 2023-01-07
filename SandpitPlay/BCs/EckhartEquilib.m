% Find an equilibrium in forced heterogeneous diffusion in
% 1D on patches as an example application of patches in
% space.  Use fsolve as quick to code, and also works for
% nonlinear, although slower to compute for this linear
% example.  Adapted from the second example of Eckhardt
% (2210.04536, sec 6.2.1).  Implement Dirichlet BCs using
% new facilities in the patch toolbox.   AJR, 8 Jan 2023
%!TEX root = doc.tex
%{
\section{\texttt{EckhartEquilib}: find an equilibrium of a 1D
heterogeneous diffusion via small patches}
\label{sec:EckhartEquilib}

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
patches.  By default this code is for \(\epsilon=0.004\)
where the microscale heterogeneity and patches are tiny.
}}
&
\def\extraAxisOptions{}
\raisebox{-\height}{\input{Figs/EckhardtEquilib}}
\end{tabular}
\end{figure}



\paragraph{First configure the patch system}
Establish the microscale heterogeneity has micro-period
\verb|mPeriod| on the lattice, and coefficients to match
Eckhardt2210.04536 \S6.2.1. 
\begin{matlab}
%}
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
rEpsilon = 250 % 25 for graphic, up to 2000 say
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
global patches
ordCC = 4
configPatches1(@heteroDiffF,[0 1],'chebyshev',nPatch ...
    ,ordCC,dx,nSubP,'EdgyInt',true,'hetCoeffs',a);
%{
\end{matlab}

Set the forcing coefficients, either the original parabolic,
or sinusoidal.  At time \(t=1\) the resultant forcing we
actually apply here is simply the sum of the two components.
\begin{matlab}
%}
if 0 % given forcing
  patches.f1 = 2*( patches.x-patches.x.^2 );
  patches.f2 = 2*0.5+0*patches.x;
else% simple sine forcing 
  patches.f1 = sin(pi*patches.x);
  patches.f2 = pi/2*sin(pi*patches.x);
end%if
%{
\end{matlab}


\paragraph{Find equilibrium with fsolve}
We seek the equilibrium for the forcing that applies at time
\(t=1\) (as if that specific forcing were applying for all
time).  Execute the function that invokes \verb|fsolve|. 
For this linear problem, it is computationally quicker using
a linear solver, but \verb|fsolve| is quicker in human time,
and generalises to nonlinear problems.   
\begin{matlab}
%}
u = squeeze(execFsolve)
%{
\end{matlab}
Then plot the equilibrium solution (\cref{fig:EckhardtEquilib}).
\begin{matlab}
%}
clf, plot(squeeze(patches.x),u,'.')
xlabel('space $x$'),ylabel('equilibrium $u(x)$')
%{
\end{matlab}
%Optionally write to graphic file.
%\begin{matlab}
%%}
%matlab2tikz('Figs/EckhardtEquilib.tex','showInfo',false ...
%,'noSize',true,'parseStrings',false,'showWarnings',false ...
%,'extraCode','\tikzsetnextfilename{Figs/EckhardtEquilib}' ...
%,'extraAxisOptions','\extraAxisOptions' )
%%{
%\end{matlab}

\paragraph{Code to execute fsolve}
We code the function \verb|execFsolve| to execute
\verb|fsolve| because easiest if a sub-function that
computes the time derivatives has access to variables
\verb|u0| and~\verb|i|.
\begin{matlab}
%}
function [u,normRes] = execFsolve
global patches
%{
\end{matlab}
Start the search from a zero field.
\begin{matlab}
%}
u0 = 0*patches.x;
%{
\end{matlab}
But set patch-edge values to \verb|Nan| in order to
use~\verb|i| to index the interior sub-patch points as they
are the variables.
\begin{matlab}
%}
u0([1 end],:,:,:) = nan;
i = find(~isnan(u0));
%{
\end{matlab}
Seek the equilibrium, and report the norm of the residual.
\begin{matlab}
%}
[u0(i),res] = fsolve(@duidt,u0(i));
normRes = norm(res)
%{
\end{matlab}

The aim is to zero the time derivatives \verb|duidt| in the
following function.  First, insert the vector of variables
into the patch-array of \verb|u0|.  Second, find the time
derivatives via the patch scheme, and finally return a
vector of those at the patch-internal points.
\begin{matlab}
%}
function res = duidt(ui)
  u = u0;   u(i) = ui;
  res = patchSys1(1,u);
  res = res(i);
end%function duidt
end%function execFsolve
%{
\end{matlab}

Fin.
%}