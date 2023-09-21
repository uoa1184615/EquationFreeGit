% Find an equilibrium in forced heterogeneous diffusion in
% 1D on patches as an example application of patches in
% space.  Use fsolve as quick to code, albeit slower to
% compute for this linear example.  Adapted from
% Hintermuller & Korolev (2309.04439, sec 3.4).   AJR, 21
% Feb 2023
%{
\section{\texttt{Hintermuller230904439Equilib}: find an
equilibrium of a 1D heterogeneous diffusion via small
patches}
\label{sec:Hintermuller230904439Equilib}

Find the equilibrium of a forced heterogeneous system. 
Computational efficiency comes from only computing the
microscale heterogeneity on small patches sparsely
distributed in space.


\paragraph{First configure the patch system}
Establish the microscale heterogeneity has micro-period
\verb|mPeriod| on the lattice, and coefficients to match
\cite{Eckhardt2022} [\S6.2.1]. 
\begin{matlab}
%}
clear all, close all
global patches
mPeriod = 6
y = linspace(0,1,mPeriod+1)';
K = 1./(1.2+sin(2*pi*y(1:mPeriod)))
%{
\end{matlab}
Specify order of coupling between patches: in this problem,
fourth-order is numerically exact!
\begin{matlab}
%}
ordCC = 2
%{
\end{matlab}
With second-order inter-patch interpolation the number of
iterations, function evaluations and errors are the
following for \(5,9,17,33,65\) patches:
\begin{verbatim}
fsolveIts =
     2     4     4     4     4
fsolveFuns =
    93        275       413       995        1565
rmsRelErrs =
    0.0656    0.0203    0.0035    0.0004         0
\end{verbatim}
For fourth-order inter-patch interpolation:
\begin{verbatim}
fsolveIts =
     2     4     4     4     4
fsolveFuns =
    93        275       515       995        1565
rmsRelErrs =
   1.0e-13 *
    0.6325    0.6574    0.2709    0.1279         0
\end{verbatim}
In contrast, H\&K (Table 2) with \(\epsilon=1/64\) achieve
relative errors 0.09--0.25 with 102,000--120,000 iterations.


Set the number of patches, the number of periods per patch,
and the spatial period~\(\epsilon\), via
integer~\(1/\epsilon\).
\begin{matlab}
%}
log2maxPatch = 6 % H&K up to 6 approx
rEpsilon = 1+2^log2maxPatch % H&K used 16, 32 or 64
dx = 1/(mPeriod*rEpsilon+1)
nSubP = mPeriod+2
%{
\end{matlab}

Loop over various number of patches for this~\(\epsilon\) up
to the full-domain solution which we take to be exact.  Put
solution into \verb|uCommon| for those patches common to
all, and~\verb|I| indexes which are common..
\begin{matlab}
%}
uCommon=[]; I=[];
fsolveIts=[]; fsolveFuns=[];
for l=2:log2maxPatch
  nPatch = 1+2^l
  if isempty(I), I=1:nPatch, else I=2*I-1, end
%{
\end{matlab}



Establish the global data struct~\verb|patches| for the
microscale heterogeneous lattice diffusion system
\cref{eq:hetroDiffF} solved on domain~\([0,1]\), with
say equispaced patches.  Use `edgy'
interpolation. 
\begin{matlab}
%}
configPatches1(@heteroDiffF,[0 1],'equispace',nPatch ...
    ,ordCC,dx,nSubP,'EdgyInt',true,'hetCoeffs',K);
%{
\end{matlab}

Set the forcing coefficient.
\begin{matlab}
%}
patches.f = -3*(patches.x-1);
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
[u(patches.i),res,flag,fout] = fsolve(@theRes,u(patches.i));
normRes = norm(res)
fsolveIts=[ fsolveIts fout.iterations ]
fsolveFuns=[ fsolveFuns fout.funcCount ]
%{
\end{matlab}

Plot the equilibrium 
\begin{matlab}
%}
clf, plot(squeeze(patches.x),squeeze(u),'.')
xlabel('space $x$'),ylabel('equilibrium $u(x)$')
drawnow,pause(1)
%ifOurCf2tex(mfilename)%optionally save
%{
\end{matlab}

Store common results, end-loop, and report errors.
\begin{matlab}
%}
uCommon = [uCommon reshape(u(2:end-1,:,:,I),[],1)];
end%for l
rmsRelErrs = rms( uCommon-uCommon(:,end) )./rms(uCommon)
%{
\end{matlab}





\subsection{\texttt{theRes()}: wrapper function to zero for equilibria}
\label{sec:theRes}
This functions converts a vector of values into the interior
values of the patches, then evaluates the time derivative of
the system at time \(t=1\), and returns the vector of
patch-interior time derivatives.
\begin{matlab}
%}
function f=theRes(u)
  global patches
  v=nan(size(patches.x));
  v(patches.i) = u;
  f = patchSys1(1,v(:),patches);
  f = f(patches.i);
end%function theRes
%{
\end{matlab}




\subsection{\texttt{heteroDiffF()}: forced heterogeneous diffusion}
\label{sec:heteroDiffF}

This function codes the lattice heterogeneous diffusion
inside the patches with forcing and with microscale boundary
conditions on the macroscale boundaries.  Computes the time
derivative at each point in the interior of a patch, output
in~\verb|ut|.  The column vector of diffusivities~\(K_i\)
has been stored in struct~\verb|patches.cs|, as has the
vector of forcing.
\begin{matlab}
%}
function ut = heteroDiffF(t,u,patches)
%{
\end{matlab}
Two basic parameters, and initialise result array to NaNs.
\begin{matlab}
%}
  dx = diff(patches.x(2:3));   % space step
  i = 2:size(u,1)-1;   % interior points in a patch
  ut = nan+u;          % preallocate output array
%{
\end{matlab}
The macroscale Dirichlet boundary conditions are zero at the
extreme edges of the two extreme patches.
\begin{matlab}
%}
  u( 1 ,:,:, 1 )=0; % left-edge of leftmost is zero
  u(end,:,:,end)=0; % right-edge of rightmost is zero
%{
\end{matlab}
Code the microscale forced diffusion.
\begin{matlab}
%}
  ut(i,:,:,:) = diff(patches.cs(:,1,:).*diff(u))/dx^2 ...
      +patches.f(i,:,:,:); 
end% function
%{
\end{matlab}

%}